macro _laserfield_struct(Name,args...)
    nfreefields = 5
    defs = Any[]
    for arg in args
        if arg isa Symbol
            nfreefields += 1
            push!(defs, :($(arg)::$(Symbol(:T,nfreefields))))
        elseif arg isa Expr && arg.head == :(::)
            push!(defs, arg)
        else
            throw(ArgumentError("invalid argument to new_lf_field"))
        end
    end

    DerivT = :LaserField
    if Name isa Expr && Name.head == :(<:)
        DerivT = Name.args[2]
        Name = Name.args[1]
    end

    Ts = Symbol.(:T,1:nfreefields)
    esc(quote
        Base.@kwdef struct $(Name){$(Ts...)} <: $(DerivT)
            is_vecpot::Bool
            E0::T1
            ω0::T2
            t0::T3
            ϕ0::T4
            chirp::T5
            $(defs...)
        end
    end)
end

# laser field with a Gaussian envelope with std dev σ
@_laserfield_struct GaussianLaserField σ

_envelope(lf::GaussianLaserField,tr) = (env = lf.E0 * exp(-tr^2/(2*lf.σ^2)); (env, -env * tr/lf.σ^2))
# F[exp(-z*t^2)] = exp(-w^2/4z)/sqrt(2z) (for real(z)>0)
_envelope_fourier(lf::GaussianLaserField,ω) = (z = 0.5/lf.σ^2 - 1im*lf.chirp; lf.E0 * exp(-ω^2/4z) / sqrt(2z))
start_time(lf::GaussianLaserField) = lf.t0 - GAUSSIAN_TIME_CUTOFF_SIGMA*lf.σ
end_time(  lf::GaussianLaserField) = lf.t0 + GAUSSIAN_TIME_CUTOFF_SIGMA*lf.σ
Teff(lf::GaussianLaserField,n_photon) = lf.σ * sqrt(π/n_photon)

@_laserfield_struct SinExpLaserField T exponent


"returns the result of the integral Int(exp(i*(a*t+b*t**2)),{t,-T/2,T/2}) / sqrt(2π)"
function expiatbt2_intT(a,b,T)
    bTsq = b * T^2
    if abs(bTsq) <= 1e-5
        # use first-order expansion for small b to avoid numerical errors
        aT = a * T
        if abs(aT) < 1e-8
            # avoid numerical errors for small aT
            return (1 + 1im/12 * bTsq - bTsq^2/160) * T / sqrt(2π)
        end
        x1 = 2im * b / a^2
        x2 = x1 * cos(aT / 2) + (1 - x1 + 0.25im * bTsq) * sinc(aT / 2π)
        return x2 * T / sqrt(2π)
    end
    t1 = inv(sqrt(complex(b))) # b might be negative
    zz1 = (1+1im)/4 * t1 # = (1+1im)/(2*sqrt(4b))
    z34 = (-1+1im)/sqrt(8) * t1 # == (-1)**(3/4) / sqrt(4b)
    (erf(z34*(a-b*T)) - erf(z34*(a+b*T))) * zz1 * cis(-a^2/4b)
end

function _envelope(lf::SinExpLaserField,tr)
    trel = tr/lf.T
    if abs(trel) > 0.5
        (0., 0.)
    else
        lf.E0 .* (cospi(trel)^lf.exponent,
                  -lf.exponent*π/lf.T * cospi(trel)^(lf.exponent-1) * sinpi(trel))
    end
end

function _envelope_fourier(lf::SinExpLaserField,ω)
    isinteger(lf.exponent) || error("sin_exp fourier transform only implemented for integer exponents")
    # rewrite the _envelope as a sum of exponentials, which are easy to Fourier transform over a limited time interval
    # cos(πt/T)^n = 1/2^n (exp(iπt/T)+ exp(-iπt/T))^n = 1/2^n sum_k=0^n binomial(n,k) exp(i(n-2k)πt/T)
    n = Int(lf.exponent)
    wd = π/lf.T
    res = 0im
    for k = 0:n
        res += binomial(n,k) * expiatbt2_intT((n-2k)*wd - ω, lf.chirp, lf.T)
    end
    return lf.E0/2^n * res
end

start_time(lf::SinExpLaserField) = lf.t0 - lf.T/2
end_time(  lf::SinExpLaserField) = lf.t0 + lf.T/2
Teff(lf::SinExpLaserField,n_photon) = lf.T * gamma(0.5 + n_photon*lf.exponent) / (sqrt(π)*gamma(1 + n_photon*lf.exponent))

abstract type FlatTopLaserField <: LaserField end

@_laserfield_struct LinearFlatTopLaserField<:FlatTopLaserField Tflat Tramp
@_laserfield_struct Linear2FlatTopLaserField<:FlatTopLaserField Tflat Tramp

start_time(lf::FlatTopLaserField) = lf.t0 - lf.Tflat/2 - lf.Tramp
end_time(  lf::FlatTopLaserField) = lf.t0 + lf.Tflat/2 + lf.Tramp

function _envelope(lf::FlatTopLaserField,tr)
    # for linear field, the peak time is the middle of the interval
    if abs(tr) > lf.Tflat/2 + lf.Tramp
        (0., 0.)
    elseif abs(tr) > lf.Tflat/2
        trel  = (lf.Tramp + lf.Tflat/2 - abs(tr))/lf.Tramp
        lf.E0 .* (ramponfunc(lf,trel), -sign(tr)*ramponfuncpr(lf,trel) / lf.Tramp)
    else
        (lf.E0, 0.)
    end
end

ramponfunc(::LinearFlatTopLaserField,trel) = trel
ramponfuncpr(::LinearFlatTopLaserField,trel) = 1.
ramponfunc(::Linear2FlatTopLaserField,trel) = sin(π/2*trel)^2
ramponfuncpr(::Linear2FlatTopLaserField,trel) = sin(π*trel) * π/2

function _envelope_fourier(lf::LinearFlatTopLaserField,ω)
    lf.chirp == 0 || error("Fourier transform of 'linear' field with chirp not implemented!")
    return lf.E0 * sqrt(8/π) * sinc(ω*lf.Tramp/2π) * sinc(ω*(lf.Tramp+lf.Tflat)/2π) * (lf.Tramp+lf.Tflat)/4
end

function _envelope_fourier(lf::Linear2FlatTopLaserField,ω)
    lf.chirp == 0 || error("Fourier transform of 'linear2' field with chirp not implemented!")
    return lf.E0 * sqrt(2π^3) * cos(ω*lf.Tramp/2) * sinc(ω*(lf.Tramp+lf.Tflat)/2π) * (lf.Tramp+lf.Tflat)/ (2π^2 - 2*lf.Tramp^2*ω^2)
end

Teff(lf::LinearFlatTopLaserField,n_photon) = lf.Tflat + 2*lf.Tramp / (1+2*n_photon)
Teff(lf::Linear2FlatTopLaserField,n_photon) = lf.Tflat + 2*lf.Tramp * gamma(0.5+2n_photon) / (sqrt(π)*gamma(1+2n_photon))

@_laserfield_struct InterpolatingLaserField duration datafile::String Efun Afun start_time end_time

function InterpolatingLaserField(datafile; is_vecpot)
    # print('# Reading laser_field from file:', datafile)
    data = readdlm(datafile)
    ndims(data) == 2 && size(data,2) == 2 || error("Laser field datafile '$datafile' must contain two columns: time and field")
    tt, ff = eachcol(data)

    # print('# Number of data points found:', len(tt))
    issorted(tt) || error("Laser field datafile '$datafile' must be sorted by time")

    # analyze the data we have read to guess some information about the field
    start_time = tt[1]
    end_time = tt[end]
    
    if is_vecpot
        Afun = CubicSpline(ff,tt)
        Efun = t -> -DataInterpolations.derivative(Afun,t)
    else
        Efun = CubicSpline(ff,tt)
        Afun = t -> -DataInterpolations.integral(Efun,start_time,t)
    end

    # guess the parameters of the field - note that this is just a simple estimation, not anything rigorous
    TX = Inf
    E0 = 0.
    Eprev = Efun(tt[1])
    t0 = 0.
    lastzerocrossing = -Inf
    for t in LinRange(tt[1],tt[end],20*length(tt))
        E = Efun(t)
        if abs(E) > E0
            E0 = abs(E)
            t0 = t
        end
        if sign(E) != sign(Eprev)
            TX = min(TX,2*(t-lastzerocrossing))
            lastzerocrossing = t
            # @show t, lastzerocrossing, TX
            Eprev = E
        end
    end

    ω0 = 2π / TX
    duration = tt[end] - tt[1]
    chirp = 0.
    ϕ0 = 0.

    InterpolatingLaserField(; is_vecpot,E0,ω0,t0,duration,chirp,ϕ0,datafile,Efun,Afun,start_time,end_time)
end

start_time(lf::InterpolatingLaserField) = lf.start_time
end_time(  lf::InterpolatingLaserField) = lf.end_time

E_field(lf::InterpolatingLaserField,t) = (start_time(lf) <= t <= end_time(lf)) ? lf.Efun(t) : 0.
A_field(lf::InterpolatingLaserField,t) = (start_time(lf) <= t <= end_time(lf)) ? lf.Afun(t) : 0.
