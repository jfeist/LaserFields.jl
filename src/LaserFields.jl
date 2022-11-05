module LaserFields

using SpecialFunctions
using DelimitedFiles
using DataInterpolations

export make_laser_field, E_field, A_field, E_fourier, A_fourier, start_time, end_time

const GAUSSIAN_TIME_CUTOFF_SIGMA = 3.5*sqrt(log(256))
const au_as   = 1/24.188843265903884 # attosecond in a.u.
const au_wcm2toel2 = 1/3.5094455205784296e16 # W/cm^2 in a.u. for electric field squared
const au_wcm2 = 1.5536611487396207e-16 # W/cm^2 in a.u
const au_m    = 1/5.29177210903e-11  # m in a.u.
const au_cm   = 1/5.29177210903e-9   # cm in a.u.
const au_nm   = 1/5.29177210903e-2   # nm in a.u.
const au_c    = 137.0359990836958    # c (speed of light) in a.u. == 1/alpha
const au_eV   = 1/27.211386245935508 # eV in a.u.
const au_m_He = 7294.299386612553    # m of He nucleus in a.u.
const au_m_n  = 1838.6836617324586   # m of neutron in a.u.

abstract type LaserField end
Base.Broadcast.broadcastable(lf::LaserField) = Ref(lf)

macro _standard_lf_props()
    esc(quote
            is_vecpot::Bool
            E0::T1
            ω0::T2
            t0::T3
            ϕ0::T4
            chirp::T5
    end)
end

Base.@kwdef struct GaussianLaserField{T1,T2,T3,T4,T5,T6} <: LaserField
    @_standard_lf_props
    σ::T6
end

Base.@kwdef struct SinExpLaserField{T1,T2,T3,T4,T5,T6,T7} <: LaserField
    @_standard_lf_props
    T::T6
    exponent::T7
end

abstract type FlatTopLaserField <: LaserField end

Base.@kwdef struct LinearFlatTopLaserField{T1,T2,T3,T4,T5,T6,T7} <: FlatTopLaserField
    @_standard_lf_props
    Tflat::T6
    Tramp::T7
end

Base.@kwdef struct Linear2FlatTopLaserField{T1,T2,T3,T4,T5,T6,T7} <: FlatTopLaserField
    @_standard_lf_props
    Tflat::T6
    Tramp::T7
end

TX(lf::LaserField) = 2π/lf.ω0

(lf::LaserField)(t) = lf.is_vecpot ? A_field(lf,t) : E_field(lf,t)

function E_field(lf::LaserField,t)
    tr = t - lf.t0
    env, envpr = envelope(lf,tr)
    phit = lf.ϕ0 + lf.ω0*tr + lf.chirp*tr^2
    osc   = sin(phit)
    if lf.is_vecpot
        # d(w(t)*(t-peak))/dt = w(t) + w'(t)*(t-peak) = ω + lf.ω0 * lf.chirp * (t-tpeak) = 2ω - lf.ω0
        oscpr = (lf.ω0 + 2lf.chirp*tr)*cos(phit)
        return -(env * oscpr + envpr * osc) / lf.ω0
    else  # describes electric field directly
        return env * osc
    end
end

function A_field(lf::LaserField,t)
    lf.is_vecpot || error("laser field is not given as a vector potential, cannot get A(t) analytically!")

    tr = t - lf.t0
    env, envpr = envelope(lf,tr)
    osc = sin(lf.ϕ0 + lf.ω0*tr + lf.chirp*tr^2)
    # Divide out derivative of oscillation to ensure peak amplitude of E0 for electric field
    return env*osc / lf.ω0
end

"""return the fourier transform of the envelope of the laser field.
we write the whole pulse as
f(t) = (env(t) exp(i*(phi0 + ω0*tp + chirp*tp^2)) + c.c. ) / 2im, where tp = t-tpeak
for the fourier transform of the envelope, we include the chirp term
exp(i chirp (t-tpeak)^2) in the envelope, so that its fourier transform is a complex function.
however, for unchirped pulses, the result will be purely real!

for the various calculations, see chirped_fourier.nb in the mathematica directory."""
function envelope_fourier end

function E_fourier(lf::LaserField,ω)
    # analytically determine the fourier transform of the defined laser fields
    # determined as Int exp(-i*ω*t) E(t) dt

    # with tp = t-tpeak, the whole pulse is
    # f(t) =  env(t) sin    (phi0 + ω0*tp + chirp*tp^2)
    #      = (env(t) exp(IU*(phi0 + ω0*tp + chirp*tp^2)) - c.c. ) / 2im
    # for the fourier transform, we include the chirp term exp(i chirp tp^2) in the envelope.
    # this part is transformed in lf_envelope_fourier.
    # exp(IU*phi0) is just a constant prefactor, and the linear phase ω0*tp just gives a shift in frequency,
    # F[f(t) exp(im ω0 t)](ω) = F[f(t)](ω-ω0)
    # complex conjugation of the transformed function gives complex conjugation + reversal of the argument in the transform, so
    # F[conjg(f(t) exp(im ω0 t))](ω) = conjg(F[f(t) exp(IU ω0 t)](-ω)) = conjg(F[f(t)](-ω-ω0))
    ELFT = (   envelope_fourier(lf, ω-lf.ω0) * cis(lf.ϕ0)
            - (envelope_fourier(lf,-ω-lf.ω0) * cis(lf.ϕ0))') / 2im

    # the fourier transform of the part was determined as if it was centered around t=0
    # shift in time now -- just adds a phase exp(-im*ω*t0), as F[f(t-a)] = exp(-im*ω*a) F[f(t)]
    ELFT *= cis(-ω*lf.t0)

    if lf.is_vecpot
        # if this laser field was defined as a vector potential, we need to multiply
        # with -im*ω to get the fourier transform of the electric field, E=-dA/dt
        # F[-dA/dt] = -iω F[A]
        # in addition, we need to take into account that A0 = E0 / lf.ω0
        ELFT *= -1im * ω / lf.ω0
    end
    return ELFT
end

A_fourier(lf::LaserField,ω) = E_fourier(lf,ω) / (-1im*ω)

function envelope(lf::GaussianLaserField,tr)
    env   = lf.E0 * exp(-tr^2/(2*lf.σ^2))
    envpr = -env * tr/lf.σ^2
    return env,envpr
end
function envelope_fourier(lf::GaussianLaserField,ω)
    # F[exp(-z*t^2)] = exp(-w^2/4z)/sqrt(2z) (for real(z)>0)
    z = 0.5/lf.σ^2 - 1im*lf.chirp
    return lf.E0 * exp(-ω^2/4z) / sqrt(2z)
end
start_time(lf::GaussianLaserField) = lf.t0 - GAUSSIAN_TIME_CUTOFF_SIGMA*lf.σ
end_time(  lf::GaussianLaserField) = lf.t0 + GAUSSIAN_TIME_CUTOFF_SIGMA*lf.σ

function expiatbt2_intT(a,b,T)
    # returns the result of the integral Int(exp(i*(a*t+b*t**2)),{t,-T/2,T/2}) / sqrt(2π)
    zz1 = (1+1im)/4
    z34 = (-1+1im)/sqrt(2) # == (-1)**(3/4)
    cb = complex(b) # we want to take the square root and b might be negative
    res = erf(z34*(a-cb*T)/sqrt(4cb)) - erf(z34*(a+cb*T)/sqrt(4cb))
    res = res * zz1 / sqrt(cb) * cis(-a^2/4cb)
    # this is surprisingly not given by mathematica - not sure yet why it misses it,
    # but it's necessary for agreement with the numerical fourier transform
    return res * sign(b)
end

function envelope(lf::SinExpLaserField,tr)
    trel = tr/lf.T
    if abs(trel) > 0.5
        env   = 0.
        envpr = 0.
    else
        env   =  lf.E0 * cospi(trel)^lf.exponent
        envpr = -lf.E0 * sinpi(trel) * lf.exponent * cospi(trel)^(lf.exponent-1) * π/lf.T
    end
    return env,envpr
end

function envelope_fourier(lf::SinExpLaserField,ω)
    if lf.exponent == 2
        if lf.chirp == 0
            # the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
            return lf.E0 * sqrt(8π^3) * sinc(ω*lf.T/2π)/(8π^2/lf.T - 2*ω^2*lf.T)
        else
            # now we use that cos(pi*t/T)**2 * exp(i*c*t**2) can be written as 0.5 exp(i*c*t**2) + 0.25 exp(i*c*t**2 - 2*i*pi*t/T) + 0.25 exp(i*c*t**2 + 2*i*pi*t/T)
            # the integral of exp(IU*(a*t+b*t**2)) from t=-T/2 to t=T/2 can be calculated analytically and is implemented in the function below
            # the arguments are a={-ω, -2*pi/T-ω, 2*pi/T-ω} and b=chirp
            wd = 2π/lf.T
            return lf.E0 * (expiatbt2_intT(    - ω, lf.chirp, lf.T)/2 +
                            expiatbt2_intT(-wd - ω, lf.chirp, lf.T)/4 +
                            expiatbt2_intT( wd - ω, lf.chirp, lf.T)/4)
        end
    elseif lf.exponent == 4
        if lf.chirp == 0
            # the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
            return lf.E0 * 24 * (sqrt(2π^7) * sinc(ω*lf.T/2π) /
                                    (128π^4/lf.T - 40π^2*ω^2*lf.T + 2ω^4*lf.T^3))
        else
            # now we use that cos(pi*t/T)**4 * exp(i*c*t**2) can be written as
            # (0.375 exp(i*c*t**2) + 0.25 exp(i*c*t**2 - 2*i*pi*t/T) + 0.25 exp(i*c*t**2 + 2*i*pi*t/T) +
            #  0.0625 exp(i*c*t**2 - 4*i*pi*t/T) + 0.0625 exp(i*c*t**2 + 4*i*pi*t/T))
            # the integral of exp(IU*(a*t+b*t**2)) from t=-T/2 to t=T/2 can be calculated analytically and is implemented in the function below
            # the arguments are a={-ω, -2*pi/T-ω, 2*pi/T-ω, -4*pi/T-ω, 4*pi/T-ω} and b=chirp
            wd = 2π/lf.T
            return lf.E0 * (expiatbt2_intT(     - ω, lf.chirp, lf.T)*0.375  +
                            expiatbt2_intT( -wd - ω, lf.chirp, lf.T)*0.25   +
                            expiatbt2_intT(  wd - ω, lf.chirp, lf.T)*0.25   +
                            expiatbt2_intT(-2wd - ω, lf.chirp, lf.T)*0.0625 +
                            expiatbt2_intT( 2wd - ω, lf.chirp, lf.T)*0.0625)
        end
    else
        if lf.chirp != 0 || !isinteger(lf.exponent)
            error("sin_exp fourier transform with exponent != 2 or 4 only implemented for integer exponents and unchirped pulses")
        end
        x = 0.5*(ω*lf.T/π - lf.exponent)
        return lf.E0 * lf.T * gamma(lf.exponent+1)*gamma(x)*sinpi(x)/(sqrt(2^(2lf.exponent+1) * π^3) * gamma(x+lf.exponent+1))
    end
end

start_time(lf::SinExpLaserField) = lf.t0 - lf.T/2
end_time(  lf::SinExpLaserField) = lf.t0 + lf.T/2

start_time(lf::FlatTopLaserField)      = lf.t0 - lf.Tflat/2 - lf.Tramp
end_time(  lf::FlatTopLaserField)      = lf.t0 + lf.Tflat/2 + lf.Tramp

function envelope(lf::FlatTopLaserField,tr)
    # for linear field, the peak time is the middle of the interval
    if abs(tr) > lf.Tflat/2 + lf.Tramp
        env   = 0.
        envpr = 0.
    elseif abs(tr) > lf.Tflat/2
        trel  = (lf.Tramp + lf.Tflat/2 - abs(tr))/lf.Tramp
        env   = ramponfunc(lf,trel)
        envpr = -sign(tr)*ramponfuncpr(lf,trel) / lf.Tramp
    else
        env   = 1.
        envpr = 0.
    end
    return lf.E0 .* (env,envpr)
end

ramponfunc(lf::LinearFlatTopLaserField,trel) = trel
ramponfuncpr(lf::LinearFlatTopLaserField,trel) = 1.
ramponfunc(lf::Linear2FlatTopLaserField,trel) = sin(π/2*trel)^2
ramponfuncpr(lf::Linear2FlatTopLaserField,trel) = sin(π*trel) * π/2

function envelope_fourier(lf::LinearFlatTopLaserField,ω)
    lf.chirp == 0 || error("Fourier transform of 'linear' field with chirp not implemented!")
    return lf.E0 * sqrt(8/π) * sinc(ω*lf.Tramp/2π) * sinc(ω*(lf.Tramp+lf.Tflat)/2π) * (lf.Tramp+lf.Tflat)/4
end

function envelope_fourier(lf::Linear2FlatTopLaserField,ω)
    lf.chirp == 0 || error("Fourier transform of 'linear2' field with chirp not implemented!")
    return lf.E0 * sqrt(2π^3) * cos(ω*lf.Tramp/2) * sinc(ω*(lf.Tramp+lf.Tflat)/2π) * (lf.Tramp+lf.Tflat)/ (2π^2 - 2*lf.Tramp^2*ω^2)
end

Base.@kwdef struct InterpolatingLaserField{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: LaserField
    @_standard_lf_props
    duration::T6
    datafile::String
    Efun::T7
    Afun::T8
    start_time::T9
    end_time::T10
end

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

function make_laser_field(; form::String, is_vecpot::Bool, phase_pi=0, pargs...)
    args = values(pargs)
    if form == "readin"
        return InterpolatingLaserField(args.datafile; is_vecpot=is_vecpot)
    end

    E0 = if haskey(args,:E0)
        haskey(args,:intensity_Wcm2) && error("Cannot specify both E0 and intensity_Wcm2")
        args.E0
    else
        sqrt(args.intensity_Wcm2 * au_wcm2toel2)
    end

    omega = if haskey(args,:omega)
        haskey(args,:lambda_nm) && error("Cannot specify both omega and lambda_nm")
        args.omega
    else
        2π*au_c / (args.lambda_nm * au_nm)
    end

    chirp = if haskey(args,:chirp)
        haskey(args,:linear_chirp_rate_w0as) && error("Cannot specify both chirp and linear_chirp_rate_w0as")
        args.chirp
    elseif haskey(args,:linear_chirp_rate_w0as)
        omega * args.linear_chirp_rate_w0as / au_as
    else
        0
    end

    peak_time = if haskey(args,:peak_time)
        haskey(args,:peak_time_as) && error("Cannot specify both peak_time and peak_time_as")
        args.peak_time
    else
        args.peak_time_as * au_as
    end

    duration = if haskey(args,:duration)
        haskey(args,:duration_as) && error("Cannot specify both duration and duration_as")
        args.duration
    else
        args.duration_as * au_as
    end

    rampon = if haskey(args,:rampon)
        haskey(args,:rampon_as) && error("Cannot specify both rampon and rampon_as")
        args.rampon
    elseif haskey(args,:rampon_as)
        args.rampon_as * au_as
    else
        0
    end
    kwargs = Dict(pairs((is_vecpot=is_vecpot, ϕ0=π*phase_pi, E0=E0, ω0=omega,
                         t0=peak_time, chirp=chirp)))
    if form in ("gaussian","gaussianF")
        # convert from FWHM of field to standard deviation of field
        kwargs[:σ] = duration / sqrt(log(256))
        return GaussianLaserField(; kwargs...)
    elseif form in ("gaussian2","gaussianI")
        # convert from FWHM of intensity to standard deviation of field
        kwargs[:σ] = duration / sqrt(log(16))
        return GaussianLaserField(; kwargs...)
    elseif form in ("sin2","sin4","sin_exp")
        kwargs[:T] = duration
        kwargs[:exponent] = form=="sin2" ? 2 : (form=="sin4" ? 4 : args.form_exponent)
        return SinExpLaserField(; kwargs...)
    elseif form in ("linear","linear2")
        kwargs[:Tflat] = duration
        kwargs[:Tramp] = rampon
        lftype = form=="linear" ? LinearFlatTopLaserField : Linear2FlatTopLaserField
        return lftype(; kwargs...)
    else
        error("Unknown laser field form '$form'")
    end
end

end