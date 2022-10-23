module LaserFields

using SpecialFunctions
using DelimitedFiles
using DataInterpolations

const GAUSSIAN_TIME_CUTOFF_SIGMA = 3.5*sqrt(log(256))

abstract type LaserField end
Base.Broadcast.broadcastable(lf::LaserField) = Ref(lf)

macro _standard_props()
    esc(quote
            is_vecpot::Bool
            E0::T1
            omega::T2
            peak_time::T3
            duration::T4
            chirp::T5
            phase_pi::T6
    end)
end

Base.@kwdef struct GaussianLaserField{T1,T2,T3,T4,T5,T6} <: LaserField
    @_standard_props
end

Base.@kwdef struct SinExpLaserField{T1,T2,T3,T4,T5,T6,T7} <: LaserField
    @_standard_props
    exponent::T7
end

abstract type FlatTopLaserField <: LaserField end

Base.@kwdef struct LinearFlatTopLaserField{T1,T2,T3,T4,T5,T6,T7} <: FlatTopLaserField
    @_standard_props
    rampon::T7
end

Base.@kwdef struct Linear2FlatTopLaserField{T1,T2,T3,T4,T5,T6,T7} <: FlatTopLaserField
    @_standard_props
    rampon::T7
end

TX(lf::LaserField) = 2π/lf.omega
get_omega(lf::LaserField,t) = lf.omega + lf.chirp*(t-lf.peak_time)

function E(lf::LaserField,t)
    env, envpr = envelope(lf,t)
    omega = get_omega(lf,t)
    osc   = sin(omega*(t-lf.peak_time) + π*lf.phase_pi)
    if lf.is_vecpot
        # d(w(t)*(t-peak))/dt = w(t) + w'(t)*(t-peak) = omega + lf%omega * lf%chirp * (t-peak) = 2 * omega - lf%omega
        oscpr = (2*omega-lf.omega)*cos(omega*(t-lf.peak_time)+π*lf.phase_pi)
        return -(env * oscpr + envpr * osc) / lf.omega
    else  # describes electric field directly
        return env * osc
    end
end

function A(lf::LaserField,t)
    lf.is_vecpot || throw(ArgumentError("laser field is not given as a vector potential, cannot get A(t) analytically!"))

    env, envpr = envelope(lf,t)
    omega = get_omega(lf,t)
    osc = sin(omega * (t-lf.peak_time) + π*lf.phase_pi)
    # Divide out derivative of oscillation to ensure peak amplitude of E0 for electric field
    return env*osc / lf.omega
end

"""return the fourier transform of the envelope of the laser field.
we write the whole pulse as
f(t) = (env(t) exp(i*(phi0 + w0*tp + chirp*tp**2)) + c.c. ) / (2*IU), where tp = t-tpeak
for the fourier transform of the envelope, we include the chirp term
exp(i chirp (t-tpeak)**2) in the envelope, so that its fourier transform is a complex function.
however, for unchirped pulses, the result will be purely real!

for the various calculations, see chirped_fourier.nb in the mathematica directory."""
function envelope_fourier end

function E_fourier(lf::LaserField,omega)
    # analytically determine the fourier transform of the defined laser fields
    # determined as Int exp(-i*omega*t) E(t) dt

    # with tp = t-tpeak, the whole pulse is
    # f(t) =  env(t) sin    (phi0 + w0*tp + chirp*tp**2)
    #      = (env(t) exp(IU*(phi0 + w0*tp + chirp*tp**2)) - c.c. ) / (2*IU)
    # for the fourier transform, we include the chirp term exp(i chirp tp**2) in the envelope.
    # this part is transformed in lf_envelope_fourier.
    # exp(IU*phi0) is just a constant prefactor, and the linear phase w0*tp just gives a shift in frequency,
    # F[f(t) exp(IU w0 t)](w) = F[f(t)](w-w0)
    # complex conjugation of the transformed function gives complex conjugation + reversal of the argument in the transform, so
    # F[conjg(f(t) exp(IU w0 t))](w) = conjg(F[f(t) exp(IU w0 t)](-w)) = conjg(F[f(t)](-w-w0))
    ELFT = (   envelope_fourier(lf, omega - lf.omega) * cispi(lf.phase_pi)
            - (envelope_fourier(lf,-omega - lf.omega) * cispi(lf.phase_pi))') / 2im

    # the fourier transform of the part was determined as if it was centered around t=0
    # shift in time now -- just adds a phase exp(-IU*omega*peak_time), as F[f(t-a)] = exp(-IU*omega*a) F[f(t)]
    ELFT *= cis(-omega*lf.peak_time)

    if lf.is_vecpot
        # if this laser field was defined as a vector potential, we need to multiply
        # with -IU*omega to get the fourier transform of the electric field, E=-dA/dt
        # F[-dA/dt] = -iw F[A]
        # in addition, we need to take into account that A0 = E0 / lf%omega
        ELFT *= -1im * omega / lf.omega
    end
    return ELFT
end

A_fourier(lf::LaserField,omega) = E_fourier(lf,omega) / (-1im*omega)

function envelope(lf::GaussianLaserField,t)
    env   = lf.E0 * exp(-(t-lf.peak_time)^2/(2*lf.duration^2))
    envpr = env * (lf.peak_time-t)/lf.duration^2
    return env,envpr
end
function envelope_fourier(lf::GaussianLaserField,omega)
    # F[exp(-z*t^2)] = exp(-w^2/4z)/sqrt(2z) (for real(z)>0)
    z = 0.5/lf.duration^2 - 1im*lf.chirp
    return lf.E0 * exp(-omega^2/4z) / sqrt(2z)
end
start_time(lf::GaussianLaserField) = lf.peak_time - GAUSSIAN_TIME_CUTOFF_SIGMA*lf.duration
end_time(  lf::GaussianLaserField) = lf.peak_time + GAUSSIAN_TIME_CUTOFF_SIGMA*lf.duration


function expiatbt2_intT(a,b,T)
    # returns the result of the integral Int(exp(i*(a*t+b*t**2)),{t,-T/2,T/2}) / sqrt(2*pi)
    zz1 = (1+1im)/4
    z34 = (-1+1im)/sqrt(2) # == (-1)**(3/4)
    cb = complex(b) # we want to take the square root and b might be negative
    res = erf(z34*(a-cb*T)/sqrt(4cb)) - erf(z34*(a+cb*T)/sqrt(4cb))
    res = res * zz1 / sqrt(cb) * cis(-a^2/4cb)
    # this is surprisingly not given by mathematica - not sure yet why it misses it,
    # but it's necessary for agreement with the numerical fourier transform
    return res * sign(b)
end

function envelope(lf::SinExpLaserField,t)
    trel = (t-lf.peak_time)/lf.duration
    if abs(trel) > 0.5
        env   = 0.
        envpr = 0.
    else
        env   =  lf.E0 * cospi(trel)^lf.exponent
        envpr = -lf.E0 * sinpi(trel) * lf.exponent * cospi(trel)^(lf.exponent-1) * π/lf.duration
    end
    return env,envpr
end

function envelope_fourier(lf::SinExpLaserField,omega)
    if lf.exponent == 2
        if lf.chirp == 0
            # the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
            return lf.E0 * sqrt(8π^3) * sinc(omega*lf.duration/2π)/(8π^2/lf.duration - 2*omega^2*lf.duration)
        else
            # now we use that cos(pi*t/T)**2 * exp(i*c*t**2) can be written as 0.5 exp(i*c*t**2) + 0.25 exp(i*c*t**2 - 2*i*pi*t/T) + 0.25 exp(i*c*t**2 + 2*i*pi*t/T)
            # the integral of exp(IU*(a*t+b*t**2)) from t=-T/2 to t=T/2 can be calculated analytically and is implemented in the function below
            # the arguments are a={-omega, -2*pi/T-omega, 2*pi/T-omega} and b=chirp
            wd = 2π/lf.duration
            return lf.E0 * (expiatbt2_intT(    - omega, lf.chirp, lf.duration)/2 +
                            expiatbt2_intT(-wd - omega, lf.chirp, lf.duration)/4 +
                            expiatbt2_intT( wd - omega, lf.chirp, lf.duration)/4)
        end
    elseif lf.exponent == 4
        if lf.chirp == 0
            # the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
            return lf.E0 * 24 * (sqrt(2π^7) * sinc(omega*lf.duration/2π) /
                                    (128π^4/lf.duration - 40π^2*omega^2*lf.duration + 2omega^4*lf.duration^3))
        else
            # now we use that cos(pi*t/T)**4 * exp(i*c*t**2) can be written as
            # (0.375 exp(i*c*t**2) + 0.25 exp(i*c*t**2 - 2*i*pi*t/T) + 0.25 exp(i*c*t**2 + 2*i*pi*t/T) +
            #  0.0625 exp(i*c*t**2 - 4*i*pi*t/T) + 0.0625 exp(i*c*t**2 + 4*i*pi*t/T))
            # the integral of exp(IU*(a*t+b*t**2)) from t=-T/2 to t=T/2 can be calculated analytically and is implemented in the function below
            # the arguments are a={-omega, -2*pi/T-omega, 2*pi/T-omega, -4*pi/T-omega, 4*pi/T-omega} and b=chirp
            wd = 2π/lf.duration
            return lf.E0 * (expiatbt2_intT(      - omega, lf.chirp, lf.duration)*0.375  +
                            expiatbt2_intT(  -wd - omega, lf.chirp, lf.duration)*0.25   +
                            expiatbt2_intT(   wd - omega, lf.chirp, lf.duration)*0.25   +
                            expiatbt2_intT(-2*wd - omega, lf.chirp, lf.duration)*0.0625 +
                            expiatbt2_intT( 2*wd - omega, lf.chirp, lf.duration)*0.0625)
        end
    else
        if lf.chirp != 0 || !isinteger(lf.exponent)
            error("sin_exp fourier transform with exponent != 2 or 4 only implemented for integer exponents and unchirped pulses")
        end
        x = 0.5*(omega*lf.duration/π - lf.exponent)
        return lf.E0 * lf.duration * gamma(lf.exponent+1)*gamma(x)*sinpi(x)/(sqrt(2^(2lf.exponent+1) * π^3) * gamma(x+lf.exponent+1))
    end
end

start_time(lf::SinExpLaserField) = lf.peak_time - lf.duration/2
end_time(  lf::SinExpLaserField) = lf.peak_time + lf.duration/2

start_time(lf::FlatTopLaserField)      = lf.peak_time - lf.duration/2 - lf.rampon
end_time(  lf::FlatTopLaserField)      = lf.peak_time + lf.duration/2 + lf.rampon
flat_start_time(lf::FlatTopLaserField) = lf.peak_time - lf.duration/2
flat_end_time(  lf::FlatTopLaserField) = lf.peak_time + lf.duration/2

function envelope(lf::FlatTopLaserField,t)
    # for linear field, the peak time is the middle of the interval
    if t < start_time(lf)
        env   = 0.
        envpr = 0.
    elseif t < flat_start_time(lf)
        trel  = (t-start_time(lf))/lf.rampon
        env   = ramponfunc(lf,trel)
        envpr = ramponfuncpr(lf,trel) / lf.rampon
    elseif t < flat_end_time(lf)
        env   = 1.
        envpr = 0.
    elseif t < end_time(lf)
        trel = (end_time(lf) - t)/lf.rampon
        env   = ramponfunc(lf,trel)
        envpr = -ramponfuncpr(lf,trel) / lf.rampon
    else
        env   = 0.
        envpr = 0.
    end
    return lf.E0 .* (env,envpr)
end

ramponfunc(lf::LinearFlatTopLaserField,trel) = trel
ramponfuncpr(lf::LinearFlatTopLaserField,trel) = 1.
ramponfunc(lf::Linear2FlatTopLaserField,trel) = sin(π/2*trel)^2
ramponfuncpr(lf::Linear2FlatTopLaserField,trel) = sin(π*trel) * π/2

function envelope_fourier(lf::LinearFlatTopLaserField,omega)
    lf.chirp == 0 || error("Fourier transform of 'linear' field with chirp not implemented!")
    return lf.E0 * sqrt(8/π) * sinc(omega*lf.rampon/2π) * sinc(omega*(lf.rampon+lf.duration)/2π) * (lf.rampon+lf.duration)/4
end

function envelope_fourier(lf::Linear2FlatTopLaserField,omega)
    lf.chirp == 0 || error("Fourier transform of 'linear2' field with chirp not implemented!")
    return lf.E0 * sqrt(2π^3) * cos(omega*lf.rampon/2) * sinc(omega*(lf.rampon+lf.duration)/2π) * (lf.rampon+lf.duration)/ (2π^2 - 2*lf.rampon^2*omega^2)
end

struct InterpolatingLaserField{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10} <: LaserField
    @_standard_props
    datafile::String
    Efun::T7
    Afun::T8
    start_time::T9
    end_time::T10
end

function InterpolatingLaserField(; is_vecpot,datafile)
    # print('# Reading laser_field from file:', datafile)
    data = readdlm(datafile)
    ndims(data) == 2 && size(data,2) == 2 || error("Laser field datafile '$datafile' must contain two columns: time and field")
    tt, ff = eachcol(data)

    # print('# Number of data points found:', len(tt))
    issorted(tt) || error("Laser field datafile '$datafile' must be sorted by time")

    if is_vecpot
        Afun = CubicSpline(ff,tt)
        Efun = t -> -DataInterpolations.derivative(Afun,t)
    else
        Efun = CubicSpline(ff,tt)
        Afun = t -> -DataInterpolations.integral(Efun,tt[1],t)
    end

    # analyze the data we have read to guess some information about the field
    start_time = tt[1]
    end_time = tt[end]

    TX = Inf
    E0 = 0.
    Eprev = 0.
    peak_time = 0.
    lastzerocrossing = start_time
    for t in LinRange(tt[1],tt[end],20*length(tt))
        E = Efun(t)
        if abs(E) > E0
            E0 = abs(E)
            peak_time = t
        end
        if E*Eprev < 0
            TX = min(TX,2*(t-lastzerocrossing))
            lastzerocrossing = t
        end
    end

    omega = 2π / TX
    duration = tt[end] - tt[1]
    chirp = 0.
    phase_pi = 0.

    InterpolatingLaserField(is_vecpot,E0,omega,peak_time,duration,chirp,phase_pi,datafile,Efun,Afun,start_time,end_time)
end

start_time(lf::InterpolatingLaserField) = lf.start_time
end_time(  lf::InterpolatingLaserField) = lf.end_time

E(lf::InterpolatingLaserField,t) = (start_time(lf) <= t <= end_time(lf)) ? lf.Efun(t) : 0.
A(lf::InterpolatingLaserField,t) = (start_time(lf) <= t <= end_time(lf)) ? lf.Afun(t) : 0.

end