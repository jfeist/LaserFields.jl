abstract type LaserField end
Base.Broadcast.broadcastable(lf::LaserField) = Ref(lf)

TX(lf::LaserField) = 2π/lf.ω0

(lf::LaserField)(t) = E_field(lf,t)

envelope(lf::LaserField,t) = _envelope(lf,t-lf.t0)[1]

function E_field(lf::LaserField,t)
    tr = t - lf.t0
    env, envpr = _envelope(lf,tr)
    ϕt = lf.ϕ0 + lf.ω0*tr + lf.chirp*tr^2
    osc   = sin(ϕt)
    if lf.is_vecpot
        # dϕ/dt = ω0 + 2 chirp tr
        oscpr = (lf.ω0 + 2lf.chirp*tr)*cos(ϕt)
        return -(env * oscpr + envpr * osc) / lf.ω0
    else  # describes electric field directly
        return env * osc
    end
end

function A_field(lf::LaserField,t)
    lf.is_vecpot || error("laser field is not given as a vector potential, cannot get A(t) analytically!")

    tr = t - lf.t0
    env, _ = _envelope(lf,tr)
    osc = sin(lf.ϕ0 + lf.ω0*tr + lf.chirp*tr^2)
    # Divide out derivative of oscillation to ensure peak amplitude of E0 for electric field
    return env*osc / lf.ω0
end

"""Return the fourier transform of the envelope of the laser field.
We write the whole pulse as
f(t) = (env(t) exp(i*(phi0 + ω0*tr + chirp*tr^2)) + c.c. ) / 2im, where tr = t-tpeak
For the fourier transform of the envelope, we include the chirp term
exp(i chirp (t-tpeak)^2) in the envelope, so that its fourier transform is a complex function.
However, for unchirped pulses, the result will be purely real."""
function _envelope_fourier end

function E_fourier(lf::LaserField,ω)
    # analytically determine the fourier transform of the defined laser fields
    # determined as Int exp(-i*ω*t) E(t) dt

    # with tr = t-tpeak, the whole pulse is
    # f(t) =  env(t) sin    (phi0 + ω0*tr + chirp*tr^2)
    #      = (env(t) exp(IU*(phi0 + ω0*tr + chirp*tr^2)) - c.c. ) / 2im
    # for the fourier transform, we include the chirp term exp(i chirp tr^2) in the envelope.
    # this part is transformed in lf_envelope_fourier.
    # exp(IU*phi0) is just a constant prefactor, and the linear phase ω0*tr just gives a shift in frequency,
    # F[f(t) exp(im ω0 t)](ω) = F[f(t)](ω-ω0)
    # complex conjugation of the transformed function gives complex conjugation + reversal of the argument in the transform, so
    # F[conjg(f(t) exp(im ω0 t))](ω) = conjg(F[f(t) exp(IU ω0 t)](-ω)) = conjg(F[f(t)](-ω-ω0))
    ELFT = (   _envelope_fourier(lf, ω-lf.ω0) * cis(lf.ϕ0)
            - (_envelope_fourier(lf,-ω-lf.ω0) * cis(lf.ϕ0))') / 2im

    # the fourier transform of the part was determined as if it was centered around t=0
    # shift in time now -- just adds a phase exp(-im*ω*t0), as F[f(t-a)] = exp(-im*ω*a) F[f(t)]
    ELFT *= cis(-ω*lf.t0)

    if lf.is_vecpot
        # if this laser field was defined as a vector potential, we need to multiply
        # with -im*ω to get the fourier transform of the electric field:
        # E = -dA/dt --> F[-dA/dt] = -iω F[A]
        # in addition, we need to take into account that A0 = E0 / lf.ω0
        ELFT *= -1im * ω / lf.ω0
    end
    return ELFT
end

A_fourier(lf::LaserField,ω) = E_fourier(lf,ω) / (-1im*ω)
