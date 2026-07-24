abstract type LaserField end
Base.Broadcast.broadcastable(lf::LaserField) = Ref(lf)

TX(lf::LaserField) = 2π/lf.ω0

(lf::LaserField)(t) = E_field(lf,t)

envelope(lf::LaserField,t) = _envelope(lf,t-lf.t0)[1]

# one "generic" field evaluation function in order to avoid unnecessary code duplication for A or E
# fields, and for positive-frequency or full fields. The actual field evaluation functions just call this.
# make the boolean arguments compile-time constants, so that the compiler can optimize away the branches.
function _field_impl(lf::LaserField,t,::Val{A},::Val{posfreq}) where {A,posfreq}
    A && !lf.is_vecpot && error("laser field is not given as a vector potential, cannot get A(t) analytically!")

    tr = t - lf.t0
    env, envpr = _envelope(lf,tr)
    ϕt = lf.ϕ0 + lf.ω0*tr + lf.chirp*tr^2
    osc = posfreq ? 0.5im * cis(-ϕt) : sin(ϕt)

    if A
        return env*osc / lf.ω0
    elseif lf.is_vecpot
        # dϕ/dt = ω0 + 2 chirp tr
        oscpr = (lf.ω0 + 2lf.chirp*tr) * (posfreq ? -1im*osc : cos(ϕt))
        return -(env * oscpr + envpr * osc) / lf.ω0
    else  # describes electric field directly
        return env * osc
    end
end

E_field(lf::LaserField,t)   = _field_impl(lf,t,Val(false),Val(false))
A_field(lf::LaserField,t)   = _field_impl(lf,t,Val(true),Val(false))
E_posfreq(lf::LaserField,t) = _field_impl(lf,t,Val(false),Val(true))
A_posfreq(lf::LaserField,t) = _field_impl(lf,t,Val(true),Val(true))

"""Return the "positive-frequency part" of the laser field, i.e. the part of the field that
oscillates as exp(-i ω0 t) with ω0>0.

NOTE: With the convention for the sign of the Fourier transform used here, the "positive-frequency
part" (which is the part responsible for photon absorption) actually corresponds to *negative*
frequency arguments in the Fourier transform.

Note also that the "positive" frequency part as defined here in any case formally has a non-zero
Fourier transform on the whole frequency axis. Since we write E(t) = F(t)
exp(-i ω0 t) + c.c., and we take the first part as the positive-frequency part. The Fourier
transform of this part will have negative frequencies as well since the envelopes are formally
nonzero at all frequencies, but the peak of the spectrum will be at ω0>0.
"""
E_posfreq, A_posfreq


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

struct LaserFieldCollection{T} <: LaserField
    lfs::T
end
start_time(lf::LaserFieldCollection) = minimum(start_time, lf.lfs)
end_time(lf::LaserFieldCollection) = maximum(end_time, lf.lfs)
TX(lf::LaserFieldCollection) = minimum(TX, lf.lfs)
for f in (:E_field, :A_field, :E_fourier, :A_fourier, :E_posfreq, :A_posfreq, :envelope)
    @eval $f(lf::LaserFieldCollection, t) = sum(Base.Fix2($f,t), lf.lfs)
end
