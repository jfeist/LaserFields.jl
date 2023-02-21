
"""Get a value from one of several possibly defined parameters, or a default value if none are defined.
Makes sure that at most one of the options are specified, and that at least one is given if there is no default value.

    sample use: `x = @select_param args (x => args.x, x_squared => sqrt(x_squared), 0.)`
    
This makes it possible to call a function with either `x`, `x_squared`, or neither, and in the last case, `x` is set to 0.
"""
macro select_param(nt,options)
    @assert options.head ∈ (:block,:tuple)
    params = []
    default = nothing
    for arg in options.args
        if arg isa Expr && arg.head == :call
            ex = arg.args
            @assert ex[1] == :(=>) && length(ex) == 3
            push!(params, ex[2] => ex[3])
        elseif !(arg isa LineNumberNode)
            @assert isnothing(default) "Cannot have more than one default argument"
            default = arg
        end
    end
    parnamestr = join(first.(params),", ")
    checkexpr(head,par) = Expr(head, :(haskey($nt,$(QuoteNode(first(par))))), last(par))
    ifexpr = currif = checkexpr(:if, params[1])
    for par in params[2:end]
        push!(currif.args, checkexpr(:elseif, par))
        currif = currif.args[end]
    end
    if isnothing(default)
        push!(currif.args, :(error("ERROR: You need to specify one out of: ($($parnamestr))!")))
    else
        push!(currif.args, default)
    end
    esc(quote
        if sum(haskey.(Ref($nt),$(first.(params)))) > 1
            error("ERROR: Cannot specify more than one out of: ($($parnamestr))\npassed arguments: $($nt)")
        end
        $ifexpr
    end)
end

# General function to make a laser field with the parameter conventions from fortran laserfields library
make_laserfield(args...; kwargs...) = LaserField(args...; kwargs...)

LaserField(d) = LaserField(; d...)

function LaserField(; form::String, is_vecpot::Bool, pargs...)
    args = values(pargs)
    if form == "readin"
        return InterpolatingLaserField(args.datafile; is_vecpot=is_vecpot)
    end

    E0 = @select_param args (E0 => args.E0, intensity_Wcm2 => sqrt(args.intensity_Wcm2 * au_wcm2toel2))
    ω0 = @select_param args (ω0 => args.ω0, omega => args.omega, lambda_nm => 2π*au_c / (args.lambda_nm * au_nm))
    ϕ0 = @select_param args (ϕ0 => args.ϕ0, phase_pi => π*args.phase_pi, 0.)
    chirp = @select_param args begin
        chirp => args.chirp
        linear_chirp_rate_w0as => ω0 * args.linear_chirp_rate_w0as / au_as
        0.
    end
    t0 = @select_param args (t0 => args.t0, peak_time => args.peak_time, peak_time_as => args.peak_time_as * au_as)
    duration = @select_param args (duration => args.duration, duration_as => args.duration_as * au_as)
    Tramp = @select_param args (Tramp => args.Tramp, rampon => args.rampon, rampon_as => args.rampon_as * au_as, 0.)

    kwargs = Dict(pairs((; is_vecpot, ϕ0, E0, ω0, t0, chirp)))
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
        kwargs[:Tramp] = Tramp
        lftype = form=="linear" ? LinearFlatTopLaserField : Linear2FlatTopLaserField
        return lftype(; kwargs...)
    else
        error("Unknown laser field form '$form'")
    end
end
