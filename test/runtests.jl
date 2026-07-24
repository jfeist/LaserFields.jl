using LaserFields
using FFTW
using Test

@testset "LaserFields.jl" begin
    using LaserFields
    using LaserFields: GaussianLaserField, SinExpLaserField, LinearFlatTopLaserField, Linear2FlatTopLaserField, InterpolatingLaserField

    general_args = (is_vecpot=true, E0=1.5, ω0=0.12, t0=500., chirp=0., ϕ0=0.8π)
    general_args_nonvec = (is_vecpot=false, E0=1.5, ω0=0.12, t0=500., chirp=0., ϕ0=0.8π)
    test_fields = [
        GaussianLaserField(;      general_args..., σ=100.),
        SinExpLaserField(;        general_args..., T=800., exponent=2),
        SinExpLaserField(;        general_args..., T=800., exponent=4),
        SinExpLaserField(;        general_args..., T=800., exponent=7),
        LinearFlatTopLaserField(; general_args_nonvec..., Tflat=400., Tramp=150),
        Linear2FlatTopLaserField(;general_args..., Tflat=400., Tramp=150),
    ]

    @testset "General arguments" begin
        expected_is_vecpot = [true, true, true, true, false, true]
        for (i, lf) in enumerate(test_fields)
            @test lf.is_vecpot == expected_is_vecpot[i]
        @test lf.E0 == 1.5
        @test lf.ω0 == 0.12
        @test lf.t0 == 500
        @test lf.chirp == 0.0
        @test lf.ϕ0 == 0.8π
        end
    end

    @testset "Time-domain field values" begin
        sample_t_shifts = (-0.2, 0.0, 0.3)
        expected_E_values = [
            [-0.4733721103340174, 1.213525491562421, 0.43939711942241655],
            [-0.46657740126715147, 1.213525491562421, 0.4560190763122897],
            [-0.4696176702845257, 1.213525491562421, 0.44856301902567197],
            [-0.47415631450844775, 1.213525491562421, 0.4374727545356417],
            [1.4265847744427302, 0.8816778784387098, -1.4265847744427274],
            [-0.4635254915624215, 1.213525491562421, 0.4635254915624302],
        ]
        expected_A_values = [
            [11.823200448244128, 7.347315653655915, -11.742442578919206],
            [11.868113281037923, 7.347315653655915, -11.843028666243946],
            [11.848054069403913, 7.347315653655915, -11.798022564282457],
            [11.818028803324419, 7.347315653655915, -11.730833895083071],
            [0.0, 0.0, 0.0],
            [11.888206453689419, 7.347315653655915, -11.888206453689396],
        ]
        expected_E_posfreq_values = [
            [ComplexF64(-0.2366860551670087, 0.7073805747087495), ComplexF64(0.6067627457812105, 0.4408389392193549), ComplexF64(0.21969855971120827, -0.7075431243057029)],
            [ComplexF64(-0.23328870063357574, 0.7114637065144673), ComplexF64(0.6067627457812105, 0.4408389392193549), ComplexF64(0.22800953815614486, -0.7115150383005721)],
            [ComplexF64(-0.23480883514226286, 0.7096391697345191), ComplexF64(0.6067627457812105, 0.4408389392193549), ComplexF64(0.22428150951283599, -0.7097408968808412)],
            [ComplexF64(-0.23707815725422388, 0.7069101152176173), ComplexF64(0.6067627457812105, 0.4408389392193549), ComplexF64(0.21873637726782086, -0.7070857016210783)],
            [ComplexF64(0.7132923872213651, 0.23176274578121075), ComplexF64(0.4408389392193549, -0.6067627457812105), ComplexF64(-0.7132923872213637, -0.23176274578121514)],
            [ComplexF64(-0.23176274578121075, 0.7132923872213651), ComplexF64(0.6067627457812105, 0.4408389392193549), ComplexF64(0.2317627457812151, -0.7132923872213638)],
        ]
        expected_A_posfreq_values = [
            [ComplexF64(5.911600224122064, 1.9207953490721235), ComplexF64(3.6736578268279576, -5.0563562148434205), ComplexF64(-5.871221289459603, -1.907675437887428)],
            [ComplexF64(5.934056640518961, 1.9280918810662835), ComplexF64(3.6736578268279576, -5.0563562148434205), ComplexF64(-5.921514333121973, -1.9240166383568338)],
            [ComplexF64(5.924027034701957, 1.924833064590886), ComplexF64(3.6736578268279576, -5.0563562148434205), ComplexF64(-5.899011282141228, -1.9167049538678569)],
            [ComplexF64(5.909014401662209, 1.9199551644239556), ComplexF64(3.6736578268279576, -5.0563562148434205), ComplexF64(-5.865416947541536, -1.9057894928745776)],
            [ComplexF64(0.0, 0.0), ComplexF64(0.0, 0.0), ComplexF64(0.0, 0.0)],
            [ComplexF64(5.944103226844709, 1.9313562148434231), ComplexF64(3.6736578268279576, -5.0563562148434205), ComplexF64(-5.944103226844698, -1.9313562148434595)],
        ]

        for (i, lf) in enumerate(test_fields)
            TX = LaserFields.TX(lf)
            sample_times = [lf.t0 + shift * TX for shift in sample_t_shifts]
            for (j, t) in enumerate(sample_times)
                @test E_field(lf, t) ≈ expected_E_values[i][j] atol=1e-12
                @test E_posfreq(lf, t) ≈ expected_E_posfreq_values[i][j] atol=1e-12
                @test E_field(lf, t) ≈ E_posfreq(lf, t) + conj(E_posfreq(lf, t)) atol=1e-12
                if lf.is_vecpot
                    @test A_field(lf, t) ≈ expected_A_values[i][j] atol=1e-12
                    @test A_posfreq(lf, t) ≈ expected_A_posfreq_values[i][j] atol=1e-12
                    @test A_field(lf, t) ≈ A_posfreq(lf, t) + conj(A_posfreq(lf, t)) atol=1e-12
                end
            end

            t_before = start_time(lf) - 0.1 * TX
            t_after = end_time(lf) + 0.1 * TX
            if lf isa GaussianLaserField
                @test abs(E_field(lf, t_before)) < 1e-12
                @test abs(E_field(lf, t_after)) < 1e-12
                @test abs(A_field(lf, t_before)) < 1e-12
                @test abs(A_field(lf, t_after)) < 1e-12
            else
                @test E_field(lf, t_before) == 0
                @test E_field(lf, t_after) == 0
                if lf.is_vecpot
                    @test A_field(lf, t_before) == 0
                    @test A_field(lf, t_after) == 0
                end
            end
        end
    end

    @testset "LaserFieldCollection" begin
        lfc = LaserFieldCollection(test_fields)
        lfc_vecpot = LaserFieldCollection((test_fields[1], test_fields[2], test_fields[3], test_fields[4], test_fields[6]))
        sample_t_shifts = (-0.2, 0.0, 0.3)
        @test lfc isa LaserFieldCollection
        @test length(lfc.lfs) == 6
        @test lfc(500.) == sum(lf(500.) for lf in test_fields)
        @test E_field(lfc, 300.) == sum(E_field(lf, 300.) for lf in test_fields)
        @test E_posfreq(lfc, 300.) == sum(E_posfreq(lf, 300.) for lf in test_fields)
        for shift in sample_t_shifts
            t = lfc.lfs[1].t0 + shift * LaserFields.TX(lfc)
            @test E_field(lfc, t) ≈ E_posfreq(lfc, t) + conj(E_posfreq(lfc, t)) atol=1e-12
        end
        @test A_field(lfc_vecpot, 300.) == sum(A_field(lf, 300.) for lf in (test_fields[1], test_fields[2], test_fields[3], test_fields[4], test_fields[6]))
        @test A_posfreq(lfc_vecpot, 300.) == sum(A_posfreq(lf, 300.) for lf in (test_fields[1], test_fields[2], test_fields[3], test_fields[4], test_fields[6]))
        for shift in sample_t_shifts
            t = lfc_vecpot.lfs[1].t0 + shift * LaserFields.TX(lfc_vecpot)
            @test A_field(lfc_vecpot, t) ≈ A_posfreq(lfc_vecpot, t) + conj(A_posfreq(lfc_vecpot, t)) atol=1e-12
        end
        @test E_fourier(lfc, 1.) == sum(E_fourier(lf, 1.) for lf in test_fields)
        @test A_fourier(lfc_vecpot, 1.) == sum(A_fourier(lf, 1.) for lf in (test_fields[1], test_fields[2], test_fields[3], test_fields[4], test_fields[6]))
        @test start_time(lfc) == minimum(start_time, test_fields)
        @test end_time(lfc) == maximum(end_time, test_fields)
    end

    @testset "read-in field vecpot" begin
        lf = InterpolatingLaserField("laserdat.dat", is_vecpot=true)
        @test lf.is_vecpot == true
        @test lf.E0 == 0.15985607526339093
        @test lf.ω0 == 0.16132596121126513
        @test lf.t0 == 353.37042585063125
        @test lf.duration == 700.0
        @test lf.ϕ0 == 0.0
        @test lf.chirp == 0.0
        @test lf.datafile == "laserdat.dat"
        @test start_time(lf) == 0.0
        @test end_time(lf) == 700.0
        @test_throws ErrorException E_posfreq(lf, 350.0)
        @test_throws ErrorException A_posfreq(lf, 350.0)
    end

    @testset "read-in field e-field" begin
        lf = InterpolatingLaserField("laserdat.dat", is_vecpot=false)
        @test lf.is_vecpot == false
        @test lf.E0 == 0.996830886761803
        @test lf.ω0 == 0.16009446532415655
        @test lf.t0 == 343.63364005991866
        @test lf.duration == 700.0
        @test lf.ϕ0 == 0.0
        @test lf.chirp == 0.0
        @test lf.datafile == "laserdat.dat"
        @test start_time(lf) == 0.0
        @test end_time(lf) == 700.0
        @test_throws ErrorException E_posfreq(lf, 350.0)
        @test_throws ErrorException A_posfreq(lf, 350.0)
    end

    @testset "LaserField" begin
        @testset "working" begin
            lf = LaserField(form="gaussianI", is_vecpot=true, phase_pi=1, duration_as=100.,
                            peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test lf isa GaussianLaserField
            @test lf.is_vecpot == true
            @test lf.σ == 100. * LaserFields.au_as / sqrt(log(16.))
            @test lf.t0 == 400. * LaserFields.au_as
            @test lf(lf.t0) == lf.E0
            @test lf.ϕ0 ≈ π
        end
        @testset "continuity constraints" begin
            @test_throws ErrorException LaserField(form="linear", is_vecpot=true, phase_pi=0.0, duration_as=100., rampon_as=20.,
                                                   peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test_throws ErrorException LaserField(form="sin_exp", is_vecpot=true, phase_pi=0.0, duration_as=100., form_exponent=1.0,
                                                   peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
        end
        @testset "overspecified parameters" begin
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration=10., duration_as=100.,
                                                   peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration_as=100.,
                                                   peak_time=0., peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration_as=100.,
                                                   peak_time_as=400, E0=0.3, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
        end
        @testset "missing parameters" begin
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5,
                                                   peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration_as=100.,
                                                   intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
            @test_throws ErrorException LaserField(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration_as=100.,
                                                   peak_time_as=400, lambda_nm=12., linear_chirp_rate_w0as=0.)
        end
    end

    @testset "Teff" begin
        refTs = Dict("gaussianI"  => [      1064.4670194312,      752.69184778925,       614.5703202121,      532.23350971561,      476.04412305096 ],
                     "gaussianF"  => [      752.69184778925,      532.23350971561,      434.56684093796,      376.34592389463,      336.61402755334 ],
                     "sin2"       => [                  375,             273.4375,          225.5859375,      196.38061523438,      176.19705200195 ],
                     "sin4"       => [             273.4375,      196.38061523438,      161.18025779724,      139.94993409142,      125.37068761958 ],
                     "linear"     => [      1066.6666666667,                 1040,      1028.5714285714,      1022.2222222222,      1018.1818181818 ],
                     "linear2"    => [                 1075,            1054.6875,         1045.1171875,      1039.2761230469,      1035.2394104004 ])
        for (form, Teffs) in refTs
            for (n_photon, T) in enumerate(Teffs)
                is_vecpot = form == "linear" ? false : true
                lf = LaserField(form=form, is_vecpot=is_vecpot, duration=1000., rampon=100., E0=1., omega=1., t0=0.)
                @test Teff(lf,n_photon) ≈ T
            end
        end
    end

    @testset "Fourier transform" begin
        # Test the Fourier transform of the laser fields
        # Compare the analytical and numerical Fourier transforms
        @testset "chirp $(chirp)" for chirp in (0.0011, -0.0009, -1e-3, -1e-20, 0, 1e-20, 1e-3, 0.0009, 0.0011)
            general_args = (is_vecpot=true, E0=1.5, ω0=0.12, t0=500., chirp=chirp, ϕ0=0.8π)
            @testset "laserfield($lf)" for lf in [
                GaussianLaserField(;      general_args..., σ=100.),
                SinExpLaserField(;        general_args..., T=100., exponent=2),
                SinExpLaserField(;        general_args..., T=100., exponent=4),
                SinExpLaserField(;        general_args..., T=100., exponent=7),
                LinearFlatTopLaserField(; merge(general_args, (is_vecpot=false,))..., Tflat=400., Tramp=150),
                Linear2FlatTopLaserField(;general_args..., Tflat=400., Tramp=150),
                ]
                if lf isa Union{LinearFlatTopLaserField,Linear2FlatTopLaserField} && chirp != 0
                    # Skip the test for LinearFlatTopLaserField with non-zero chirp
                    # because the analytical Fourier transform is not implemented for this case
                    continue
                end
                T = end_time(lf) - start_time(lf)
                t0 = start_time(lf) - 5T
                t1 = end_time(lf) + 5T
                dt = LaserFields.TX(lf) / 100
                ts = t0:dt:t1
                ωs = 2π * fftfreq(length(ts), 1/dt)
                Eω = E_fourier.(lf,ωs)
                Eω2 = fft(lf.(ts)) .* dt ./ sqrt(2π)
                # FFT acts as if ts[1] was t=0, shift to the correct value
                @. Eω2 *= exp(-1im * ts[1] * ωs)

                atol = lf isa LinearFlatTopLaserField ? 0.02 : 1e-3
                @test all(isapprox.(Eω, Eω2, atol=atol))
            end
        end
    end
end
