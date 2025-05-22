using LaserFields
using FFTW
using Test

@testset "LaserFields.jl" begin
    using LaserFields
    using LaserFields: GaussianLaserField, SinExpLaserField, LinearFlatTopLaserField, Linear2FlatTopLaserField, InterpolatingLaserField

    general_args = (is_vecpot=true, E0=1.5, ω0=0.12, t0=500., chirp=0., ϕ0=0.8π)
    test_fields = [
        GaussianLaserField(;      general_args..., σ=100.),
        SinExpLaserField(;        general_args..., T=800., exponent=2),
        SinExpLaserField(;        general_args..., T=800., exponent=4),
        SinExpLaserField(;        general_args..., T=800., exponent=7),
        LinearFlatTopLaserField(; general_args..., Tflat=400., Tramp=150),
        Linear2FlatTopLaserField(;general_args..., Tflat=400., Tramp=150),
    ]

    @testset "General arguments" for lf in test_fields
        @test lf.is_vecpot == true
        @test lf.E0 == 1.5
        @test lf.ω0 == 0.12
        @test lf.t0 == 500
        @test lf.chirp == 0.0
        @test lf.ϕ0 == 0.8π
    end

    @testset "LaserFieldCollection" begin
        lfc = LaserFieldCollection(test_fields)
        @test lfc isa LaserFieldCollection
        @test length(lfc.lfs) == 6
        @test lfc(500.) == sum(lf(500.) for lf in test_fields)
        @test E_field(lfc, 300.) == sum(E_field(lf, 300.) for lf in test_fields)
        @test A_field(lfc, 300.) == sum(A_field(lf, 300.) for lf in test_fields)
        @test E_fourier(lfc, 1.) == sum(E_fourier(lf, 1.) for lf in test_fields)
        @test A_fourier(lfc, 1.) == sum(A_fourier(lf, 1.) for lf in test_fields)
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
                lf = LaserField(form=form, is_vecpot=true, duration=1000., rampon=100., E0=1., omega=1., t0=0.)
                @test Teff(lf,n_photon) ≈ T
            end
        end
    end

    @testset "Fourier transform" begin
        # Test the Fourier transform of the laser fields
        # Compare the analytical and numerical Fourier transforms

        @testset "chirp $(chirp)" for chirp in (-1e-3, -1e-20, 0, 1e-20, 1e-3)
            general_args = (is_vecpot=true, E0=1.5, ω0=0.12, t0=500., chirp=chirp, ϕ0=0.8π)
            @testset "laserfield($lf)" for lf in [
                GaussianLaserField(;      general_args..., σ=100.),
                SinExpLaserField(;        general_args..., T=800., exponent=2),
                SinExpLaserField(;        general_args..., T=800., exponent=4),
                SinExpLaserField(;        general_args..., T=800., exponent=7),
                LinearFlatTopLaserField(; general_args..., Tflat=400., Tramp=150),
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

                atol = lf isa LinearFlatTopLaserField ? 0.02 : 1e-4
                @test all(isapprox.(Eω, Eω2, atol=atol))
            end
        end
    end
end
