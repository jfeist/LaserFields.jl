using LaserFields
using Test

@testset "LaserFields.jl" begin
    using LaserFields
    using LaserFields: GaussianLaserField, SinExpLaserField, LinearFlatTopLaserField, Linear2FlatTopLaserField, InterpolatingLaserField

    general_args = (is_vecpot=true, E0=1.5, ω0=0.12, t0=500., chirp=0., ϕ0=0.8π)
    test_fields = [
        GaussianLaserField(;      general_args..., duration=100.),
        SinExpLaserField(;        general_args..., duration=800., exponent=2),
        SinExpLaserField(;        general_args..., duration=800., exponent=4),
        SinExpLaserField(;        general_args..., duration=800., exponent=7),
        LinearFlatTopLaserField(; general_args..., duration=400., rampon=150),
        Linear2FlatTopLaserField(;general_args..., duration=400., rampon=150),
    ]

    for lf in test_fields
        @test lf.is_vecpot == true
        @test lf.E0 == 1.5
        @test lf.ω0 == 0.12
        @test lf.t0 == 500
        @test lf.chirp == 0.0
        @test lf.ϕ0 == 0.8π
    end

    let lf = InterpolatingLaserField("laserdat.dat", is_vecpot=true)
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

    let lf = InterpolatingLaserField("laserdat.dat", is_vecpot=false)
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

    let lf = make_laser_field(form="gaussianI", is_vecpot=true, phase_pi=0.5, duration_as=100.,
                              peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
        @test lf isa GaussianLaserField
        @test lf.is_vecpot == true
        @test lf.duration == 100. * LaserFields.au_as / sqrt(log(16.))
        @test lf.t0 == 400. * LaserFields.au_as
        @test lf(lf.t0) == lf.E0/lf.ω0
        @test lf.ϕ0 == 0.5π
    end
end
