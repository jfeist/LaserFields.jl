using LaserFields
using LaserFields: get_ωt
using Test

@testset "LaserFields.jl" begin
    # Write your tests here.
    using LaserFields
    using LaserFields: GaussianLaserField, SinExpLaserField, LinearFlatTopLaserField, Linear2FlatTopLaserField, InterpolatingLaserField

    general_args = (is_vecpot=true,E0=1.5,ω0=0.12,peak_time=500.,chirp=0.,phase_pi=0.8)
    test_fields = [
        GaussianLaserField(;     general_args...,duration=100.),
        SinExpLaserField(;       general_args...,duration=800.,exponent=2),
        SinExpLaserField(;       general_args...,duration=800.,exponent=4),
        SinExpLaserField(;       general_args...,duration=800.,exponent=7),
        LinearFlatTopLaserField(;general_args...,duration=400.,rampon=150),
        Linear2FlatTopLaserField(;general_args...,duration=400.,rampon=150),
        InterpolatingLaserField("laserdat.dat", is_vecpot=true),
        InterpolatingLaserField("laserdat.dat", is_vecpot=false)
    ];

    for lf in test_fields[begin:end-2]
        @test lf.is_vecpot == true
        @test lf.E0 == 1.5
        @test lf.ω0 == 0.12
        @test lf.peak_time == 500
        @test lf.chirp == 0.0
        @test lf.phase_pi == 0.8

        @test get_ωt(lf, 0.)   == 0.12
        @test get_ωt(lf, 500.) == 0.12
        @test get_ωt(lf, 600.) == 0.12
    end

    let lf = test_fields[end-1]
        @test lf.is_vecpot == true
        @test lf.E0 == 0.1598555396981788
        @test lf.ω0 == 0.16126036751849762
        @test lf.peak_time == 353.36430391128425
        @test lf.duration == 1000.0
        @test lf.phase_pi == 0.0
        @test lf.chirp == 0.0
        @test lf.datafile == "laserdat.dat"
        @test start_time(lf) == 0.0
        @test end_time(lf) == 1000.0
    end
    let lf = test_fields[end]
        @test lf.is_vecpot == false
        @test lf.E0 == 0.9968278450466621
        @test lf.ω0 == 0.16002937234660114
        @test lf.peak_time == 343.6735101653429
        @test lf.duration == 1000.0
        @test lf.phase_pi == 0.0
        @test lf.chirp == 0.0
        @test lf.datafile == "laserdat.dat"
        @test start_time(lf) == 0.0
        @test end_time(lf) == 1000.0
    end
end
