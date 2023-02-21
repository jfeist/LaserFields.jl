using SnoopPrecompile

@precompile_setup begin
    general_args = (is_vecpot=true,E0=1.5,ω0=0.12,t0=500.,chirp=0.,ϕ0=0.8π)
    test_fields = [
        GaussianLaserField(;      general_args...,σ=100.),
        SinExpLaserField(;        general_args...,T=800.,exponent=2),
        SinExpLaserField(;        general_args...,T=800.,exponent=4),
        SinExpLaserField(;        general_args...,T=800.,exponent=7),
        LinearFlatTopLaserField(; general_args...,Tflat=400.,Tramp=150),
        Linear2FlatTopLaserField(;general_args...,Tflat=400.,Tramp=150),
    ]

    @precompile_all_calls begin    
        ts = LinRange(0,1000,1001)
        for lf in test_fields
            lf.(ts)
            A_field.(lf,ts)
            try
                envelope.(lf,ts)
            catch
            end
        end        
        test_fields_fourier = filter(test_fields) do lf
            try
                E_fourier(lf,1.)
                return true
            catch
                return false
            end
        end
        
        ts = LinRange(0,10000,10001)
        ωs = LinRange(0,0.2,1001)[2:end]
        for lf in test_fields_fourier
            EF = E_fourier.(lf,ωs)
            AF = A_fourier.(lf,ωs)
            E = E_field.(lf,ts)
            A = A_field.(lf,ts)
        end
        
        lf = LaserField(form="gaussianI",is_vecpot=true,intensity_Wcm2=1e14,lambda_nm=45.,
                        peak_time_as=0,duration_as=1000.,ϕ0=0.3π,linear_chirp_rate_w0as=1e-4)
        ts = LinRange(-100,100,2001)
        A_field.(lf,ts)
        A_fourier.(lf,ωs)
        A_fourier.(lf,ωs)
        
        lf = make_laserfield(form="gaussianI", is_vecpot=true, phase_pi=1, duration_as=100.,
                            peak_time_as=400, intensity_Wcm2=1e14, lambda_nm=12., linear_chirp_rate_w0as=0.)
        lf = GaussianLaserField(; is_vecpot=true, ϕ0=π, E0=1., ω0=1., t0=0., σ=1., chirp=0.)
        ts = start_time(lf):TX(lf)/100:end_time(lf)
        lf.(ts)
    end
end