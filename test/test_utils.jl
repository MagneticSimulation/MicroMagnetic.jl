using Printf

function test_functions(test_name, funcs...; platforms=["CPU", "CUDA", "AMDGPU", "oneAPI", "Metal"], precisions=[Float32, Float64])
    for platform in platforms
        if Base.find_package(platform) === nothing && platform != "CPU"
            continue
        end
        
        # The runner deliberately manages the global backend/precision between
        # groups; zero the Sim counter around the switches so that the
        # "existing Sim(s)" warning does not fire (same trick as in
        # test_global_state.jl).
        saved_n_sims = MicroMagnetic._n_sims[]
        MicroMagnetic._n_sims[] = 0
        if !set_backend(platform)
            MicroMagnetic._n_sims[] = saved_n_sims
            continue
        end

        for precision in precisions
            if platform == "Metal" && precision == Float64
                continue
            end

            name = @sprintf("%s %s %s", test_name, platform, precision)
            @testset "$name" begin

                set_precision(precision)
                for func in funcs
                    func()
                end
            end
        end
        MicroMagnetic._n_sims[] = saved_n_sims
    end
end