using Printf

function test_functions(test_name, funcs...; platforms=["CPU", "CUDA", "AMDGPU", "oneAPI", "Metal"], precisions=[Float32, Float64])
    for platform in platforms
        if Base.find_package(platform) === nothing && platform != "CPU"
            continue
        end
        
        if !set_backend(platform)
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
    end
end