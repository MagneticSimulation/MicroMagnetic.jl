using JuMag
using Test

mesh = FDMesh(nx = 2, ny = 1, nz = 1, dx = 1e-9, dy = 1e-9, dz = 1e-9)
sim = Sim(mesh, name = "test_ovf1")
set_Ms(sim, 8.0e5)
m_ovf = [0.6,0.8,0,0.6,0.8,0]
init_m0(sim, m_ovf)
save_ovf(sim, "test_ovf_bin", dataformat = "Binary 8")
save_ovf(sim, "test_ovf_text", dataformat = "Text")

function test_read_ovf_with_sim(ovf_name)
    sim = Sim(mesh, name = "test_ovf3")
    set_Ms(sim, 8.0e5)
    read_ovf(sim, ovf_name)
    println("spin:",sim.spin)
    for i =1:6
        @test sim.spin[i] == m_ovf[i]
    end
end

function test_read_ovf(ovf_name)
    ovf = read_ovf(ovf_name)
    @test length(ovf.data) == 6
    for i =1:6
        @test ovf.data[i] ==m_ovf[i]
    end
end

function test_save_ovf_without_sim(ovfname)
    ovf=read_ovf(ovfname)
    save_ovf(ovf,ovfname*"_binary")
    ovf.data_type="Text"
    save_ovf(ovf,ovfname*"text")
    ovf_bin=read_ovf(ovfname*"_binary")
    ovf_text=read_ovf(ovfname*"text")
    @test length(ovf_bin.data) == 6
    @test length(ovf_text.data) == 6
    for i =1:6
        @test ovf_bin.data[i] ==m_ovf[i]
        @test ovf_text.data[i] ==m_ovf[i]
    end
end

@testset "Test ovfs" begin
    test_read_ovf_with_sim("test_ovf_bin")
    test_read_ovf_with_sim("test_ovf_text")
    test_read_ovf("test_ovf_bin")
    test_read_ovf("test_ovf_text")
    test_save_ovf_without_sim("test_ovf_bin")
    test_save_ovf_without_sim("test_ovf_text")
end

#JuMag.cuda_using_double(true)
