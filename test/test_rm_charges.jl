using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function relax_rm_charges(;nx=30, ny=8, nz=3)
    I = argmax((nx, ny, nz))
    mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=2e-9, dy=2e-9, dz=2e-9)
    sim = Sim(mesh; name="relax", driver="SD")
    set_Ms(sim, 8e5)

    add_exch(sim, 1.3e-12)
    add_demag(sim)
    if I == 1
        rm_demag_charges(sim, 8e5, x=(-1, 1))
    elseif I == 2
        rm_demag_charges(sim, 8e5, y=(-1, 1))
    elseif I == 3
        rm_demag_charges(sim, 8e5, z=(-1, 1))
    end

    init_m0(sim, (1,1,1))
    relax(sim; max_steps=20000, stopping_dmdt=0.01)
    m = Array(sim.spin)
    m = reshape(m, 3, nx, ny, nz)

    mp1 = sqrt(1-m[I, 1, 1, 1]^2)
    mp2 = sqrt(1-m[I, 1, end, 1]^2)
    mp3 = sqrt(1-m[I, end, 1, 1]^2)
    mp4 = sqrt(1-m[I, end, end, 1]^2)

    println(mp1, " ", mp3, " ", mp3, " ", mp4) 

    @test max(mp1, mp2, mp3, mp4) < 5e-3
    #m = MicroMagnetic.average_m(sim)
    #save_vtk(sim, "3.vts", fields=["demagc"])
end

function test_rm_charges_xyz()
    relax_rm_charges(nx=30, ny=8, nz=3)
    relax_rm_charges(nx=9, ny=30, nz=3)
    relax_rm_charges(nx=3, ny=8, nz=30)
end


@using_gpu()
test_functions("DemagCharges", test_rm_charges_xyz)