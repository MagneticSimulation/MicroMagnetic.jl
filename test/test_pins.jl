using MicroMagnetic
using Test

include("test_utils.jl")

function init_mm(i,j,k,dx,dy,dz)
    if i == 1
        return (0,0,1)
    elseif  i == 100
        return (0,1,0)
    else
        return (1,1,1)
    end
end

function pinning_boundary(i,j,k,dx,dy,dz)
    if i == 1 || i == 100
        return true
    end
    return false
end

function relax_system(;driver="SD")
    mesh =  FDMesh(nx=100, ny=1, nz=1, dx=2e-9, dy=2e-9, dz=1e-9)
	sim = Sim(mesh, name="test_pinning", driver=driver)
    if driver == "LLG"
        sim.driver.precession = false
        sim.driver.alpha = 0.5
    end
    set_Ms(sim, 8.6e5)
    set_pinning(sim, pinning_boundary)

    add_exch(sim, 1.3e-11)
    #add_anis(sim, 1e5, axis=(1,0,0))

    init_m0(sim, init_mm)
    relax(sim, maxsteps=10000, stopping_dmdt=0.1)
    save_vtk(sim, "pinning.vts")
    m = Array(sim.spin)

    @test abs(m[1])<1e-15
    @test abs(m[2])<1e-15
    @test abs(m[3]-1)<1e-15

end

@using_gpu()
test_functions("Pinning", relax_system)