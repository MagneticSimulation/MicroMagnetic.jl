using JuMag
using Printf
using Test

function init_dw(i, j, k, dx, dy, dz)
  if k < 150
    return (0,0.1,1)
elseif k<160
   return (0,1,0.1)
  else
   return (0,0.1,-1)
  end
end

function relax_system()
    mesh =  TriangularMesh(nx=1, ny=1, nz=300, dz=1e-9)
    sim = Sim(mesh, name="dw_atomic")
    set_mu_s(sim, 2.0*mu_B)
    sim.driver.alpha = 0.5
    sim.driver.gamma = 1.76e11
    sim.driver.precession = false

    add_exch(sim, 40*meV)
    add_anis(sim, 0.005*40*meV, axis=(0,0,1))

    init_m0(sim, init_dw)

    dmdt_factor = (2 * pi / 360) * 1e9
    println(dmdt_factor)
    relax(sim, stopping_dmdt=0.1*dmdt_factor, maxsteps=2000, save_vtk_every=100, save_ovf_every=10000)

end


function run_dynamics(;alpha=0.1, beta=0.0, u=50)
    mesh =  TriangularMesh(nx=1, ny=1, nz=300, dz=1e-9)
    sim = Sim(mesh, name="dw_dyn", driver="LLG_STT")
    set_mu_s(sim, 2.0*mu_B)
    sim.driver.alpha = alpha
    sim.driver.beta = beta
    sim.driver.gamma = 1.76e11

    set_uz(sim, u)

    add_exch(sim, 40*meV)
    add_anis(sim, 0.005*40*meV, axis=(0,0,1))

    ovf = read_ovf("ovfs/dw_atomic_1199.ovf")
    init_m0(sim, ovf.data)

    for i=1:20
        run_until(sim, 1e-11*i)
        println("t = ", 1e-11*i)
        save_vtk(sim, "stt_dyn_"*string(i))
    end
end

relax_system()
run_dynamics()
