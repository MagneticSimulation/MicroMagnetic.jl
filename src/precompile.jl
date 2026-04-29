using PrecompileTools
import MicroMagnetic

function setup_workload()
    # atomistic
    @setup_workload begin
        mesh = MicroMagnetic.CubicMesh(nx=16, ny=16, nz=16, dx=1e-9, dy=1e-9, dz=1e-9, pbc="xyz")
        @compile_workload begin
            sim = MicroMagnetic.Sim(mesh; driver="SD", name="_precompile")
            MicroMagnetic.set_mu_s(sim, 1.0)
            MicroMagnetic.init_m0(sim, (1, 0, 0))
            MicroMagnetic.add_exch(sim, 1; name="exch", J2=-0.164, J3=0, J4=-0.082)
            MicroMagnetic.add_zeeman(sim, (1.0, 0, 0))
            MicroMagnetic.relax(sim; max_steps=5, stopping_dmdt=0.0)
        end
    end

    # finite difference
    @setup_workload begin
        mesh = MicroMagnetic.FDMesh(; dx=4e-9, dy=4e-9, dz=4e-9, nx=8, ny=8, nz=8)
        @compile_workload begin
            
            sim = MicroMagnetic.Sim(mesh; driver="SD", name="_precompile")
            MicroMagnetic.set_Ms(sim, 8e5)
            MicroMagnetic.init_m0(sim, (1, 0.25, 0.1))
            MicroMagnetic.add_exch(sim, 1.3e-11)
            MicroMagnetic.add_demag(sim)
            MicroMagnetic.add_anis(sim, 5e2, axis=(0, 1, 0))
            MicroMagnetic.add_dmi(sim, 3e-3; type="interfacial")
            MicroMagnetic.add_zeeman(sim, (0, 0, 1e4))
            MicroMagnetic.relax(sim; max_steps=3, stopping_dmdt=0.0)

            MicroMagnetic.set_driver(sim; driver="LLG", alpha=0.02)
            MicroMagnetic.add_sot(sim, 1e11, 0, (0, 1, 0))
            MicroMagnetic.add_stt(sim, model=:zhang_li, P=0.5, Ms=8e5, xi=0.05, J=(1e12,0,0))
            MicroMagnetic.add_thermal_noise(sim, 300.0)
            
            MicroMagnetic.run_sim(sim; steps=1, dt=1e-12, save_m_every=1)
        end
    end
end