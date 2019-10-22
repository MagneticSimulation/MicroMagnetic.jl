using Dates

mutable struct NEB_LLG_Driver <: NEBDriver
    nsteps::Int64
    ode::Integrator
end

function neb_llg_call_back(neb::NEB, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)
  effective_field_NEB(neb, spin)
  neb.field[:, 1] .= 0
  neb.field[:, neb.N] .= 0
  field = reshape(neb.field, 3*neb.nxyz)
  llg_rhs(dm_dt, spin, field, 1.0, 2.21e5, false, neb.nxyz)
  return nothing
end

function relax_NEB_LLG(neb::NEB, maxsteps::Int64, stopping_dmdt::Float64, save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String,save_ovf_every::Int64,ovf_folder::String)
    N = neb.N
    sim = neb.sim

    dmdt_factor = (2 * pi / 360) * 1e9

    driver = neb.driver

    if save_m_every>0
        compute_system_energy(neb)
        write_data(neb, neb.saver_energy)
        write_data(neb, neb.saver_distance)
    end
    if save_vtk_every > 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, 0)))
    end
    if save_ovf_every > 0
        save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, 0)))
    end

    rk_data = neb.driver.ode
    for i=1:maxsteps

        if !advance_step(neb, rk_data)
            @info("Advance step failed, end.")
            break
        end

        neb.driver.nsteps = rk_data.nsteps

        step_size = rk_data.step
        max_dmdt = compute_dmdt(neb.prespin, neb.spin, neb.nxyz, step_size)
        @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e  time=%s",
                     i, rk_data.step, rk_data.t, max_dmdt/dmdt_factor, Dates.now())

        if save_m_every>0 && i%save_m_every == 0
            compute_system_energy(neb)
            write_data(neb, neb.saver_energy)
            write_data(neb, neb.saver_distance)
        end
        if save_vtk_every > 0 && i%save_vtk_every == 0
            save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
        end
        if save_ovf_every > 0 && i%save_ovf_every == 0
            save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
        end

        if max_dmdt < stopping_dmdt*dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)

            if save_m_every>0
                compute_system_energy(neb)
                write_data(neb, neb.saver_energy)
                write_data(neb, neb.saver_distance)
            end
            if save_vtk_every > 0
                save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
            end
            if save_ovf_every > 0
                save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
            end
            break
        end
    end
    return nothing
end
