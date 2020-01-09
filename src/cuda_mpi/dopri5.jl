#implement the  Runge-Kutta method with adaptive stepsize using Dormand–Prince pair
#https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Dormand%E2%80%93Prince

#To use this integrator, one needs to provide:
#   - A simulation instance with field "spin" and "prespin"
#   - A call back function
using MPI
using CuArrays

function dopri5_step_inner_GPU_MPI(neb::NEB_GPU_MPI, step::Float64, t::Float64)
    max_error = dopri5_step_inner_GPU(neb, step, t)
    all_max_error = [0.0]
    #println("rank = $(neb.comm_rank), max_error= $max_error")
    if neb.comm_rank == 0
        all_max_error[1] = MPI.Reduce(max_error, MPI.MAX, 0, MPI.COMM_WORLD)
    else
        MPI.Reduce(max_error, MPI.MAX, 0, MPI.COMM_WORLD)
    end
    MPI.Bcast!(all_max_error, 0, MPI.COMM_WORLD)
    #println("rank = $(neb.comm_rank), all_max_error= $all_max_error")
    return all_max_error[1]
end

function compute_init_step_GPU_MPI(neb::NEB_GPU_MPI, dt::Float64)
    step_next = compute_init_step_DP(neb, dt)
    all_step_next = [0.0]

    if neb.comm_rank == 0
        all_step_next[1] = MPI.Reduce(step_next, MPI.MIN, 0, MPI.COMM_WORLD)
    else
        MPI.Reduce(step_next, MPI.MIN, 0, MPI.COMM_WORLD)
    end
    MPI.Bcast!(all_step_next, 0, MPI.COMM_WORLD)
    return all_step_next[1]
end


function advance_step(sim::NEB_GPU_MPI, integrator::DormandPrinceGPU)

    max_nsteps = 100

    t = integrator.t

    sim.prespin .= sim.spin

    if integrator.step_next <= 0
        integrator.step_next = compute_init_step_GPU_MPI(sim, 1e-12)
        if integrator.step_next<1e-15
            #integrator.step_next = 1e-15
        end
    end

    step_next = integrator.step_next

    nstep = 1
    while true
        max_error = dopri5_step_inner_GPU_MPI(sim, step_next, t)/integrator.tol
        if isnan(max_error)
            step_next = 1e-14
        end
        nstep += 1
        if nstep > max_nsteps
            @info("Too many inner steps, the system may not converge!", step_next, max_error)
            return false
        end
        integrator.succeed = (max_error <= 1)

        if integrator.succeed
            integrator.nsteps += 1
            integrator.step = step_next
            integrator.t += integrator.step
            factor =  integrator.safety*(1.0/max_error)^0.2
            integrator.step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
            break
        else
            factor =  integrator.safety*(1.0/max_error)^0.25
            step_next = step_next*min(integrator.facmax, max(integrator.facmin, factor))
        end
    end
    return true
end