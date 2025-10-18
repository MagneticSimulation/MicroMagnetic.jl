function effective_field(zee::Zeeman, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    factor = sim.mesh.volume

    back = default_backend[]
    zeeman_kernel!(back, groupsize[])(spin, zee.field, zee.energy, sim.mu0_Ms, T(factor);
                                      ndrange=N)
    return nothing
end

function effective_field(zee::TimeZeeman, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    factor = sim.mesh.volume
    if zee.is_scalar
       tx = zee.time_fun(t)
       zee.time_fx = tx
       zee.time_fy = tx
       zee.time_fz = tx
    else
       tx, ty, tz = zee.time_fun(t)
       zee.time_fx = tx
       zee.time_fy = ty
       zee.time_fz = tz
    end

    kernel = time_zeeman_kernel!(default_backend[], groupsize[])
    kernel(spin, zee.field, zee.init_field, zee.energy, sim.mu0_Ms, T(factor), T(zee.time_fx), T(zee.time_fy), T(zee.time_fz);
            ndrange=N)
    return nothing
end

function effective_field(anis::Anisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    axis = anis.axis
    volume = T(sim.mesh.volume)

    back = default_backend[]
    anisotropy_kernel!(back, groupsize[])(spin, anis.field, anis.energy, anis.Ku, axis[1],
                                          axis[2], axis[3], sim.mu0_Ms, volume; ndrange=N)

    return nothing
end

function effective_field(anis::SpatialAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    back = default_backend[]
    spatial_anisotropy_kernel!(back, groupsize[])(spin, anis.field, anis.energy, anis.Ku, anis.axes, sim.mu0_Ms, volume; ndrange=N)

    return nothing
end

function effective_field(anis::CubicAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    a1x, a1y, a1z = anis.axis1
    a2x, a2y, a2z = anis.axis2
    a3x, a3y, a3z = anis.axis3
    heff = output === nothing ? anis.field : output
    cubic_anisotropy_kernel!(default_backend[])(spin, heff, anis.energy, anis.Kc,
                                                T(a1x), T(a1y), T(a1z), T(a2x), T(a2y),
                                                T(a2z), T(a3x), T(a3y), T(a3z), sim.mu0_Ms,
                                                volume; ndrange=N)

    return nothing
end

function effective_field(anis::HexagonalAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
            t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    hexagonal_anisotropy_kernel!(default_backend[])(spin, anis.field, anis.energy, anis.K1, 
            anis.K2, anis.K3, sim.mu0_Ms, volume; ndrange=N)

    return nothing
end

function effective_field(exch::SpatialExchange, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    n_total = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    exchange_kernel!(back, groupsize[])(spin, exch.field, exch.energy, sim.mu0_Ms, exch.A,
                                        dx, dy, dz, mesh.ngbs, volume; ndrange=n_total)
    
    return nothing
end

function effective_field(exch::UniformExchange, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    uniform_exchange_kernel!(back, groupsize[])(spin, exch.field, exch.energy, sim.mu0_Ms,
                                                exch.Ax, exch.Ay, exch.Az, dx, dy, dz,
                                                mesh.ngbs, volume; ndrange=N)
    
    return nothing
end

function effective_field(dmi::BulkDMI, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    bulkdmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy, sim.mu0_Ms, dmi.Dx,
                                       dmi.Dy, dmi.Dz, dx, dy, dz, mesh.ngbs, volume;
                                       ndrange=N)

    return nothing
end

function effective_field(dmi::SpatialBulkDMI, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    spatial_bulkdmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy, sim.mu0_Ms,
                                               dmi.D, dx, dy, dz, mesh.ngbs, volume;
                                               ndrange=N)

    
    return nothing
end

function effective_field(dmi::InterfacialDMI, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    interfacial_dmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy, sim.mu0_Ms,
                                               dmi.D, dx, dy, dz, mesh.ngbs, volume;
                                               ndrange=N)

    return nothing
end

function effective_field(dmi::InterlayerDMI, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)
    nx, ny = mesh.nx, mesh.ny
    dz = T(mesh.dz)

    back = default_backend[]
    interlayer_dmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy, sim.mu0_Ms,
                                              dmi.Dx, dmi.Dy, dmi.Dz, dmi.k1, dmi.k2,
                                              Int32(nx), Int32(ny), dz, volume;
                                              ndrange=(nx, ny))

    return nothing
end

function effective_field(exch::InterlayerExchange, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)
    nx, ny = mesh.nx, mesh.ny
    dz = T(mesh.dz)

    back = default_backend[]
    interlayer_exch_kernel!(back, groupsize[])(spin, exch.field, exch.energy, sim.mu0_Ms,
                                               exch.Js, exch.k1, exch.k2, Int32(nx),
                                               Int32(ny), dz, volume; ndrange=(nx, ny))

    return nothing
end

function effective_field(stochastic::StochasticField, sim::MicroSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = sim.mesh.volume
    integrator = sim.driver.integrator

    if integrator.nsteps > stochastic.nsteps
        randn!(stochastic.eta)
        stochastic.nsteps = integrator.nsteps
    end

    dt = integrator.step
    gamma = sim.driver.gamma
    alpha = sim.driver.alpha
    k_B = stochastic.k_B
    scaling_factor = stochastic.scaling_fun(t)
    factor = 2 * alpha * k_B / (volume * gamma * dt) * scaling_factor

    stochastic.scaling_factor = scaling_factor

    back = default_backend[]
    stochastic_field_kernel!(back, groupsize[])(spin, stochastic.field, stochastic.energy,
                                                sim.mu0_Ms, stochastic.eta,
                                                stochastic.temperature, factor, T(volume);
                                                ndrange=N)


    return nothing
end

function effective_field(torque::DFTorqueField, sim::AbstractSim, spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma
    
    back = default_backend[]
    ms = isa(sim, MicroSim) ? sim.mu0_Ms : sim.mu_s
    df_torque_kernel!(back, groupsize[])(spin, torque.field, ms, gamma, torque.aj, 
                      torque.bj, torque.px, torque.py, torque.pz; ndrange=N)

    return nothing
end


function effective_field(torque::SAHETorqueField, sim::AbstractSim, spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma
    
    back = default_backend[]
    c2x, c2y, c2z = torque.c2[1], torque.c2[2], torque.c2[3]
    ms = isa(sim, MicroSim) ? sim.mu0_Ms : sim.mu_s
    sahe_torque_kernel!(back, groupsize[])(spin, torque.field, ms, gamma, torque.beta, 
                      torque.c1, T(c2x), T(c2y), T(c2z), torque.c3; ndrange=N)

    return nothing
end

#we keep this function for debug and testing purpose, only works on CPU
function effective_field_debug(exch::SpatialExchange, sim::MicroSim, spin::Array{Float64,1},
                               t::Float64)
    mesh = sim.mesh
    dx = mesh.dx
    dy = mesh.dy
    dz = mesh.dz
    ngbs = mesh.ngbs
    n_total = sim.n_total
    field = exch.field
    energy = exch.energy
    A = exch.A
    Ms = sim.mu0_Ms
    ax = 2.0 / (dx * dx)
    ay = 2.0 / (dy * dy)
    az = 2.0 / (dz * dz)
    nabla = (ax, ax, ay, ay, az, az)

    for index in 1:n_total
        i = 3 * index - 2
        if Ms[index] == 0.0
            energy[index] = 0.0
            field[i] = 0.0
            field[i + 1] = 0.0
            field[i + 2] = 0.0
            continue
        end
        fx, fy, fz = 0.0, 0.0, 0.0
        for j in 1:6
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                k = 3 * id - 2
                fx += A[index] * nabla[j] * (spin[k] - spin[i])
                fy += A[index] * nabla[j] * (spin[k + 1] - spin[i + 1])
                fz += A[index] * nabla[j] * (spin[k + 2] - spin[i + 2])
            end
        end
        Ms_inv = 1.0 / Ms[index]
        energy[index] = -0.5 *
                        (fx * spin[i] + fy * spin[i + 1] + fz * spin[i + 2]) *
                        mesh.volume
        field[i] = fx * Ms_inv
        field[i + 1] = fy * Ms_inv
        field[i + 2] = fz * Ms_inv
    end
end

#we keep this function for debug and testing purpose, only works on CPU
function effective_field_debug(dmi::BulkDMI, sim::MicroSim, spin::Array{Float64,1},
                               t::Float64)
    mesh = sim.mesh
    dx = mesh.dx
    dy = mesh.dy
    dz = mesh.dz
    ngbs = mesh.ngbs
    n_total = sim.n_total
    field = dmi.field
    energy = dmi.energy
    Ms = sim.mu0_Ms
    Dx, Dy, Dz = dmi.Dx, dmi.Dy, dmi.Dz
    Ds = (Dx / dx, Dx / dx, Dy / dy, Dy / dy, Dz / dz, Dz / dz)
    ax = (1.0, -1.0, 0.0, 0.0, 0.0, 0.0)
    ay = (0.0, 0.0, 1.0, -1.0, 0.0, 0.0)
    az = (0.0, 0.0, 0.0, 0.0, 1.0, -1.0)

    for index in 1:n_total
        i = 3 * index - 2
        if Ms[index] == 0.0
            energy[index] = 0.0
            field[i] = 0.0
            field[i + 1] = 0.0
            field[i + 2] = 0.0
            continue
        end
        fx, fy, fz = 0.0, 0.0, 0.0

        for j in 1:6
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                k = 3 * (id - 1) + 1
                fx += Ds[j] *
                      cross_x(ax[j], ay[j], az[j], spin[k], spin[k + 1], spin[k + 2])
                fy += Ds[j] *
                      cross_y(ax[j], ay[j], az[j], spin[k], spin[k + 1], spin[k + 2])
                fz += Ds[j] *
                      cross_z(ax[j], ay[j], az[j], spin[k], spin[k + 1], spin[k + 2])
            end
        end

        Ms_inv = 1.0 / Ms[index]
        energy[index] = -0.5 *
                        (fx * spin[i] + fy * spin[i + 1] + fz * spin[i + 2]) *
                        mesh.volume
        field[i] = fx * Ms_inv
        field[i + 1] = fy * Ms_inv
        field[i + 2] = fz * Ms_inv
    end
end

"""
    effective_field(sim::AbstractSim, spin::Array{Float64, 1}, t::Float64)

Calculate the effective field of a CPU Sim.

Parameters:

  sim : AbstractSim struct whose field is to be calculated.

  spin : 1-d array that matches sim.mesh. 

  t : Time in second. This term is used when TimeZeeman is added into sim.

For example:

  ```julia
    #sim = Sim(CPUmesh)
    effective_field(sim, sim.spin, 0.0)
  ```

After running this function, the effective field is calculated and saved in sim.field.
"""
function effective_field(sim::AbstractSim, spin, t::Float64=0.0)
    fill!(sim.field, 0.0)
    for interaction in sim.interactions
        if !isa(interaction, Zeeman)
            effective_field(interaction, sim, spin, t)
        end
        vector_add(sim.field, interaction.field)
    end
    return nothing
end

"""
function effective_field(sim::AbstractSim, spin, output)
    fill!(output, 0.0)
    for interaction in sim.interactions
        if !isa(interaction, Zeeman)
            effective_field(interaction, sim, spin, 0.0)
        end
        vector_add(output, interaction.field)
    end
    return nothing
end
"""

function compute_system_energy(sim::AbstractSim, spin::AbstractArray, t::Float64)
    @timeit timer "compute_system_energy" begin
        fill!(sim.energy, 0.0)
        for interaction in sim.interactions
            if hasproperty(interaction, :energy)
                effective_field(interaction, sim, spin, t)
                sim.energy .+= interaction.energy
            end
        end
    end
    return 0
end

function effective_field_energy(sim::AbstractSim, spin, t::Float64=0.0)
    fill!(sim.field, 0.0)
    fill!(sim.energy, 0.0)
    for interaction in sim.interactions
        effective_field(interaction, sim, spin, t)
        sim.field .+= interaction.field
        sim.energy .+= interaction.energy
    end
    return nothing
end
