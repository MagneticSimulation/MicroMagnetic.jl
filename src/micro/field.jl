# Unpack a Zeeman time modulation into three scalar factors (a scalar return
# applies to all components) and record them for the saver display.
function _zeeman_time_factors(zee::Zeeman, t::Float64)
    f = zee.ft(t)
    if f isa Tuple
        zee.time_fx = f[1]
        zee.time_fy = f[2]
        zee.time_fz = f[3]
    else
        zee.time_fx = f
        zee.time_fy = f
        zee.time_fz = f
    end
    return zee.time_fx, zee.time_fy, zee.time_fz
end

function effective_field(zee::Zeeman, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    factor = sim.mesh.volume
    tx, ty, tz = _zeeman_time_factors(zee, t)

    back = default_backend[]
    if zee.ft === _static_time
        # the interaction field was materialised at add/update time and never
        # changes, so only the energy needs refreshing
        zeeman_energy_kernel!(back, groupsize[])(spin, zee.field, zee.energy,
                                                 sim.mu0_Ms, T(factor); ndrange=N)
    else
        zeeman_field_kernel!(back, groupsize[])(spin, zee.field, zee.energy,
                                                sim.mu0_Ms, zee.Hx, zee.Hy, zee.Hz,
                                                T(factor), T(tx), T(ty), T(tz);
                                                ndrange=N)
    end
    return nothing
end

function effective_field(anis::Anisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    back = default_backend[]
    anisotropy_kernel!(back, groupsize[])(spin, anis.field, anis.energy, anis.Ku,
                                          anis.axis_x, anis.axis_y, anis.axis_z,
                                          sim.mu0_Ms, volume; ndrange=N)

    return nothing
end

function effective_field(anis::CubicAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    heff = output === nothing ? anis.field : output
    cubic_anisotropy_kernel!(default_backend[])(spin, heff, anis.energy, anis.Kc,
                                                anis.axis1x, anis.axis1y, anis.axis1z,
                                                anis.axis2x, anis.axis2y, anis.axis2z,
                                                anis.axis3x, anis.axis3y, anis.axis3z,
                                                sim.mu0_Ms, volume;
                                                ndrange=N)

    return nothing
end

function effective_field(anis::HexagonalAnisotropy, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    hexagonal_anisotropy_kernel!(default_backend[])(spin, anis.field, anis.energy, anis.K1,
                                                    anis.K2, anis.K3, sim.mu0_Ms, volume;
                                                    ndrange=N)

    return nothing
end

function effective_field(anis::TwinMonoclinicAnisotropy, sim::MicroSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)

    twin_monoclinic_anisotropy_kernel!(default_backend[], groupsize[])(
        spin, anis.field, anis.energy, anis.axis_a, anis.axis_b, anis.axis_u111,
        sim.mu0_Ms, T(anis.Ka), T(anis.Kb), T(anis.Kaa), T(anis.Kbb),
        T(anis.Kab), T(anis.Ku), volume; ndrange=N)

    return nothing
end

function effective_field(exch::Exchange, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    n_total = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    back = default_backend[]
    exchange_kernel!(back, groupsize[])(spin, exch.field, exch.energy, sim.mu0_Ms,
                                        exch.Ax, exch.Ay, exch.Az,
                                        dx, dy, dz, mesh.ngbs, volume; ndrange=n_total)

    return nothing
end

function effective_field(dmi::BulkDMI, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    volume = T(mesh.volume)

    dx, dy, dz = T(mesh.dx), T(mesh.dy), T(mesh.dz)
    tfac = T(dmi.ft(t))
    back = default_backend[]
    if dmi.type === :bulk
        bulkdmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy, sim.mu0_Ms,
                                           dmi.Dx, dmi.Dy, dmi.Dz, dx, dy, dz,
                                           mesh.ngbs, volume, tfac; ndrange=N)
    else
        interfacial_dmi_kernel!(back, groupsize[])(spin, dmi.field, dmi.energy,
                                                   sim.mu0_Ms, dmi.Dx, dx, dy, dz,
                                                   mesh.ngbs, volume, tfac; ndrange=N)
    end

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

function effective_field(st::StochasticField, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = sim.mesh.volume
    integrator = sim.driver.integrator

    if integrator.nsteps > st.nsteps
        randn!(st.eta)
        st.nsteps = integrator.nsteps
    end

    dt = integrator.step
    gamma = sim.driver.gamma
    k_B = st.k_B

    back = default_backend[]
    # alpha enters per spin inside the kernel (uniform damping arrives as a Fill)
    factor = 2 * k_B / (volume * gamma * dt)
    if st.spatiotemporal_mode
        mesh = sim.mesh
        Nx, Ny, Nz = mesh.nx, mesh.ny, mesh.nz
        temp = reshape(st.temperature, Nx, Ny, Nz)
        spatiotemporal_kernel!(back, groupsize[])(temp, mesh.dx, mesh.dy, mesh.dz, mesh.x0,
                                                  mesh.y0, mesh.z0, T(t), st.scaling_fun;
                                                  ndrange=(Nx, Ny, Nz))
    else
        scaling_factor = st.scaling_fun(t)
        factor = factor * scaling_factor
        st.scaling_factor = scaling_factor
    end

    stochastic_field_kernel!(back, groupsize[])(spin, st.field, st.energy, sim.mu0_Ms,
                                                st.eta, st.temperature, st.offset_temp,
                                                as_param_array(sim.driver.alpha, N),
                                                T(factor), T(volume); ndrange=N)

    return nothing
end

function effective_field(torque::TorqueField, sim::AbstractSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma

    torque.torque_fun(torque.field, spin, t)

    back = default_backend[]
    torque_kernel!(back, groupsize[])(spin, torque.field, gamma; ndrange=N)

    return nothing
end

function effective_field(me::Magnetoelastic, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    volume = T(sim.mesh.volume)
    
    if me.model == :tensor
        back = default_backend[]
        magnetoelastic_tensor_kernel!(back, groupsize[])(spin, me.field, me.energy, me.field_data,
                                                       T(me.lambda_s),
                                                       sim.mu0_Ms, volume; ndrange=N)
    elseif me.model == :cubic
        back = default_backend[]
        magnetoelastic_cubic_kernel!(back, groupsize[])(spin, me.field, me.energy, me.field_data,
                                                  T(me.B1), T(me.B2),
                                                  sim.mu0_Ms, volume; ndrange=N)
    else
        throw(ArgumentError("Unknown magnetoelastic model: $(me.model)"))
    end

    return nothing
end

function effective_field(torque::DFTorqueField, sim::AbstractSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma

    back = default_backend[]
    df_torque_kernel!(back, groupsize[])(spin, torque.field, gamma, torque.aj, torque.bj,
                                         torque.px, torque.py, torque.pz; ndrange=N)

    return nothing
end

function effective_field(torque::ZhangLiTorque, sim::AbstractSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma
    mesh = sim.mesh

    ut = T(torque.ufun(t)/gamma)

    back = default_backend[]
    zhangli_torque_kernel!(back, groupsize[])(spin, torque.field, torque.bJ, mesh.ngbs,
                                              torque.xi, ut, T(mesh.dx), T(mesh.dy),
                                              T(mesh.dz); ndrange=N)

    return nothing
end

function effective_field(torque::SlonczewskiTorque, sim::AbstractSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total

    lambda_sq = T(torque.Lambda^2)
    ft = T(torque.ufun(t)*torque.beta)

    back = default_backend[]
    slonczewski_torque_kernel!(back, groupsize[])(spin, torque.field, torque.J, lambda_sq,
                                                  torque.P, torque.xi, ft, torque.px,
                                                  torque.py, torque.pz; ndrange=N)

    return nothing
end

function effective_field(torque::SAHETorqueField, sim::AbstractSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    gamma = sim.driver.gamma

    back = default_backend[]
    ms = isa(sim, MicroSim) ? sim.mu0_Ms : sim.mu_s
    sahe_torque_kernel!(back, groupsize[])(spin, torque.field, ms, torque.sigma_s_a1, torque.sigma_sa1_a1,
                                           torque.sigma_sa2, torque.a2, torque.a3,
                                           gamma, torque.beta;
                                           ndrange=N)

    return nothing
end

#we keep this function for debug and testing purpose, only works on CPU.
#It mirrors exchange_kernel! line by line so the two can be compared tightly.
function effective_field_debug(exch::Exchange, sim::MicroSim, spin::Array{Float64,1},
                               t::Float64)
    mesh = sim.mesh
    dx = mesh.dx
    dy = mesh.dy
    dz = mesh.dz
    ngbs = mesh.ngbs
    n_total = sim.n_total
    field = exch.field
    energy = exch.energy
    Axs = Array(exch.Ax)
    Ays = Array(exch.Ay)
    Azs = Array(exch.Az)
    Ms = sim.mu0_Ms
    nx = 2.0 / (dx * dx)
    ny = 2.0 / (dy * dy)
    nz = 2.0 / (dz * dz)

    for index in 1:n_total
        i = 3 * index - 2
        if Ms[index] == 0.0
            energy[index] = 0.0
            field[i] = 0.0
            field[i + 1] = 0.0
            field[i + 2] = 0.0
            continue
        end
        Ax_I, Ay_I, Az_I = Axs[index], Ays[index], Azs[index]
        fx, fy, fz = 0.0, 0.0, 0.0

        # ---- (±x) ----
        for j in 1:2
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Ax_n = Axs[id]
                if Ax_I != 0 && Ax_n != 0
                    k = 3 * id - 2
                    A_eff = 2 * Ax_I * Ax_n / (Ax_I + Ax_n)
                    fx += A_eff * nx * (spin[k] - spin[i])
                    fy += A_eff * nx * (spin[k + 1] - spin[i + 1])
                    fz += A_eff * nx * (spin[k + 2] - spin[i + 2])
                end
            end
        end

        # ---- (±y) ----
        for j in 3:4
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Ay_n = Ays[id]
                if Ay_I != 0 && Ay_n != 0
                    k = 3 * id - 2
                    A_eff = 2 * Ay_I * Ay_n / (Ay_I + Ay_n)
                    fx += A_eff * ny * (spin[k] - spin[i])
                    fy += A_eff * ny * (spin[k + 1] - spin[i + 1])
                    fz += A_eff * ny * (spin[k + 2] - spin[i + 2])
                end
            end
        end

        # ---- (±z) ----
        for j in 5:6
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Az_n = Azs[id]
                if Az_I != 0 && Az_n != 0
                    k = 3 * id - 2
                    A_eff = 2 * Az_I * Az_n / (Az_I + Az_n)
                    fx += A_eff * nz * (spin[k] - spin[i])
                    fy += A_eff * nz * (spin[k + 1] - spin[i + 1])
                    fz += A_eff * nz * (spin[k + 2] - spin[i + 2])
                end
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

#we keep this function for debug and testing purpose, only works on CPU.
#It mirrors bulkdmi_kernel! line by line so the two can be compared tightly.
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
    Dxs = Array(dmi.Dx)
    Dys = Array(dmi.Dy)
    Dzs = Array(dmi.Dz)
    axes6 = (1.0 / dx, -1.0 / dx, 1.0 / dy, -1.0 / dy, 1.0 / dz, -1.0 / dz)

    for index in 1:n_total
        i = 3 * index - 2
        if Ms[index] == 0.0
            energy[index] = 0.0
            field[i] = 0.0
            field[i + 1] = 0.0
            field[i + 2] = 0.0
            continue
        end
        Dx_I = Dxs[index]
        Dy_I = Dys[index]
        Dz_I = Dzs[index]
        fx, fy, fz = 0.0, 0.0, 0.0

        # ---- (±x) ----
        for j in 1:2
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Dx_n = Dxs[id]
                if Dx_I != 0 && Dx_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dx_I * Dx_n / (Dx_I + Dx_n)
                    fy += D_eff * (-axes6[j] * spin[k + 2])
                    fz += D_eff * ( axes6[j] * spin[k + 1])
                end
            end
        end

        # ---- (±y) ----
        for j in 3:4
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Dy_n = Dys[id]
                if Dy_I != 0 && Dy_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dy_I * Dy_n / (Dy_I + Dy_n)
                    fx += D_eff * ( axes6[j] * spin[k + 2])
                    fz += D_eff * (-axes6[j] * spin[k])
                end
            end
        end

        # ---- (±z) ----
        for j in 5:6
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                Dz_n = Dzs[id]
                if Dz_I != 0 && Dz_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dz_I * Dz_n / (Dz_I + Dz_n)
                    fx += D_eff * (-axes6[j] * spin[k + 1])
                    fy += D_eff * ( axes6[j] * spin[k])
                end
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

  t : Time in second. This term is used when a time-modulated Zeeman is added to sim.

For example:

  ```julia
    #sim = Sim(CPUmesh)
    effective_field(sim, sim.spin, 0.0)
  ```

After running this function, the effective field is calculated and saved in sim.field.
"""
function effective_field(sim::AbstractSim, spin, t::Float64=0.0;)
    fill!(sim.field, 0.0)
    for interaction in sim.interactions
        # a static Zeeman contributes its precomputed field only; a time-modulated
        # one must refresh interaction.field every step
        if !(interaction isa Zeeman && interaction.ft === _static_time)
            @timeit timer interaction.name effective_field(interaction, sim, spin, t)
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
        if isa(sim.driver, LLG) && isa(sim.driver.integrator, GPSM)
            exch = sim.driver.integrator.exch
            effective_field(exch, sim, spin, t)
            sim.energy .+= exch.energy
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

function build_exch_matrix(exch::Exchange, sim::MicroSim)
    mesh = sim.mesh
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    ngbs = Array(mesh.ngbs)
    Ms = Array(sim.mu0_Ms)
    Axs = Array(exch.Ax)
    Ays = Array(exch.Ay)
    Azs = Array(exch.Az)

    nx = 2 / (dx * dx)
    ny = 2 / (dy * dy)
    nz = 2 / (dz * dz)

    n_total = sim.n_total

    T = Float[]
    I = Int64[]
    J = Int64[]
    V = T[]

    safe_div(a, b) = b == 0.0 ? 0.0 : a / b

    for index in 1:n_total
        if Ms[index] == 0.0
            continue
        end

        Ms_inv = 1.0 / Ms[index]
        diag_sum = 0.0
        for j in 1:6
            id = ngbs[j, index]
            if id > 0 && Ms[id] > 0
                # Effective exchange stiffness for the pair (index, id) along the
                # bond direction (1:2 → x, 3:4 → y, 5:6 → z)
                A_I, A_nb, nabla = j <= 2 ? (Axs[index], Axs[id], nx) :
                                   j <= 4 ? (Ays[index], Ays[id], ny) :
                                   (Azs[index], Azs[id], nz)
                A_eff = safe_div(2 * A_I * A_nb, A_I + A_nb)
                coeff = A_eff * nabla * Ms_inv

                diag_sum -= coeff
                push!(I, index)
                push!(J, id)
                push!(V, coeff)
            end
        end

        # Diagonal entry
        push!(I, index)
        push!(J, index)
        push!(V, diag_sum)
    end

    Laplaian = sparse(I, J, V, n_total, n_total)

    if default_backend[] != CPU()
        Laplaian = GPUSparseMatrixCSR[](Laplaian)
    end

    return Laplaian
end
