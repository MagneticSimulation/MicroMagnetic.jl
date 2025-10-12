using Printf

function set_Ms(sim::MicroSimFE, Ms::NumberOrArray)
    sim.mu0_Ms .= mu_0 * Ms
    compute_L_Ms!(sim.L_mu, sim.mesh, sim.mu0_Ms)
    return true
end

function set_Ms(sim::MicroSimFE, Ms::Number; region_id=1)
    mesh = sim.mesh
    mu0_Ms = Array(sim.mu0_Ms)
    for i in 1:mesh.number_cells
        if mesh.region_ids[i] == region_id
            mu0_Ms[i] = mu_0 * Ms
        end
    end
    sim.mu0_Ms .= mu0_Ms
    compute_L_Ms!(sim.L_mu, sim.mesh, sim.mu0_Ms)
    return true
end

function init_m0(sim::MicroSimFE, m0::TupleOrArrayOrFunction; norm=true)
    init_vector!(sim.prespin, sim.mesh, m0)
    if norm
        normalise(sim.prespin, sim.n_total)
    end

    if any(isnan, sim.prespin)
        error("NaN is given by the input m0!")
    end

    sim.spin[:] .= sim.prespin[:]
    if isdefined(sim.driver, :integrator) && sim.driver.integrator isa AdaptiveRK
        set_initial_condition!(sim, sim.driver.integrator)
    end
end

function average_m(sim::MicroSimFE)
    b = sim.spin .* sim.L_mu

    b = reshape(b, 3, sim.n_total)

    mx = 3 * sum(b[1, :]) / sum(sim.L_mu)
    my = 3 * sum(b[2, :]) / sum(sim.L_mu)
    mz = 3 * sum(b[3, :]) / sum(sim.L_mu)

    return (mx, my, mz)
end

"""
    add_anis(sim::MicroSimFE, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::MicroSimFE, Ku::NumberOrArrayOrFunction; axis=(0, 0, 1),
                  name="anis")
    mesh = sim.mesh
    Kus = zeros(Float64, sim.n_cells)
    init_scalar!(Kus, sim.mesh, Ku)

    N = sim.n_total
    T = Float[]
    field = create_zeros(3*N)
    energy = create_zeros(N)

    if isa(axis, Tuple)
        @info "The uniform anisotropy axis is used in the simulation!"
        lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
        axes = zeros(Float64, 3 * sim.n_cells)
        for i in 1:(sim.n_cells)
            axes[3 * (i - 1) + 1] = axis[1] / lt
            axes[3 * (i - 1) + 2] = axis[2] / lt
            axes[3 * (i - 1) + 3] = axis[3] / lt
        end
    elseif isa(axis, Array)
        @info "The material-dependent anisotropy axes are used in the simulation!"
        axes = zeros(Float64, 3 * sim.n_cells)
        for i in 1:(sim.n_cells)
            id = mesh.material_ids[i]
            axes[3 * (i - 1) + 1] = axis[1, id]
            axes[3 * (i - 1) + 2] = axis[2, id]
            axes[3 * (i - 1) + 3] = axis[3, id]
        end
    end

    K_matrix = spzeros(3 * N, 3 * N)

    anis = AnisotropyFE(Kus, axes, field, energy, K_matrix, name)
    push!(sim.interactions, anis)

    assemble_anis_matirx(anis, sim)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(o.interactions[id].energy)))
    end
    return anis
end


"""
    add_exch(sim::MicroSimFEM, A::NumberOrArrayOrFunction; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::MicroSimFE, A::NumberOrArrayOrFunction; name="exch", method=:Direct)
    N = sim.n_total

    Spatial_A = zeros(Float64, sim.n_cells)

    #initialize a sparse matrix on CPU
    K_mat = spzeros(3 * N, 3 * N)

    T = Float[]
    field = create_zeros(3*N)
    energy = create_zeros(N)

    mass_matrix = assemble_mass_matirx(sim.mesh)

    exch = ExchangeFE(Spatial_A, field, energy, K_mat, mass_matrix, method, name)

    init_scalar!(Spatial_A, sim.mesh, A)

    # note that exch.K_matrix will be come to the GPU version if necessary
    assemble_exch_matirx(exch, sim)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(o.interactions[id].energy)))
    end

    return exch
end

"""
method could be "bem", "bem_hmatrix", "direct", "direct_hmatrix"
"""
function add_demag(sim::MicroSimFE; name="demag", method="bem", kwargs...)
    if method == "direct" || method == "direct_hmatrix"
        #demag = init_hmatix_demag(sim, method; kwargs...)
    else
        demag = init_demag(sim, method; kwargs...)
    end

    push!(sim.interactions, demag)
    demag.name = name

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(o.interactions[id].energy)))
    end

    return demag
end


function save_inp(sim::MicroSimFE, fname::String)
    mesh = sim.mesh

    #create the parent folder if necessary. TODO: copy this line to other places.
    mkpath(dirname(fname))

    f = open(fname, "w")

    data_number = 1 # FIX ME
    write(f, @sprintf("%d %d %d 0 0\n", mesh.number_nodes, mesh.number_cells, data_number))# data number
    for n in 1:(mesh.number_nodes)
        write(f,
              @sprintf("%d %0.12g %0.12g %0.12g\n", n, mesh.coordinates[1, n],
                       mesh.coordinates[2, n], mesh.coordinates[3, n]))
    end

    for c in 1:(mesh.number_cells)
        write(f,
              @sprintf("%d %d tet %d %d %d %d\n", c, 1, #FIX ME: material id
                       mesh.cell_verts[1, c], mesh.cell_verts[2, c], mesh.cell_verts[3, c],
                       mesh.cell_verts[4, c]))
    end

    write(f, @sprintf("%d", data_number))
    for i in 1:data_number
        write(f, " 3")
    end
    write(f, "\nm, none")

    m = reshape(sim.spin, 3, mesh.number_nodes)
    # FIX ME: to add more data
    for n in 1:(mesh.number_nodes)
        write(f, @sprintf("\n%d %0.12g %0.12g %0.12g", n, m[1, n], m[2, n], m[3, n]))
    end

    return close(f)
end

function save_inp_field(sim::MicroSimFE, fname::String, field::Array)
    mesh = sim.mesh

    #create the parent folder if necessary. TODO: copy this line to other places.
    mkpath(dirname(fname))

    f = open(fname, "w")

    data_number = 1 # FIX ME
    write(f, @sprintf("%d %d %d 0 0\n", mesh.number_nodes, mesh.number_cells, data_number))# data number
    for n in 1:(mesh.number_nodes)
        write(f,
              @sprintf("%d %0.12g %0.12g %0.12g\n", n, mesh.coordinates[1, n],
                       mesh.coordinates[2, n], mesh.coordinates[3, n]))
    end

    for c in 1:(mesh.number_cells)
        write(f,
              @sprintf("%d %d tet %d %d %d %d\n", c, 1, #FIX ME: material id
                       mesh.cell_verts[1, c], mesh.cell_verts[2, c], mesh.cell_verts[3, c],
                       mesh.cell_verts[4, c]))
    end

    write(f, @sprintf("%d", data_number))
    for i in 1:data_number
        write(f, " 3")
    end
    write(f, "\nm, none")

    m = reshape(field, 3, mesh.number_nodes)
    # FIX ME: to add more data
    for n in 1:(mesh.number_nodes)
        write(f, @sprintf("\n%d %0.12g %0.12g %0.12g", n, m[1, n], m[2, n], m[3, n]))
    end

    return close(f)
end