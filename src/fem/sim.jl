using Printf
export interpolate_field

function set_Ms(sim::MicroSimFE, Ms::NumberOrArrayOrFunction; region_id=-1)
    Ms_array = zeros(Float64, sim.n_cells)
    init_scalar!(Ms_array, sim.mesh, Ms)

    mu0_Ms = Array(sim.mu0_Ms)
    if region_id >= 0
        for i in 1:mesh.number_cells
            if mesh.region_ids[i] == region_id
                mu0_Ms[i] = mu_0 * Ms_array[i]
            end
        end
    else
        mu0_Ms .= mu_0 .* Ms_array
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
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(anis.energy)))
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

    assemble_exch_matirx(exch, sim)

    push!(sim.interactions, exch)

    if sim.save_data
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(exch.energy)))
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
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::MicroSimFE -> sum(demag.energy)))
    end

    return demag
end

@doc raw"""
    add_exch_int(sim::MicroSimFE, J::Number; k1=1, k2=2, name="rkky")

Add an RKKY-type exchange interaction between layers. The energy of RKKY-type exchange is defined as

```math
E_\mathrm{rkky} =  - \int_\Gamma J_\mathrm{rkky} \mathbf{m}_{i} \cdot \mathbf{m}_{j} dA
```
where $\Gamma$ is the interface between magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$,
and $J_\mathrm{rkky}$ is the coupling constant related to the spacer layer thickness.

k1 and k2 are the surface IDs of the two layers, which can be specified in Netgen's .geo file, for example:
```
solid p1 = plane (0, 0, 1; 0, 0, -1) -bc=1;
solid p2 = plane (0, 0, 2; 0, 0, 1) -bc=5;
```

Suppose the RKKY coupling between p1 and p2 is J = -1e-4 J/m^2. We can use the following example to add the RKKY interaction:

**Examaple**
```julia
add_exch_int(sim, -1e-4; k1=1, k2=5)
```
"""
function add_exch_int(sim::MicroSimFE, J::Number; k1=1, k2=2, name="exch_int")
    n_total = sim.n_total
    field = create_zeros(3 * n_total)
    energy = create_zeros(n_total)

    N = sim.n_total
    #initialize a sparse matrix on CPU
    K_mat = spzeros(3 * N, 3 * N)

    rkky = InterlayerExchangeFE(J, Int32(k1), Int32(k2), field, energy, K_mat, name)

    assemble_rkky_matirx(rkky, sim)

    push!(sim.interactions, rkky)

    if sim.save_data
        push!(sim.saver.items, SaverItem(string("E_", name), "<J>", o::AbstractSim -> sum(rkky.energy)))
    end
    return rkky
end


"""
    interpolate_field(mesh::FEMesh, field::AbstractVector{T}, points::AbstractMatrix{T}) where T 

Interpolate a field (vector or scalar) at multiple points.
    
Parameters:
- mesh: The FEMesh
- field: The field data (vector: size 3*N, scalar: size N)
- points: 3×M matrix of point coordinates

**Examples:**
```julia
mesh = FEMesh("meshes/cylinder.mesh")
spin = Array(sim.spin)
points = hcat([0.5, 0.4, 1.1], [1.6, 0.5, 2.4])
f = interpolate_field(mesh, spin, points)
```
"""
function interpolate_field(mesh::FEMesh, field::AbstractVector{T}, 
                          points::AbstractMatrix{T}) where T

    # Precompute tetrahedron centers once (always needed)
    centers = compute_tetrahedron_centers(mesh)
    
    # Determine if field is vector or scalar
    is_vector_field = length(field) == 3 * mesh.number_nodes
    M = size(points, 2)
    
    # Initialize output
    if is_vector_field
        values = zeros(3, M)
    else
        values = zeros(M)
    end
    
    # Build KDTree once for all points
    kdtree = KDTree(centers)
    
    # Process each point
    for i in 1:M
        point = points[:, i]
        
        # Find containing tetrahedron using KDTree
        cell_id, bary = find_containing_tetrahedron(mesh, point, kdtree)
        
        if cell_id == -1
            cell_id, bary = find_containing_tetrahedron_all_cells(mesh, point)
            if cell_id == -1
                error("Point (x=", point[1], ", y=", point[2], ", z=", point[3], ") is outside the mesh")
            end
        end
        
        α, β, γ, δ = bary
        v_ids = mesh.cell_verts[:, cell_id]
        
        if is_vector_field
            # Vector field interpolation
            for comp in 1:3
                f1 = field[3*(v_ids[1]-1) + comp]
                f2 = field[3*(v_ids[2]-1) + comp]
                f3 = field[3*(v_ids[3]-1) + comp]
                f4 = field[3*(v_ids[4]-1) + comp]
                values[comp, i] = α*f1 + β*f2 + γ*f3 + δ*f4
            end
        else
            # Scalar field interpolation
            f1 = field[v_ids[1]]
            f2 = field[v_ids[2]]
            f3 = field[v_ids[3]]
            f4 = field[v_ids[4]]
            values[i] = α*f1 + β*f2 + γ*f3 + δ*f4
        end
    end
    
    return values
end

function MicroMagnetic.interpolate_field(mesh::FEMesh, field::AbstractVector{T}, points::Vector{Vector{T}}) where T
    M = length(points)
    point_matrix = Matrix{T}(undef, 3, M)
    for i in 1:M
        point_matrix[:, i] = points[i]
    end
    return interpolate_field(mesh, field, point_matrix)
end

function interpolate_field(mesh::FEMesh, field::AbstractVector{T}, point::AbstractVector{T}) where T
    points = reshape(point, 3, 1)
    result = interpolate_field(mesh, field, points)
    return result[:, 1]
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