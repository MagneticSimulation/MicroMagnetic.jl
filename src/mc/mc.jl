using Random
using WriteVTK

export MonteCarlo

function MonteCarlo(mesh::Mesh; name="mc", mc_2d=false)
    if !isa(mesh, CubicMesh) && (!isa(mesh, TriangularMesh))
        error("Only support CubicMesh and TriangularMesh.")
    end

    sim = MonteCarlo{Float[]}()
    sim.mc_2d = mc_2d
    sim.mesh = mesh
    n_total = mesh.n_total
    sim.n_total = n_total
    sim.shape = create_ones(Bool, n_total)
    sim.spin = create_zeros(3 * n_total)
    sim.nextspin = create_zeros(3 * n_total)
    sim.rnd = create_zeros(3 * n_total)
    sim.energy = create_zeros(n_total)
    sim.delta_E = create_zeros(n_total)
    sim.steps = 0
    sim.name = name
    sim.T = 300

    F = Float[]
    D0 = create_zeros(3, mesh.n_ngbs)
    sim.exch = ExchangeDMI(F(0), F(0), F(0), D0)
    sim.zeeman = ZeemanMC(F(0), F(0), F(0))
    axis = (F(1), F(0), F(0))
    sim.anis = AnisotropyMC(F(0), axis, F(0))

    saver = init_saver(string(name, ".txt"), "MC")
    push!(saver.items, SaverItem("E_total", "<J>", o::AbstractSim -> sum(o.energy)))
    push!(saver.items,
          SaverItem(("m_x", "m_y", "m_z"), ("<>", "<>", "<>"), average_m))

    sim.saver = saver

    return sim
end

function add_exch(sim::MonteCarlo; J=0) 
    if isa(J, Tuple)
        sim.exch.Jx = J[1] / k_B
        sim.exch.Jy = J[2] / k_B
        sim.exch.Jz = J[3] / k_B
    else
        sim.exch.Jx = J / k_B
        sim.exch.Jy = J / k_B
        sim.exch.Jz = J / k_B
    end
    return nothing
end

function add_dmi(sim::MonteCarlo; D=0, Dz=0, type="bulk")
    mesh = sim.mesh
    D = D / k_B
    Dz = Dz / k_B

    Dij = zeros(Float[], (3, mesh.n_ngbs))
    if isa(mesh, TriangularMesh)
        for i in 1:6
            theta = (i - 1) * 2 * pi / 6
            Dij[:, i] .= [D * cos(theta), D * sin[theta], 0]
        end

        if type == "interfacial"
            for i in 1:6
                Dij[:, i] .= cross_product(Dij[:, i], [0, 0, 1.0])
            end
            @info("Interfacial DMI for TriangularMesh has been added!")
        else
            @info("Bulk DMI for TriangularMesh has been added!")
        end
    elseif isa(sim.mesh, CubicMesh)
        Dij[:, 1] .= [-D, 0, 0]
        Dij[:, 2] .= [D, 0, 0]
        Dij[:, 3] .= [0, -D, 0]
        Dij[:, 4] .= [0, D, 0]
        Dij[:, 5] .= [0, 0, -D]
        Dij[:, 6] .= [0, 0, D]

        if type == "interfacial"
            for i in 1:6
                Dij[:, i] .= cross_product(Dij[:, i], [0, 0, 1.0])
            end
            @info("Interfacial DMI for CubicMesh has been added!")
        else
            @info("Bulk DMI for CubicMesh has been added!")
        end
    end

    if type == "interfacial"
        Dij[3, :] .+= Dz
    end
    copyto!(sim.exch.D, Dij)
    return nothing
end

#Hx, Hy, Hz in energy unit， just as J and D
function add_zeeman(sim::MonteCarlo; Hx=0, Hy=0, Hz=0)
    zeeman = sim.zeeman
    zeeman.Hx = Hx / k_B
    zeeman.Hy = Hy / k_B
    zeeman.Hz = Hz / k_B
    return nothing
end

function update_zeeman(sim::MonteCarlo; Hx=0, Hy=0, Hz=0)
    add_zeeman(sim; Hx=Hx, Hy=Hy, Hz=Hz)
    return nothing
end

function add_anis(sim::MonteCarlo; Ku=1, Kc=0, axis=(0, 0, 1))
    anis = sim.anis
    anis.Ku = Ku / k_B
    anis.axis = normalize(axis)
    anis.Kc = Kc / k_B
    return nothing
end

"""
    add_anis_kagome(sim::MonteCarlo, Ku::Float64)

Add Anisotropy for kagome system, where the energy density is given by

```math
    E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
and u is one of ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0) and ax3=(-0.5,sqrt(3)/2,0).
"""
function add_anis_kagome(sim::MonteCarlo; Ku=0)
    T = Float[]

    sim.anis = KagomeAnisotropyMC(T(0))
    sim.anis.Ku = Ku / k_B

    return nothing
end


function add_anis_kagome_6fold(sim::MonteCarlo; K1=0, K2=0)
    T = Float[]

    sim.anis = KagomeAnisotropy6FoldMC(T(0), T(0))
    sim.anis.K1 = K1 / k_B
    sim.anis.K2 = K2 / k_B

    return nothing
end

function init_m0(sim::MonteCarlo, m0::TupleOrArrayOrFunction; norm=true)
    spin = zeros(Float[], 3 * sim.n_total)
    init_vector!(spin, sim.mesh, m0)
    if norm
        normalise(spin, sim.n_total)
    end
    shape = Array(sim.shape)
    for i in 1:(sim.n_total)
        if !shape[i]
            spin[3 * i - 2] = NaN32
            spin[3 * i - 1] = NaN32
            spin[3 * i] = NaN32
        end
    end
    copyto!(sim.spin, spin)
    return true
end

function update_ngbs(mesh, shape::Array{Bool})
    ngbs = Array(mesh.ngbs)
    for i in 1:(mesh.n_total)
        for j in 1:(mesh.n_ngbs)
            id = ngbs[j, i]
            if id > 0 && ((!shape[id]) || (!shape[i]))
                ngbs[j, i] = -1
            end
        end
    end
    copyto!(mesh.ngbs, ngbs)
    return nothing
end

function set_shape(sim::MonteCarlo, fun_Ms::Function)
    mesh = sim.mesh
    shape = ones(Bool, mesh.n_total)
    for k in 1:(mesh.nz), j in 1:(mesh.ny), i in 1:(mesh.nx)
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        shape[id] = fun_Ms(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    copyto!(sim.shape, shape)
    update_ngbs(sim.mesh, shape)
    return true
end

function set_shape_to_kagome(sim::MonteCarlo)
    mesh = sim.mesh
    shape = ones(Bool, mesh.n_total)
    for k in 1:(mesh.nz), j in 1:(mesh.ny), i in 1:(mesh.nx)
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        #shape[id] = true
        if i % 2 == 0 && j % 2 == 0
            shape[id] = false
        end
    end
    copyto!(sim.shape, shape)
    update_ngbs(sim.mesh, shape)
    return true
end

function run_step(sim::MonteCarlo)
    if sim.mc_2d
        uniform_random_circle_xy(sim.nextspin, sim.rnd, sim.n_total)
    else
        uniform_random_sphere(sim.nextspin, sim.rnd, sim.n_total)
    end

    run_step_bias(sim, 0)
    run_step_bias(sim, 1)
    run_step_bias(sim, 2)
    sim.steps += 1

    return nothing
end

function average_m(sim::MonteCarlo)
    shape = sim.shape
    b = reshape(sim.spin, 3, sim.n_total)
    return Tuple(sum(b; dims=2) ./ sum(shape))
end


function run_sim(sim::MonteCarlo; maxsteps=10000, save_m_every=10, save_vtk_every=-1,
                 save_ovf_every=-1, ovf_format="binary8")

    for i in 1:maxsteps
        if save_m_every > 0
            if sim.steps % save_m_every == 0
                energy = compute_system_energy(sim)
                @info @sprintf("step = %5d  total_energy=%g k_B", sim.steps, energy)
                sim.saver.nsteps += 1
                write_data(sim)
            end
        end

        if save_ovf_every > 0
            if sim.steps % save_ovf_every == 0
                save_ovf(sim, @sprintf("%s_%d", sim.name, sim.steps); dataformat=ovf_format)
            end
        end

        if save_vtk_every > 0
            if sim.steps % save_vtk_every == 0
                save_vtk(sim, @sprintf("%s_%d", sim.name, sim.steps))
            end
        end

        run_step(sim)
    end
end

function compute_clock_number(m::Array{T,1}, cn::Array{Float32,1}, shape::Array{Bool},
                              mesh::TriangularMesh) where {T<:AbstractFloat}
    # nx,ny = mesh.nx, mesh.ny
    # v = zeros(T, nx*ny)
    # compute_skyrmion_number(v, m, mesh)
    # return sum(v)
    signA = zeros(typeof(0.0), 6)
    ngbs = Array(mesh.ngbs)
    for i in 1:(mesh.n_total)
        mx = m[i * 3 - 2]
        my = m[i * 3 - 1]
        if shape[i]
            if ngbs[1, i] == -1
                #     -m×n2+m×n3     -m×n5+m×n6
                signA .= (0, -1, 1, 0, -1, 1)
                # cn[i]=1
            elseif ngbs[2, i] == -1
                # m×n1     -m×n3+m×n4     -m×n6
                signA .= (1, 0, -1, 1, 0, -1)
                # cn[i]=2
            elseif ngbs[3, i] == -1
                #-m×n1+m×n2     -m×n4+m×n5
                signA .= (-1, 1, 0, -1, 1, 0)
                # cn[i]=3
            end
            for j in 1:6
                id = ngbs[j, i]
                if id > 0
                    nx = m[id * 3 - 2]
                    ny = m[id * 3 - 1]
                    cn[i] += signA[j] * (mx * ny - my * nx)
                end
            end
            cn[i] = cn[i] * 2 / sqrt(3)
        else
            cn[i] = NaN32
        end#if
    end
end

function save_vtk_clocknum(sim::AbstractSim, cn::Array{Float32,1}, fname::String;
                           fields::Array{String,1}=String[])
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx, ny, nz)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    for k in 1:nz, j in 1:ny, i in 1:nx
        xyz[1, i, j, k] = (i - 0.5) * dx - (j - 1) * dx / 2
        xyz[2, i, j, k] = (j - 0.5) * dy
        xyz[3, i, j, k] = (k - 0.5) * dz
    end

    vtk = vtk_grid(fname, xyz)
    b = reshape(Array(spin), (3, nx, ny, nz))
    vtk_point_data(vtk, b, "m")
    vtk_point_data(vtk, reshape(cn, (nx, ny, nz)), "clocknum")
    if length(fields) > 0
        compute_fields_to_gpu(sim, sim.spin, 0.0)
        fields = Set(fields)
        for i in sim.interactions
            if i.name in fields
                b = reshape(i.field, (3, nx, ny, nz))
                vtk_point_data(vtk, b, i.name)
            end
        end
    end
    return vtk_save(vtk)
end
