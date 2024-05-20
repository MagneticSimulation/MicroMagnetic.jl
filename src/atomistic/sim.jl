
export set_mu_s, set_mu_s_kagome

"""
    set_mu_s(sim::AtomisticSim, Ms::NumberOrArrayOrFunction)

Set magnetic moment mu_s of the studied system. For example,

```julia
   set_mu_s(sim, 2*mu_B)
```
or
```julia
function circular_shape(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 2*mu_B
    end
    return 0.0
end
set_mu_s(sim, circular_shape)
```
"""
function set_mu_s(sim::AtomisticSim, init::NumberOrArrayOrFunction)
    T = Float[]
    Ms = zeros(T, sim.n_total)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.mu_s, Ms)
    return true
end

function set_mu_s_kagome(sim::AtomisticSim, Ms::Number)
    mesh = sim.mesh
    T = Float[]
    mu_s = zeros(T, sim.n_total)
    for k in 1:(mesh.nz), j in 1:(mesh.ny), i in 1:(mesh.nx)
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        mu_s[id] = Ms
        if i % 2 == 0 && j % 2 == 0
            mu_s[id] = 0.0
        end
    end
    copyto!(sim.mu_s, mu_s)
    return true
end

"""
    add_exch(sim::AtomisticSim, J::NumberOrArray; name="exch")

Add exchange energy to the system. The length of J should be equal to the length of neigbours.
"""
function add_exch(sim::AtomisticSim, J1::NumberOrArray; name="exch", J2=0, J3=0, J4=0)
    mesh = sim.mesh
    T = Float[]

    Js = create_zeros(mesh.n_ngbs)
    if isa(J1, Number)
        Js .= J1
    elseif length(J1) == length(Js)
        copyto!(Js, J1)
    else
        @error("J1 should be a number or an array with length $(mesh.n_ngbs)!")
    end

    Js2 = T[]
    if J2 != 0
        Js2 = create_zeros(mesh.n_ngbs2)
        if isa(J2, Number)
            Js2 .= J2
        elseif length(J2) == length(Js2)
            copyto!(Js2, J2)
        else
            @error("J2 should be a number or an array with length $(mesh.n_ngbs2)!")
        end
    end

    Js3 = T[]
    if J3 != 0
        Js3 = create_zeros(mesh.n_ngbs3)
        if isa(J3, Number)
            Js3 .= J3
        elseif length(J3) == length(Js3)
            copyto!(Js3, J3)
        else
            @error("J3 should be a number or an array with length $(mesh.n_ngbs3)!")
        end
    end

    Js4 = T[]
    if J4 != 0
        Js4 = create_zeros(mesh.n_ngbs4)
        if isa(J4, Number)
            Js4 .= J4
        elseif length(J4) == length(Js4)
            copyto!(Js4, J4)
        else
            @error("J4 should be a number or an array with length $(mesh.n_ngbs4)!")
        end
    end

    N = sim.n_total
    field = create_zeros(3 * N)
    energy = create_zeros(N)

    exch = HeisenbergExchange(Js, Js2, Js3, Js4, field, energy, name)
    push!(sim.interactions, exch)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return exch
end

@doc raw"""
    add_dmi(sim::AtomisticSim, D::Real; name="dmi", type="bulk")

Add bulk dmi energy to the system. The DMI is defined as
```math
\mathcal{H}_\mathrm{dmi} = \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
```
where $\mathbf{D}_{i j}$ is the DM vector. For the bulk dmi $\mathbf{D}_{i j} = D \hat{r}_{ij}$ and for interfacial dmi 
$\mathbf{D}_{i j} = D \hat{r}_{ij} \times \hat{z}$
"""
function add_dmi(sim::AtomisticSim, D::Real; name="dmi", type="bulk")
    N = sim.n_total
    mesh = sim.mesh
    T = Float[]
    field = create_zeros(3 * N)
    energy = create_zeros(N)

    Dij = zeros(T, (3, mesh.n_ngbs))

    if isa(sim.mesh, TriangularMesh)
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
            @info("Bulk DMI  for TriangularMesh has been added!")
        end

        Dij = kernel_array(Dij)
        dmi = HeisenbergDMI(Dij, field, energy, name)
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
        Dij = kernel_array(Dij)
        dmi = HeisenbergDMI(Dij, field, energy, name)
    elseif isa(sim.mesh, CylindricalTubeMesh)
        nr = sim.mesh.nr
        coords = sim.mesh.coordinates
        ngbs = Array(sim.mesh.ngbs)
        #Dij stores D*r_{ij}, i.e., r_{n_r, 1}, r12, r23, ..., r_{(n_r-1)n_r},  r_{n_r, 1}
        Dij = zeros(T, 3, nr + 1)
        for i in 1:nr
            j = ngbs[1, i]  # the left one
            rx = coords[1, i] - coords[1, j]
            ry = coords[2, i] - coords[2, j]
            rz = coords[3, i] - coords[3, j]
            r = sqrt(rx * rx + ry * ry + rz * rz)
            Dij[1, i] = D * rx / r
            Dij[2, i] = D * ry / r
            Dij[3, i] = D * rz / r
        end
        Dij[:, nr + 1] .= Dij[:, 1]
        Dij = kernel_array(Dij)
        dmi = HeisenbergTubeBulkDMI(T(D), Dij, field, energy, name)
        @info("Bulk DMI for CylindricalTubeMesh has been added!")
    end

    push!(sim.interactions, dmi)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return dmi
end

"""
    add_exch_kagome(sim::AtomisticSim, Jxy::Number, Jz::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch_kagome(sim::AtomisticSim, Jxy::Number, Jz::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(n_ngbs)

    if n_ngbs != 8
        error("The number of neigbours is not 8.")
    end

    Js[1:6] .= Jxy
    Js[7:8] .= Jz
    return add_exch(sim, Js; name=name)
end

"""
    add_anis_kagome(sim::AtomisticSim, Ku::Float64; ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0), ax3=(-0.5,sqrt(3)/2,0), name="anis")

Add Anisotropy for kagome system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis_kagome(sim::AtomisticSim, Ku::Float64; ax1=(-0.5, -sqrt(3) / 2, 0),
                         ax2=(1, 0, 0), ax3=(-0.5, sqrt(3) / 2, 0), name="anis")
    n_total = sim.n_total
    T = _cuda_using_double.x ? Float64 : Float32
    field = zeros(T, 3 * n_total)
    energy = zeros(T, n_total)
    anis = KagomeAnisotropy(T(Ku), ax1, ax2, ax3, field, energy, T(0.0), name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return anis
end

@doc raw"""
    add_anis_tube(sim::AtomisticSim, Ku::Float64; name="anis")

add anisotropy to the system when the tube mesh is used. The anisotropy axis $u$ 
is along with the radial direction.

```math
E_\mathrm{anis} = - K_{u} (\vec{m} \cdot \hat{u})^2
```
"""
function add_anis_tube(sim::AtomisticSim, Ku::Float64; name="anis")
    n_total = sim.n_total
    T = _cuda_using_double.x ? Float64 : Float32
    field = zeros(T, 3 * n_total)
    energy = zeros(T, n_total)
    axes = zeros(T, 3, n_total)

    nr = sim.mesh.nr
    for i in 1:n_total
        theta = 2 * pi * (i - 1) / nr
        axes[1, i] = cos(theta)
        axes[2, i] = sin(theta)
    end

    anis = TubeAnisotropy(T(Ku), CuArray(axes), field, energy, T(0.0), name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return anis
end

"""
    add_magnetoelectric_laser(sim::AtomisticSim, lambda::Float64, E::Float64, B::Float64, omega::Float64; delta=0, direction=001, name="lasers")

Add the interaction of high-frequency lasers to multiferroic insulator Cu2OSeO3. 
The Hamiltonian is given by 

```math
\\mathcal{H}_\\mathrm{laser} =  -\\sum_{i} \\mu_s \\mathbf{m}_i \\cdot \\mathbf{B}(t) - \\sum_{i} \\mathbf{P}_i \\cdot \\mathbf{E}(t)
```

The high-frequency laser is described as 

```math
\\mathbf{E}(t) =  E ( \\sin (\\omega t + \\delta), \\cos \\omega t, 0) \\qquad
\\mathbf{B}(t) =  B ( \\cos \\omega t, -\\sin(\\omega t + \\delta), 0)
```
where δ determines the laser polarization, i.e., δ = 0 for right-circularly polarized (RCP), 
δ=π/2 for linearly polarized and δ=π for left-circularly polarized (LCP).
"""
function add_magnetoelectric_laser(sim::AtomisticSim, lambda::Float64, E::Float64,
                                   B::Float64, omega::Float64; delta=0, direction=001,
                                   name="lasers")
    n_total = sim.n_total
    F = _cuda_using_double.x ? Float64 : Float32
    field = zeros(F, 3 * n_total)
    energy = zeros(F, n_total)

    laser = MagnetoelectricLaser(F(lambda), F(E), F(B), F(omega), F(delta), direction,
                                 field, energy, F(0.0), name)

    push!(sim.interactions, laser)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return laser
end

@doc raw"""
    add_demag(sim::AtomisticSim; name="demag", Nx=0, Ny=0, Nz=0 )

add dipolar interaction into the system.

```math
\mathcal{H}_{\mathrm{d}}=-\frac{\mu_0 \mu_s^2}{4 \pi} \sum_{i<j} \frac{3\left(\mathbf{m}_i \cdot \hat{\mathbf{r}}_{i j}\right)\left(\mathbf{m}_j \cdot \hat{\mathbf{r}}_{i j}\right)-\mathbf{m}_i \cdot \mathbf{m}_j}{r_{i j}^3}
```
"""
function add_demag(sim::AtomisticSim; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J",
                         o::AbstractSim -> sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return demag
end
