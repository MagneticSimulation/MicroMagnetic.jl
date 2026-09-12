# Adaptive Mesh Refinement for micromagnetics simulations, following
# C. J. García-Cervera and A. M. Roma, "Adaptive Mesh Refinement for
# Micromagnetics Simulations", IEEE Trans. Magn. 42(6), 1648 (2006).
#
# The sample is covered by composite (block-structured) grids: a base uniform
# grid plus nested rectangular patches at a refinement ratio of two, which
# follow sharp transitions of the magnetization (walls, vortices) and the
# domain boundary. All levels advance with the same time step using the
# Gauss-Seidel projection method (GPSM) integrator, in which the exchange term
# is treated implicitly through a constant sparse Cholesky factorization per
# region (re-built only at remeshing), while anisotropy/Zeeman are pointwise
# and the stray field is assembled from a base-grid FFT plus local correction
# boxes around the patches (amr/demag.jl, demag v3).
#
# v1 simplifications relative to the paper
#   * demag uses per-level full-box FFT solves with delta-magnetization
#     corrections instead of a composite multigrid solver;
#   * ghost cells are always filled by linear interpolation from the next
#     coarser composite level (no same-level copy at patch seams);
#   * open boundaries only, uniform Ms, CPU backend, exchange + uniaxial
#     anisotropy + static Zeeman + demag energy terms.

"""
    AMRRegion

One integration region of the composite grid: the base level (`level == 1`,
covering the whole domain) or a rectangular patch (`level >= 2`).

A patch owns a `MicroSim` built on a *ghosted* `FDMesh`: the patch interior is
embedded in a box with `ghost` layers of ghost cells on every side, so the
package's exchange kernels and `build_exch_matrix` work verbatim; the ghost
values are refreshed from the next coarser composite level at every step. The
physical-domain neighbours of a patch at the sample edge are filled by
constant extension, which reproduces the free (ngbs = -1) boundary used by
uniform sims.
"""
mutable struct AMRRegion{T<:AbstractFloat}
    level::Int                        # 1 = base level, >= 2 = patch level
    sim::MicroSim{T}                  # light sim: mesh + Exchange interaction
    # interior cell range in the coordinates of this level (base: 1:n)
    ir::UnitRange{Int}
    jr::UnitRange{Int}
    kr::UnitRange{Int}
    # GPSM per-region context (exchange treated implicitly). Cells that are
    # NOT authoritative for this region (ghost ring, cells under finer
    # patches) are Dirichlet: their auxiliary field g is prescribed from the
    # composite box Gc[level] and injected into the neighboring equations.
    G::Any                            # cholesky(I - dtg*L_dirichlet), CPU sparse
    g1::AbstractArray{T,1}
    g2::AbstractArray{T,1}
    g3::AbstractArray{T,1}
    rhs::AbstractArray{T,1}
    hs_a::AbstractArray{T,1}          # anisotropy field buffer (3N, whole box)
    hs_z::AbstractArray{T,1}          # Zeeman field buffer (3N, whole box)
    prespin::AbstractArray{T,1}       # state before the step (for max dm/dt)
    dm::AbstractArray{T,1}            # per-cell |dm| scratch
    isdir::Vector{Bool}               # box cells with Dirichlet g
    cx::T                             # dtg * exchange coupling per x bond
    cy::T                             # ... per y bond
    cz::T                             # ... per z bond
end

"""
    DemagBox

One local correction box at level `level`: the bounding rectangle of a patch
cluster expanded by `buffer` cells (in level-`level` coordinates), carrying a
lightweight sim whose FDMesh has exactly the box dimensions and a single
Demag interaction. `patches` lists the indices (into `amr.patches[level-1]`)
of the patches whose Δ sources live in this box. See `src/amr/demag.jl` for
the demag-v3 scheme.
"""
mutable struct DemagBox{T<:AbstractFloat}
    level::Int
    ir::UnitRange{Int}
    jr::UnitRange{Int}
    kr::UnitRange{Int}
    sim::MicroSim{T}
    patches::Vector{Int}
end

"""
    AMRSim{T}

Composite-grid (block-structured AMR) micromagnetics simulation following
García-Cervera & Roma (2006). See `AMRSim(::FDMesh; ...)` for construction.
"""
mutable struct AMRSim{T<:AbstractFloat}
    base_mesh::FDMesh
    levels::Int                       # number of levels (>= 1)
    ghost::Int                        # ghost layers around patches
    refined::NTuple{3,Bool}           # which dimensions are refined (n0 > 1)
    dims::Vector{Tuple{Int,Int,Int}}  # global grid dims per level (1..levels)
    base::AMRRegion{T}
    patches::Vector{Vector{AMRRegion{T}}}   # patches[level-1] for level 2..levels
    C::Vector{AbstractArray{T,1}}     # composite magnetization box per level
    Gc::Vector{AbstractArray{T,1}}    # composite GPSM auxiliary field per level
    Phi::Vector{AbstractArray{T,1}}   # assembled demag field box per level
    shadow::Vector{MicroSim{T}}       # per-level full-box sim (replicated-parent FFT)
    scratch_a::Vector{AbstractArray{T,1}}   # per-level scratch
    scratch_b::Vector{AbstractArray{T,1}}   # per-level scratch
    dboxes::Vector{Vector{DemagBox{T}}}     # correction boxes per level 2..levels
    box_buffer::Int                   # correction-box buffer width (level-ℓ cells)
    auth::Vector{Vector{Bool}}        # auth[level]: cell is covered by a patch
    covered::Vector{Vector{Bool}}     # covered[level]: cell lies under finer patches
    # material / driver parameters (uniform in v1)
    Ms::T
    A::T
    Ku::T
    axis::Tuple{T,T,T}
    H0::Union{Nothing,Tuple{T,T,T}}
    alpha::T
    gamma::T
    precession::Bool
    demag::Bool
    dt::Float64
    time::Float64
    nsteps::Int
    maxdmdt::Float64
    g_initialized::Bool
    # remesh controls
    remesh_interval::Int
    refine_threshold::T
    grid_efficiency::T
    refine_boundary::Bool
    boundary_layers::Int
    max_coverage::T
    # bookkeeping
    name::String
    saver::DataSaver
    save_data::Bool
    n_remesh::Int
    ecache::Any                       # cached amr_energies for the saver (nothing = stale)
end

# ---------------------------------------------------------------- helpers

"""Dims of level `ell` (`refined` dims double per level)."""
_level_dims(base_mesh::FDMesh, refined::NTuple{3,Bool}, ell::Int) =
    (refined[1] ? base_mesh.nx << (ell - 1) : base_mesh.nx,
     refined[2] ? base_mesh.ny << (ell - 1) : base_mesh.ny,
     refined[3] ? base_mesh.nz << (ell - 1) : base_mesh.nz)

"""Cell size at level `ell` (`refined` dims halve per level)."""
_level_h(base_mesh::FDMesh, refined::NTuple{3,Bool}, ell::Int) =
    (refined[1] ? base_mesh.dx / 2^(ell - 1) : base_mesh.dx,
     refined[2] ? base_mesh.dy / 2^(ell - 1) : base_mesh.dy,
     refined[3] ? base_mesh.dz / 2^(ell - 1) : base_mesh.dz)

@inline _cell_index(i::Int, j::Int, k::Int, nx::Int, ny::Int) =
    (k - 1) * nx * ny + (j - 1) * nx + i

"""
    _region_sim(mesh, Ms, name) -> MicroSim

A minimal `MicroSim` without the `Sim` factory side effects (`@info`, saver
files, server state): only the array fields plus the parameter caches are
initialised, so `set_Ms`, `make_param` and `build_exch_matrix` work on it.
"""
function _region_sim(mesh::FDMesh, Ms::Number, name::String)
    T = Float[]
    sim = MicroSim{T}()
    sim.time = 0.0
    sim.name = name
    sim.mesh = mesh
    sim.n_total = mesh.n_total
    sim.spin = create_zeros(3 * mesh.n_total)
    sim.prespin = create_zeros(3 * mesh.n_total)
    sim.field = create_zeros(3 * mesh.n_total)
    sim.energy = create_zeros(mesh.n_total)
    sim.pins = Fill(false, mesh.n_total)
    sim.mu0_Ms = Fill(T(0), mesh.n_total)
    sim.mat_class = nothing
    sim.n_classes = 0
    sim.mat_class_layout = -1
    sim.inv_ms = Fill(T(0), mesh.n_total)
    sim.driver_name = "None"
    sim.driver = EmptyDriver()
    sim.interactions = []
    sim.save_data = false
    sim.saver = DataSaver(name, true, 0.0, 0, [])   # header_saved: never written
    set_Ms(sim, Ms)
    return sim
end

"""Exchange interaction with uniform stiffness on a region sim."""
function _add_region_exch!(sim::MicroSim, A::Real)
    T = eltype(sim.spin)
    n = sim.n_total
    A_kb = make_param(T, A, sim.mesh, n)
    exch = Exchange(A_kb, A_kb, A_kb, create_zeros(3 * n), create_zeros(n),
                    "exch", nothing, nothing, nothing, -1)
    push!(sim.interactions, exch)
    return exch
end

"""Demag interaction (open boundaries) on a region sim."""
function _add_region_demag!(sim::MicroSim)
    demag = init_demag(sim, 0, 0, 0)
    push!(sim.interactions, demag)
    return demag
end

function _region_buffers(::Type{T}, n::Int) where {T<:AbstractFloat}
    return (create_zeros(T, n), create_zeros(T, n), create_zeros(T, n),
            create_zeros(T, n))
end

"""Build the GPSM context of a region: the exchange Laplacian pruned to the
region's authoritative cells (ghost ring and cells under finer patches become
Dirichlet rows whose g is prescribed from the composite box `Gc[level]`), the
constant Cholesky factorization of  I - dtg*L  (rebuilt at remesh only), and
the per-bond coupling scalars used to inject the Dirichlet g into the right
hand side."""
function _build_gpsm!(amr::AMRSim{T}, r::AMRRegion{T}) where {T<:AbstractFloat}
    idx = findfirst(x -> isa(x, Exchange), r.sim.interactions)
    idx === nothing && error("GPSM region requires an Exchange interaction")
    L = build_exch_matrix(r.sim.interactions[idx], r.sim)
    dtg = amr.dt * amr.gamma / (1 + amr.alpha^2)

    # Dirichlet mask over the region box
    n = r.sim.n_total
    isdir = zeros(Bool, n)
    ngx, ngy, ngz = r.sim.mesh.nx, r.sim.mesh.ny, r.sim.mesh.nz
    if r.level == 1
        isdir .= amr.covered[1]     # identity mapping on the base box
    else
        gx, gy, gz = _ghosts(amr)
        lbx, lby = amr.dims[r.level][1], amr.dims[r.level][2]
        cov = amr.covered[r.level]
        for c in 1:ngz, b in 1:ngy, a in 1:ngx
            I = _cell_index(a, b, c, ngx, ngy)
            outside = !(gx < a <= ngx - gx && gy < b <= ngy - gy && gz < c <= ngz - gz)
            if outside
                isdir[I] = true
            else
                Ig = _cell_index(first(r.ir) - 1 + a - gx, first(r.jr) - 1 + b - gy,
                                 first(r.kr) - 1 + c - gz, lbx, lby)
                cov[Ig] && (isdir[I] = true)
            end
        end
    end

    # prune the Dirichlet rows AND columns (symmetric): G[d,d] = 1 then returns
    # the prescribed g at those cells, and the pruned interior block stays SPD
    didx = findall(isdir)
    if !isempty(didx)
        L[didx, :] .= 0
        L[:, didx] .= 0
        dropzeros!(L)
    end

    r.isdir = isdir
    # uniform exchange coupling per bond: dtg * A * (2/h^2) / (mu0 * Ms)
    hx, hy, hz = _level_h(amr.base_mesh, amr.refined, r.level)
    mscale = 1 / (mu_0 * amr.Ms)
    r.cx = T(dtg * amr.A * 2 / hx^2 * mscale)
    r.cy = T(dtg * amr.A * 2 / hy^2 * mscale)
    r.cz = T(dtg * amr.A * 2 / hz^2 * mscale)
    r.G = cholesky(I - dtg * L)
    n = r.sim.n_total
    r.g1, r.g2, r.g3, r.rhs = _region_buffers(T, n)
    r.hs_a = create_zeros(3 * n)
    r.hs_z = create_zeros(3 * n)
    r.prespin = create_zeros(3 * n)
    r.dm = create_zeros(T, n)
    return r
end

# ---------------------------------------------------------------- patches

"""Create a patch region at level `ell` covering `ir/jr/kr` (level-`ell`
coords), initialise its interior from the current composite state and build
its GPSM context."""
function _make_patch(amr::AMRSim{T}, ell::Int, ir, jr, kr, old_auth) where {T<:AbstractFloat}
    gx, gy, gz = _ghosts(amr)
    hx, hy, hz = _level_h(amr.base_mesh, amr.refined, ell)
    nxi, nyi, nzi = length(ir), length(jr), length(kr)
    ngx, ngy, ngz = nxi + 2gx, nyi + 2gy, nzi + 2gz
    # origin of the ghosted box in level-ell coordinates (1-based)
    i0 = first(ir) - gx
    j0 = first(jr) - gy
    k0 = first(kr) - gz
    mesh = FDMesh(dx=hx, dy=hy, dz=hz, nx=ngx, ny=ngy, nz=ngz,
                  x0=amr.base_mesh.x0 + (i0 - 1) * hx,
                  y0=amr.base_mesh.y0 + (j0 - 1) * hy,
                  z0=amr.base_mesh.z0 + (k0 - 1) * hz)
    sim = _region_sim(mesh, amr.Ms, "$(amr.name)_L$(ell)")
    _add_region_exch!(sim, amr.A)
    r = AMRRegion{T}(ell, sim, ir, jr, kr, nothing,
                     create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                     create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                     create_zeros(T, 1), create_zeros(T, 1),
                     zeros(Bool, 1), T(0), T(0), T(0))
    _init_patch_interior!(amr, r, old_auth)
    _build_gpsm!(amr, r)
    return r
end

# ---------------------------------------------------------------- constructor

"""
    AMRSim(mesh::FDMesh; levels=2, Ms, A, Ku=0.0, axis=(0,0,1), H0=nothing,
           demag=true, alpha=0.1, gamma=2.21e5, precession=true, dt=1e-12,
           remesh_interval=50, refine_threshold=0.1, grid_efficiency=0.7,
           refine_boundary=true, boundary_layers=2, ghost_width=2,
           max_coverage=0.85, name="amr", save_data=true)

Create an adaptive-mesh-refinement (composite grid) simulation following
García-Cervera & Roma, *IEEE Trans. Magn.* **42**, 1648 (2006).

The base uniform grid is `mesh`; `levels` is the total number of grid levels
(the base level plus `levels - 1` refinement levels, refinement ratio 2 per
refined dimension; dimensions with a single cell are never refined, so a
2D `nz=1` mesh refines only in-plane). The composite grid is regenerated
every `remesh_interval` steps from flags where the magnetization divergence
exceeds `refine_threshold` times its global maximum and, when
`refine_boundary`, within `boundary_layers` cells of the sample boundary
(Sections III-A/III-B of the paper). Flag cells are clustered into
rectangular patches by the Berger-Rigoutsos algorithm with grid efficiency
`grid_efficiency`; the finest-level coverage is capped at `max_coverage`.

All levels advance with the same time step `dt` (seconds) using the
Gauss-Seidel projection method, so the scheme is unconditionally stable and
`dt` on the order of a picosecond can be used (Section II of the paper).
Material parameters `Ms` (A/m), `A` (J/m, exchange; required), `Ku` (J/m³,
uniaxial anisotropy along `axis`) and the static field `H0` (A/m) are uniform
in this version. Call `init_m0(amr, m0)` before `relax`/`run_sim`.

!!! note "v1 limitations"
    CPU backend and Float64 only; uniform `Ms`; open boundaries; demag is
    treated by per-level FFT solves with delta corrections (not composite
    multigrid); ghost cells are linear interpolations from the coarser
    composite level.
"""
function AMRSim(mesh::FDMesh; levels::Int=2, Ms::Number, A::Real, Ku::Number=0.0,
                axis::Tuple=(0, 0, 1), H0=nothing, demag::Bool=true,
                alpha::Real=0.1, gamma::Real=2.21e5, precession::Bool=true,
                dt::Real=1e-12, remesh_interval::Int=50, refine_threshold::Real=0.1,
                grid_efficiency::Real=0.7, refine_boundary::Bool=true,
                boundary_layers::Int=2, ghost_width::Int=2, max_coverage::Real=0.85,
                box_buffer::Int=8, name::String="amr", save_data::Bool=true)
    Float[] === Float64 ||
        error("AMRSim v1 requires Float64 precision (set_precision(Float64)); " *
              "the GPSM Cholesky solver needs CHOLMOD/Float64.")
    default_backend[] isa KernelAbstractions.CPU ||
        error("AMRSim v1 supports the CPU backend only; call set_backend(\"cpu\").")
    levels >= 1 || error("levels must be >= 1")
    ghost_width >= 1 || error("ghost_width must be >= 1")
    remesh_interval >= 1 || error("remesh_interval must be >= 1")
    T = Float64
    refined = (mesh.nx > 1, mesh.ny > 1, mesh.nz > 1)
    # normalized easy axis
    alt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
    alt > 0 || error("the anisotropy axis must be nonzero")
    axis_n = (axis[1] / alt, axis[2] / alt, axis[3] / alt)

    amr = AMRSim{T}(mesh, levels, ghost_width, refined,
                    [_level_dims(mesh, refined, l) for l in 1:levels],
                    _region(T), Vector{Vector{AMRRegion{T}}}(),
                    Vector{AbstractArray{T,1}}(), Vector{AbstractArray{T,1}}(),
                    Vector{AbstractArray{T,1}}(), Vector{MicroSim{T}}(),
                    Vector{AbstractArray{T,1}}(), Vector{AbstractArray{T,1}}(),
                    Vector{Vector{DemagBox{T}}}(), box_buffer,
                    Vector{Vector{Bool}}(),
                    Vector{Vector{Bool}}(),
                    T(Ms), T(A), T(Ku),
                    (T(axis_n[1]), T(axis_n[2]), T(axis_n[3])),
                    H0 === nothing ? nothing : Tuple(T.(H0)),
                    T(alpha), T(gamma), precession, demag, Float64(dt), 0.0, 0, 0.0,
                    false,
                    remesh_interval, T(refine_threshold), T(grid_efficiency),
                    refine_boundary, boundary_layers, T(max_coverage),
                    name, DataSaver("", true, 0.0, 0, []), save_data, 0, nothing)

    # base region: whole domain, exchange (+ demag) interactions
    bx, by, bz = amr.dims[1]
    base_mesh = mesh
    base = _region_sim(base_mesh, Ms, name * "_L1")
    _add_region_exch!(base, A)
    amr.demag && _add_region_demag!(base)
    amr.base = AMRRegion{T}(1, base, 1:bx, 1:by, 1:bz, nothing,
                            create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                            create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                            create_zeros(T, 1), create_zeros(T, 1),
                            zeros(Bool, bx * by * bz), T(0), T(0), T(0))

    # per-level work arrays: C[1] aliases the base spin (the level-1
    # authoritative state), Phi[1] aliases the base Demag field buffer
    push!(amr.C, base.spin)
    push!(amr.Gc, create_zeros(3 * mesh.n_total))
    if demag
        di = findfirst(x -> isa(x, Demag), base.interactions)
        demag_field = base.interactions[di].field
    else
        demag_field = create_zeros(3 * mesh.n_total)
    end
    push!(amr.Phi, demag_field)
    for l in 2:levels
        nx, ny, nz = amr.dims[l]
        push!(amr.C, create_zeros(3 * nx * ny * nz))
        push!(amr.Gc, create_zeros(3 * nx * ny * nz))
        push!(amr.Phi, create_zeros(3 * nx * ny * nz))
        push!(amr.scratch_a, create_zeros(3 * nx * ny * nz))
        push!(amr.scratch_b, create_zeros(3 * nx * ny * nz))
        # per-level shadow sim: the replicated-parent composite source's
        # full-level FFT (the level-(l-1) context field at h_l resolution)
        hx, hy, hz = _level_h(mesh, refined, l)
        shmesh = FDMesh(dx=hx, dy=hy, dz=hz, nx=nx, ny=ny, nz=nz,
                        x0=mesh.x0, y0=mesh.y0, z0=mesh.z0)
        sh = _region_sim(shmesh, Ms, "$(name)_rep_L$l")
        _add_region_demag!(sh)
        push!(amr.shadow, sh)
    end
    empty!(amr.patches)
    for l in 2:levels
        push!(amr.patches, Vector{AMRRegion{T}}())
        push!(amr.dboxes, Vector{DemagBox{T}}())
    end
    for l in 1:levels
        n = prod(amr.dims[l])
        push!(amr.auth, ones(Bool, n))
        push!(amr.covered, zeros(Bool, n))
    end
    _build_gpsm!(amr, amr.base)

    if save_data
        _init_amr_saver!(amr)
    end
    @info "AMRSim created: $(amr.dims[1]) base grid, $levels levels, " *
          "refinement ratio 2, ghost $ghost_width"
    return amr
end

function _region(::Type{T}) where {T<:AbstractFloat}
    # placeholder replaced right after in the constructor; keeps the struct
    # instantiation type-stable
    return AMRRegion{T}(1, MicroSim{T}(), 1:1, 1:1, 1:1, nothing,
                        create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                        create_zeros(T, 1), create_zeros(T, 1), create_zeros(T, 1),
                        create_zeros(T, 1), create_zeros(T, 1),
                        zeros(Bool, 1), T(0), T(0), T(0))
end

"""All integration regions: the base level followed by every patch."""
regions(amr::AMRSim) = (amr.base, (p for ps in amr.patches for p in ps)...)

"""Number of composite-grid cells (sum over levels of authoritative cells)."""
function composite_ncells(amr::AMRSim)
    n = sum(!amr.covered[1][i] for i in 1:length(amr.covered[1]); init=0)
    for l in 2:amr.levels
        n += sum(!amr.covered[l][i] for i in 1:length(amr.covered[l]); init=0)
    end
    return n
end
