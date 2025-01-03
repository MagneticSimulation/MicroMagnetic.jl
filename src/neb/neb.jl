using LinearAlgebra
using Printf

export NEB

mutable struct NEB{T<:AbstractFloat} <: AbstractSim
    sim::AbstractSim # we need the sim instance to compute the effect fields and related functions
    driver::Driver # SD driver or LLG driver
    n_images::Int64 # number of free images
    n_total::Int64 # n_total = neb.sim.n_total * n_images
    dof::Int64 # dof = 3*sim.n_total
    clib_image::Int64
    nsteps::Int64 # number of steps
    image_l::AbstractArray{T,1} # left image
    image_r::AbstractArray{T,1} # right image
    spin::AbstractArray{T,1}
    prespin::AbstractArray{T,1}
    field::AbstractArray{T,1}
    tangent::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    energy_cpu::Array{T,1}
    distance::Array{T,1}
    name::String
    saver_energy::DataSaver
    saver_distance::DataSaver
    spring_constant::Float64
    pins::AbstractArray{Bool,1}
    NEB{T}() where {T<:AbstractFloat} = new()
end

"""
    NEB(sim::AbstractSim, init_images::TupleOrArray, frames_between_images::TupleOrArray; 
        name="NEB", spring_constant=1.0e5, driver="LLG", clib_image=-1)

Create a NEB instance.

# Arguments

- `sim::AbstractSim`: Simulation instance that includes the interactions.
- `init_images::TupleOrArray`: Tuple or array of initial and final state images. Intermediate state images can also be included.
- `frames_between_images::TupleOrArray`: Tuple or array specifying the number of frames between each pair of images.
- `name::String="NEB"`: Name for the NEB instance.
- `spring_constant::Float64=1.0e5`: Spring constant used in the NEB simulation.
- `driver::String="LLG"`: Driver for the NEB simulation. Options are `"SD"` or `"LLG"`.
- `clib_image::Int=-1`: Optional parameter for specifying a particular image in the library (default is -1).

# Returns

A NEB instance configured with the provided parameters.

# Example

```julia
# define mesh and the corresponding parameters
mesh = FDMesh(nx=60, ny=60, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")
params = Dict(
    :Ms => 3.84e5,
    :A => 3.25e-12,
    :D => 5.83e-4,
    :H => (0, 0, 120 * mT)
)

# Define the initial and final state, stored in the init_images list.
# Any acceptable object, such as a function, a tuple, or an array, can be used.
# The init_images list can also contain the intermediate state if you have one.
init_images = [read_vtk("skx.vts"), (0, 0, 1)]

# Define the interpolation array to specify the number of images used in the NEB simulation.
# The length of the interpolation array should be the length of init_images minus one.
# For example, if init_images = [read_vtk("skx.vts"), read_vtk("skx2.vts"), (0, 0, 1)], 
# the length of interpolation should be 2, i.e., something like interpolation = [5,5].
interpolation = [6]

# Use the create_sim method to create a Sim instance.
sim = create_sim(mesh; params...)

# Create the NEB instance and set the spring_constant. The driver can be "SD" or "LLG".
neb = NEB(sim, init_images, interpolation; name="skx_fm", driver="SD")
# neb.spring_constant = 1e7

# Relax the entire system.
relax(neb; stopping_dmdt=0.1, save_vtk_every=1000, max_steps=5000)
```
"""
function NEB(sim::AbstractSim, given_images::TupleOrArray,
             frames_between_images::TupleOrArray; name="NEB", spring_constant=1.0e5,
             driver="LLG", clib_image=-1)
    n_total = sim.mesh.nx * sim.mesh.ny * sim.mesh.nz
    n_images = sum(frames_between_images) + length(given_images) - 2

    if n_images <= 0
        error("The number of free images should larger than 1 (typically 10~20).")
    end

    F = Float[]
    neb = NEB{F}()
    neb.n_images = n_images
    neb.n_total = n_total * n_images
    neb.dof = 3 * sim.n_total
    neb.sim = sim
    neb.name = name
    neb.image_l = create_zeros(3 * sim.n_total)
    neb.image_r = create_zeros(3 * sim.n_total)
    neb.spin = create_zeros(3 * neb.n_total)
    neb.prespin = create_zeros(3 * neb.n_total)
    neb.field = create_zeros(3 * neb.n_total)
    neb.tangent = create_zeros(3 * neb.n_total)
    neb.energy = create_zeros(neb.n_images + 2)
    neb.energy_cpu = zeros(F, neb.n_images + 2)
    neb.distance = zeros(F, neb.n_images + 1)

    neb.spring_constant = spring_constant
    neb.clib_image = clib_image
    neb.nsteps = 0

    neb.pins = repeat(sim.pins, n_images)
    given_images_aligned = []
    for i in 1:length(given_images)
        init_m0(sim, given_images[i])
        push!(given_images_aligned, Array(sim.spin))
    end

    images = init_images(given_images_aligned, frames_between_images, sim.n_total)
    if size(images)[2] != n_images + 2
        msg = @sprintf("The number of init_images is not correct: Expected: %d images but got %d.",
                       n_images, size(images)[2])
        error(msg)
    end

    copyto!(neb.image_l, images[:, 1])
    copyto!(neb.image_r, images[:, end])
    copyto!(neb.spin, images[:, 2:(neb.n_images + 1)])
    copyto!(neb.prespin, neb.spin)

    compute_system_energy(sim, neb.image_l, 0.0)
    neb.energy_cpu[1] = sum(sim.energy)
    compute_system_energy(sim, neb.image_r, 0.0)
    neb.energy_cpu[end] = sum(sim.energy)

    effective_field_energy(neb, neb.spin)

    init_saver(neb)
    neb.driver = create_driver(driver, "DormandPrince", neb.n_total)
    if driver == "LLG"
        neb.driver.tol = 1e-5
        neb.driver.precession = false
        neb.driver.alpha = 0.2
    end

    return neb
end

function init_images(given_images::TupleOrArray, frames_between_images::TupleOrArray,
                     n_total)
    n_images = length(given_images) + sum(frames_between_images)
    images = zeros(3 * n_total, n_images)
    images[:, 1] .= given_images[1]
    image_count = 2
    for i in 1:(length(given_images) - 1)
        if frames_between_images[i] > 0
            for j in 1:frames_between_images[i]
                images[:, image_count] .= slerp(given_images[i], given_images[i + 1],
                                                j / (frames_between_images[i] + 1), n_total)
                image_count += 1
            end
            images[:, image_count] .= given_images[i + 1]
            image_count += 1
        end
    end
    return images
end

function init_saver(neb::NEB)
    step_item = SaverItem("step", "", o::NEB -> o.nsteps)
    neb.saver_energy = DataSaver(@sprintf("%s_energy.txt", neb.name), false, 0.0, 0,
                                 [step_item])
    neb.saver_distance = DataSaver(@sprintf("%s_distance.txt", neb.name), false, 0.0, 0,
                                   [step_item])
    N = neb.n_images
    for n in 1:(N + 2)
        name = @sprintf("E_total_%g", n)
        item = SaverItem(name, "J", o::NEB -> o.energy_cpu[n])
        push!(neb.saver_energy.items, item)
    end

    for n in 1:(N + 1)  #TODO: using a single function?
        name = @sprintf("distance_%d", n)
        item = SaverItem(name, "", o::NEB -> o.distance[n])
        push!(neb.saver_distance.items, item)
    end
end

function effective_field(neb::NEB, spin::AbstractArray, t::Float64)
    dof = neb.dof

    effective_field_energy(neb, spin)
    compute_tangents(neb, spin)
    compute_distance(neb, spin)

    for n in 1:(neb.n_images)
        f = view(neb.field, ((n - 1) * dof + 1):(n * dof))
        t = view(neb.tangent, ((n - 1) * dof + 1):(n * dof))

        if n == neb.clib_image
            ft = -2 * LinearAlgebra.dot(f, t)
        else
            ft = neb.spring_constant * (neb.distance[n + 1] - neb.distance[n]) -
                 LinearAlgebra.dot(f, t)
        end
        f .+= ft .* t
    end
    return nothing
end

# compute the micromagnetic effective field and energy
function effective_field_energy(neb::NEB, spin::AbstractArray)
    sim = neb.sim
    dof = neb.dof

    for n in 1:(neb.n_images)
        f = view(neb.field, ((n - 1) * dof + 1):(n * dof))
        m = view(neb.spin, ((n - 1) * dof + 1):(n * dof))
        effective_field_energy(sim, m, 0.0) # will compute the effective field and energy
        f .= sim.field
        neb.energy_cpu[n + 1] = sum(sim.energy)
    end
    return copyto!(neb.energy, neb.energy_cpu)
end

function compute_distance(neb::NEB, spin::AbstractArray)
    n_total = neb.sim.n_total
    N = neb.n_images
    dof = neb.dof
    ds = neb.sim.energy #we borrow the sim.energy

    kernel! = compute_distance_kernel!(default_backend[], groupsize[])

    for n in 0:N
        m1 = n == 0 ? neb.image_l : view(neb.spin, ((n - 1) * dof + 1):(n * dof))
        m2 = n == N ? neb.image_r : view(neb.spin, (n * dof + 1):((n + 1) * dof))

        kernel!(m1, m2, ds; ndrange=n_total)
        neb.distance[n + 1] = LinearAlgebra.norm(ds)
    end

    return nothing
end

function compute_tangents(neb::NEB, spin::AbstractArray)
    N = neb.n_images

    kernel! = compute_tangents_kernel!(default_backend[], groupsize[])
    kernel!(neb.tangent, spin, neb.image_l, neb.image_r, neb.energy, N, neb.dof;
            ndrange=neb.dof)

    kernel! = reduce_tangent_kernel!(default_backend[], groupsize[])
    kernel!(neb.tangent, spin; ndrange=neb.n_total)

    dof = neb.dof
    for n in 1:N
        t = view(neb.tangent, (dof * (n - 1) + 1):(dof * n))
        norm_t = LinearAlgebra.norm(t)
        t .= t / norm_t
    end

    return nothing
end

function relax(sim::NEB; max_steps=10000, stopping_dmdt=0.01, save_data_every=1,
               save_vtk_every=-1, using_time_factor=true, vtk_folder="vtks")

    # to dertermine which driver is used.
    llg_driver = isa(sim.driver, LLG)

    time_factor = using_time_factor ? 2.21e5 / 2 : 1.0
    dmdt_factor = using_time_factor ? (2 * pi / 360) * 1e9 : 1

    if save_vtk_every > 0 && !isdir(vtk_folder)
        mkdir(vtk_folder)
    end

    if save_vtk_every > 0
        save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, 0)))
    end

    N_spins = sim.n_total
    dm = create_zeros(3 * N_spins)

    driver = sim.driver
    @info @sprintf("Running Driver : %s.", typeof(driver))
    for i in 0:max_steps
        @timeit timer "run_step" run_step(sim, driver)

        sim.nsteps += 1

        step_size = llg_driver ? driver.integrator.step : driver.tau / time_factor

        compute_dm!(dm, sim.prespin, sim.spin, N_spins)
        max_dmdt = maximum(dm) / step_size

        t = llg_driver ? sim.driver.integrator.t : 0.0
        if llg_driver
            @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e", i,
                           step_size, t, max_dmdt / dmdt_factor)
        else
            @info @sprintf("step =%5d  step_size=%10.6e  max_dmdt=%10.6e", i, step_size,
                           max_dmdt / dmdt_factor)
        end

        if save_data_every > 0 && i % save_data_every == 0
            effective_field_energy(sim, sim.spin)
            write_data(sim, sim.saver_energy)
            write_data(sim, sim.saver_distance)
        end

        if save_vtk_every > 0 && i % save_vtk_every == 0
            save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
        end

        if max_dmdt < stopping_dmdt * dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
            if save_data_every > 0 || save_data_every == -1
                effective_field_energy(sim, sim.spin)
                write_data(sim, sim.saver_energy)
                write_data(sim, sim.saver_distance)
            end

            if save_vtk_every > 0
                save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)))
            end

            break
        end
    end

    return nothing
end

#TODO: merge save_ovf and save_vtk???
function save_ovf(neb::NEB, fname::String; type::DataType=Float64)
    sim = neb.sim
    dof = neb.dof

    id_start = 0
    sim.spin .= neb.image_l
    save_ovf(sim, @sprintf("%s_0", fname); type=type)

    N = neb.n_images
    for n in 1:N
        b = view(neb.spin, ((n - 1) * dof + 1):(n * dof))
        sim.spin .= b
        save_ovf(sim, @sprintf("%s_%d", fname, n); type=type)
    end

    sim.spin .= neb.image_r
    return save_ovf(sim, @sprintf("%s_%d", fname, N + 1); type=type)
end

function save_vtk(neb::NEB, fname::String)
    sim = neb.sim
    dof = neb.dof
    N = neb.n_images

    sim.spin .= neb.image_l
    save_vtk(sim, @sprintf("%s_0", fname))
    for n in 1:N
        b = view(neb.spin, ((n - 1) * dof + 1):(n * dof))
        sim.spin .= b
        save_vtk(sim, @sprintf("%s_%d", fname, n))
    end
    sim.spin .= neb.image_r
    return save_vtk(sim, @sprintf("%s_%d", fname, N + 1))
end
