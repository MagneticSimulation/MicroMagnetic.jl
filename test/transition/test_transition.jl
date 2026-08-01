using LinearAlgebra
using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

function transition_macrospin()
    mesh = CubicMesh(; nx=1, ny=1, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="transition_macrospin", save_data=false)
    set_mu_s(sim, mu_B)
    add_anis(sim, 1.0meV; axis=(1, 0, 0), name="easy_x")
    add_anis(sim, -0.5meV; axis=(0, 0, 1), name="hard_z")
    init_m0(sim, (1, 0, 0))
    return sim
end

function transition_affine_atomistic()
    mesh = CubicMesh(; nx=3, ny=3, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="transition_affine",
              save_data=false)
    set_mu_s(sim, mu_B)
    add_exch(sim, 1.0meV)
    add_dmi(sim, 0.45meV)
    add_zeeman(sim, (0, 0, 2.8))
    pins = fill(false, 9)
    pins[1] = true
    set_pinning(sim, pins)
    init_m0(sim, (0, 0, 1))
    return sim
end

function transition_symmetric_atomistic()
    mesh = CubicMesh(; nx=3, ny=3, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="transition_symmetric",
              save_data=false)
    set_mu_s(sim, mu_B)
    add_exch(sim, 1.0meV)
    add_zeeman(sim, (0, 0, 1.0))
    init_m0(sim, (0, 0, 1))
    return sim
end

function test_transition_hessian()
    modes = compute_hessian_modes(
        transition_macrospin, [0.0, 1.0, 0.0];
        n_modes=2, finite_difference_step=1e-4,
        krylov_dimension=3, maximum_restarts=3,
        residual_tolerance=1e-27)
    @test modes.eigenvalues[1] < 0
    @test modes.eigenvalues[2] > 0
    @test all(isfinite, modes.eigenvalues)

    backend = MicroMagnetic._TransitionBackend(
        transition_affine_atomistic)
    @test !backend.active_cpu[1]
    spin = MicroMagnetic._device_vector(repeat([0.0, 0.0, 1.0], 9))
    tangent = MicroMagnetic._device_vector(
        collect(range(-0.7, 0.7; length=27)))
    context = MicroMagnetic._prepare_transition_hessian(backend, spin)
    @test !isnothing(context)
    analytic = MicroMagnetic.create_zeros(27)
    finite_difference = MicroMagnetic.create_zeros(27)
    MicroMagnetic._apply_transition_hessian!(
        analytic, backend, spin, tangent, 1e-5, context)
    MicroMagnetic._apply_transition_hessian!(
        finite_difference, backend, spin, tangent, 1e-5, nothing)
    @test Array(analytic) ≈ Array(finite_difference) rtol=1e-7 atol=1e-28
    @test all(iszero, Array(analytic)[1:3])

    symmetry_backend = MicroMagnetic._TransitionBackend(
        transition_symmetric_atomistic)
    symmetric_spin = repeat([0.0, 0.0, 1.0], 9)
    symmetries = MicroMagnetic._detect_transition_symmetries(
        symmetry_backend, (symmetric_spin, symmetric_spin))
    @test length(symmetries) == 1
    perturbed = copy(symmetric_spin)
    perturbed[4:6] .= [0.1, 0.0, sqrt(0.99)]
    projected = MicroMagnetic._device_vector(perturbed)
    source = copy(projected)
    MicroMagnetic._project_transition_spin!(
        projected, source, first(symmetries), symmetry_backend)
    @test MicroMagnetic._transition_symmetry_error(
        projected, first(symmetries), symmetry_backend.active_cpu) < 1e-12
end

function test_transition_saddle_and_path()
    saddle = find_saddle(
        transition_macrospin, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0];
        direction=1, maximum_steps=120, step_size=0.06,
        mode_update_every=3,
        gradient_tolerance=1e-27, eigenvalue_tolerance=1e-27,
        transverse_relaxation=0.0, n_modes=2,
        finite_difference_step=1e-4, krylov_dimension=3,
        maximum_restarts=3, residual_tolerance=1e-27)
    @test saddle.converged
    @test saddle.first_order
    @test saddle.negative_mode_count == 1
    @test abs(saddle.spin[1]) < 1e-4

    backend = MicroMagnetic._TransitionBackend(transition_macrospin)
    same = MicroMagnetic._safe_interpolate(
        backend, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5)
    opposite = MicroMagnetic._safe_interpolate(
        backend, [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5)
    @test Array(same) ≈ [1.0, 0.0, 0.0]
    @test norm(Array(opposite)) ≈ 1.0
    @test abs(dot(Array(opposite), [1.0, 0.0, 0.0])) < 1e-10
    interpolated_path = MicroMagnetic._safe_interpolated_path(
        backend, [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 7)
    @test interpolated_path[:, 1] ≈ [1.0, 0.0, 0.0]
    @test interpolated_path[:, end] ≈ [-1.0, 0.0, 0.0]
    @test maximum(abs.(vec(sum(interpolated_path.^2; dims=1)) .- 1)) <
          1e-12

    path = hcat([1.0, 0.0, 0.0], saddle.spin, [-1.0, 0.0, 0.0])
    band = GNEB(
        transition_macrospin, path; name="",
        spring_constant=1.0meV, stationary_images=[2])
    @test gneb_images(band) ≈ path
    @test all(isfinite, band.energies)
    band.stationary[2] = false
    set_climbing_image!(band, :automatic)
    relax_gneb!(
        band; climbing=true, maximum_steps=2,
        force_tolerance=1e-30, save_data_every=-1)
    @test all(isfinite, gneb_images(band))

    angles = range(0, pi; length=7)
    displaced_path = reduce(
        hcat,
        ([cos(angle), sin(angle), 0.25sin(angle)] ./
         norm([cos(angle), sin(angle), 0.25sin(angle)])
         for angle in angles))
    displaced_band = GNEB(
        transition_macrospin, displaced_path;
        name="", spring_constant=1.0meV)
    relax_gneb!(
        displaced_band; maximum_steps=500, force_tolerance=1e-27,
        save_data_every=-1, log_every=500)
    @test displaced_band.converged
    set_climbing_image!(displaced_band, :automatic)
    relax_gneb!(
        displaced_band; climbing=true, maximum_steps=500,
        force_tolerance=1e-28, save_data_every=-1, log_every=500)
    @test displaced_band.converged
    @test (maximum(displaced_band.energies) -
           displaced_band.energies[1]) / meV ≈ 1.0 atol=1e-10
    @test maximum(abs.(vec(sum(gneb_images(displaced_band).^2; dims=1)) .-
                       1)) < 1e-12

    fire_band = GNEB(
        transition_macrospin, displaced_path;
        name="", spring_constant=1.0meV)
    relax_gneb!(
        fire_band; maximum_steps=500, force_tolerance=1e-27,
        solver=:fire, save_data_every=-1, log_every=500)
    @test fire_band.converged
    @test (maximum(fire_band.energies) - fire_band.energies[1]) / meV ≈
          1.0 atol=1e-9
end

function test_automatic_transition()
    folder = mktempdir()
    result = find_transitions(
        transition_macrospin, [1.0, 0.0, 0.0];
        n_modes=1, initial_modes=reshape([0.0, 1.0, 0.0], 3, 1),
        directions=(1,), n_transitions=1, images=5,
        output_folder=folder, finite_difference_step=1e-4,
        krylov_dimension=3, maximum_restarts=3,
        residual_tolerance=1e-27, saddle_maximum_steps=120,
        saddle_step_size=0.06, gradient_tolerance=1e-27,
        eigenvalue_tolerance=1e-27, transverse_relaxation=0.0,
        branch_maximum_steps=100, branch_stopping_dmdt=1e-6,
        spring_constant=1.0meV, ordinary_maximum_steps=2,
        ordinary_force_tolerance=1e-30, climbing_maximum_steps=2,
        climbing_force_tolerance=1e-30, save_data_every=-1,
        using_time_factor=false, require_gneb_convergence=false)
    @test length(result.transitions) == 1
    @test result.attempts[1].status === :accepted
    @test result.transitions[1].valid
    @test result.transitions[1].zero_vector_count == 0
    @test !result.transitions[1].ordinary_converged
    @test !result.transitions[1].climbing_converged
    @test isfile(joinpath(folder, "minimum.ovf"))
    @test isfile(joinpath(folder, "hessian_modes.txt"))
    @test isfile(joinpath(folder, "attempts.txt"))
    @test isfile(joinpath(folder, "transitions_energy.txt"))
    @test isfile(joinpath(folder, "transition_01", "energy.txt"))
    @test isfile(joinpath(folder, "transition_01", "distance.txt"))
    @test isfile(joinpath(folder, "transition_01", "saddle.ovf"))
    @test !isdir(joinpath(folder, "candidates"))
    @test !isdir(joinpath(folder, "diagnostics"))

    resumed = find_transitions(
        transition_macrospin, [1.0, 0.0, 0.0];
        resume_result=result, directions=(1,), n_transitions=1)
    @test length(resumed.transitions) == 1
    @test length(resumed.attempts) == 1
end

@using_gpu()
test_functions("Transition Hessian", test_transition_hessian;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("Transition MMF and GNEB", test_transition_saddle_and_path;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("Automatic Transition", test_automatic_transition;
               platforms=["CPU"], precisions=[Float64])

if Base.find_package("CairoMakie") !== nothing
    @eval using CairoMakie
    @testset "Transition path plotting" begin
        saddle = SaddlePoint(
            [0.0, 1.0, 0.0], 1.0, [1.0, 0.0, 0.0],
            [-1.0, 1.0], [0.0, 0.0], 0.0, 1, true, true, 1,
            NamedTuple[])
        path = TransitionPath(
            hcat([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]),
            [0.0, 1.0, 0.0] .* meV, [0.0, 0.5, 1.0], 2,
            meV, meV, 0.0, 0.0, 0.0, true, true,
            saddle, 1, 1, true, 0, 0.0)
        output = joinpath(mktempdir(), "transition.png")
        figure = plot_transition_paths(["escape" => path]; output=output)
        @test figure !== nothing
        @test isfile(output)
        @test figure.content[1].ylabel[] == "Energy (meV)"
        figure_ev = plot_transition_paths(path; energy_unit=eV)
        @test figure_ev.content[1].ylabel[] == "Energy (eV)"
        figure_custom = plot_transition_paths(
            path; energy_unit=meV, energy_label="Relative energy (meV)")
        @test figure_custom.content[1].ylabel[] == "Relative energy (meV)"
    end
end
