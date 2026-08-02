using LinearAlgebra
using MicroMagnetic
using Test

if !isdefined(Main, :test_functions)
    include("../test_utils.jl")
end

# Single-spin easy-axis system with a known 1 meV barrier between (1,0,0) and
# (-1,0,0) through the saddle at (0,+/-1,0). Self-contained so this file can
# be included independently of test_transition.jl.
function _solver_macrospin()
    mesh = CubicMesh(; nx=1, ny=1, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="solver_macrospin", save_data=false)
    set_mu_s(sim, mu_B)
    add_anis(sim, 1.0meV; axis=(1, 0, 0), name="easy_x")
    add_anis(sim, -0.5meV; axis=(0, 0, 1), name="hard_z")
    init_m0(sim, (1, 0, 0))
    return sim
end

# 3x3 affine atomistic system (exchange + DMI + zeeman, one pinned site) used
# to exercise the 2N basis on a multi-site Hamiltonian with an exact
# Hessian-vector product.
function _solver_affine_atomistic()
    mesh = CubicMesh(; nx=3, ny=3, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="solver_affine", save_data=false)
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

# 7-image path from (1,0,0) to (-1,0,0) with a small out-of-plane kick so the
# band relaxes onto the great-circle minimum-energy path.
function _solver_displaced_path()
    angles = range(0, pi; length=7)
    return reduce(
        hcat,
        ([cos(a), sin(a), 0.25sin(a)] /
         norm([cos(a), sin(a), 0.25sin(a)]) for a in angles))
end

# Spirit Atlas must reach the same analytic 1 meV barrier as VP/FIRE.
function test_gneb_solvers()
    path = _solver_displaced_path()
    for solver in (:lbfgs_atlas,)
        max_steps = 500
        band = GNEB(_solver_macrospin, path;
                    name="", spring_constant=1.0meV)
        relax_gneb!(band; maximum_steps=max_steps, force_tolerance=1e-27,
                    solver=solver, save_data_every=-1, log_every=max_steps)
        @test band.converged
        set_climbing_image!(band, :automatic)
        relax_gneb!(band; climbing=true, maximum_steps=500,
                    force_tolerance=1e-28, solver=solver,
                    save_data_every=-1, log_every=500)
        @test band.converged
        @test first(band.history).time_step < 0.05
        @test first(band.history).time_step >= 1e-10
        @test (maximum(band.energies) - band.energies[1]) / meV ≈ 1.0 atol=1e-6
        @test maximum(abs.(vec(sum(gneb_images(band).^2; dims=1)) .- 1)) < 1e-12
    end
end

# --- MMF saddle search: the climbing phase must converge to the index-one ---
# saddle at (0,+/-1,0).
function test_saddle_solvers()
    for solver in (:lbfgs_atlas, :fire)
        saddle = find_saddle(
            _solver_macrospin, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0];
            direction=1, maximum_steps=300, step_size=0.06,
            mode_update_every=3, gradient_tolerance=1e-27,
            eigenvalue_tolerance=1e-27, transverse_relaxation=0.0,
            n_modes=2, finite_difference_step=1e-4, krylov_dimension=3,
            maximum_restarts=3, residual_tolerance=1e-27,
            solver=solver)
        @test saddle.converged
        @test saddle.first_order
        @test saddle.negative_mode_count == 1
        @test abs(saddle.spin[1]) < 1e-4
    end
end

function test_spirit_gneb_formulas()
    path = _solver_displaced_path()
    band = GNEB(_solver_macrospin, path; name="", spring_constant=1.0meV)
    expected_rx = cumsum([0.0; [acos(clamp(dot(path[:, i - 1], path[:, i]), -1, 1))
                                for i in 2:size(path, 2)]])
    @test band.reaction_coordinate ≈ expected_rx atol=1e-14

    image = 3
    forward = path[:, image + 1] - path[:, image]
    backward = path[:, image] - path[:, image - 1]
    e, ep, em = band.energies[image], band.energies[image + 1], band.energies[image - 1]
    dmax = max(abs(ep - e), abs(em - e))
    dmin = min(abs(ep - e), abs(em - e))
    expected = ep > em ? dmax * forward + dmin * backward :
                         dmin * forward + dmax * backward
    spin = path[:, image]
    expected .-= dot(expected, spin) .* spin
    expected ./= norm(expected)
    @test Array(view(band.tangents, :, image)) ≈ expected atol=1e-13
    lengths = MicroMagnetic._gneb_hermite_lengths(band, 0.5)
    @test all(isfinite, lengths)
    @test any(abs.(lengths[2:end] .- diff(band.reaction_coordinate)) .> 1e-12)

    set_image_type!(band, 3, :climbing)
    set_image_type!(band, 4, :climbing)
    @test band.climbing_image == 4
    @test count(==(:climbing), band.image_types) == 1
    set_image_type!(band, 4, :falling)
    @test band.climbing_image == 0

    MicroMagnetic._update_gneb!(band; spring_force_ratio=0.5,
                                path_shortening_constant=1e-4meV)
    @test all(isfinite, band.forces)
    @test band.maximum_force >= band.maximum_component
    @test_throws ErrorException relax_gneb!(band; solver=:ncg_oso,
                                             maximum_steps=1)

    folder = mktempdir()
    band.name = joinpath(folder, "units")
    push!(band.history,
          (step=0, maximum_force=band.maximum_force,
           maximum_component=band.maximum_component,
           climbing_image=band.climbing_image, energies=copy(band.energies),
           reaction_coordinate=copy(band.reaction_coordinate)))
    MicroMagnetic._write_gneb_history(band)
    energy_lines = readlines(band.name * "_energy.txt")
    distance_lines = readlines(band.name * "_distance.txt")
    diagnostic_lines = readlines(band.name * "_gneb_history.txt")
    @test startswith(energy_lines[1], "# steps E_total_1")
    @test count(==("<J>"), split(energy_lines[2])) == size(path, 2)
    @test startswith(distance_lines[1], "# steps distance_1")
    @test all(==("<>"), split(distance_lines[2])[2:end])
    @test occursin("maximum_site_force_J", diagnostic_lines[1])
end

# --- 2N tangent basis: eigenvalues and eigenvectors must match the ---
# matrix-free :projective path. Exercises both the single-site and the
# multi-site affine Hamiltonian (exact Hessian-vector action).
function test_tangent_basis_2N()
    # Macrospin: dim_2N = 2.
    modes_p = compute_hessian_modes(
        _solver_macrospin, [0.0, 1.0, 0.0];
        n_modes=2, finite_difference_step=1e-4, krylov_dimension=3,
        maximum_restarts=3, residual_tolerance=1e-27,
        tangent_basis=:projective)
    modes_2 = compute_hessian_modes(
        _solver_macrospin, [0.0, 1.0, 0.0];
        n_modes=2, finite_difference_step=1e-4, krylov_dimension=3,
        maximum_restarts=3, residual_tolerance=1e-27,
        tangent_basis=:full_2N)
    @test modes_2.eigenvalues ≈ modes_p.eigenvalues rtol=1e-6 atol=1e-30
    @test all(isfinite, modes_2.eigenvalues)
    for j in 1:2
        @test abs(dot(modes_2.eigenvectors[:, j],
                      modes_p.eigenvectors[:, j])) ≈ 1.0 atol=1e-6
    end

    # Affine atomistic 3x3 (8 active sites -> dim_2N = 16).
    spin = repeat([0.0, 0.0, 1.0], 9)
    modes_p2 = compute_hessian_modes(
        _solver_affine_atomistic, spin;
        n_modes=4, finite_difference_step=1e-4, krylov_dimension=28,
        maximum_restarts=5, residual_tolerance=1e-27,
        tangent_basis=:projective)
    modes_2_2 = compute_hessian_modes(
        _solver_affine_atomistic, spin;
        n_modes=4, finite_difference_step=1e-4, krylov_dimension=28,
        maximum_restarts=5, residual_tolerance=1e-27,
        tangent_basis=:full_2N)
    @test modes_2_2.eigenvalues ≈ modes_p2.eigenvalues rtol=1e-5 atol=1e-30
    @test all(isfinite, modes_2_2.eigenvalues)
    for j in 1:4
        # For near-degenerate modes the 2N and projective paths may pick
        # different eigenvectors within the degenerate subspace. Check that
        # the 2N eigenvector overlaps with at least one projective eigenvector
        # of the same Hessian, which is the correct equivalence test.
        overlaps = [abs(dot(modes_2_2.eigenvectors[:, j],
                            modes_p2.eigenvectors[:, k])) for k in 1:4]
        @test maximum(overlaps) > 0.99
    end
end

# --- Spirit-compatible force enhancements must not break convergence. ---
function test_gneb_force_enhancements()
    path = _solver_displaced_path()
    band = GNEB(_solver_macrospin, path;
                name="", spring_constant=1.0meV)
    relax_gneb!(band; maximum_steps=800, force_tolerance=1e-27,
                spring_force_ratio=0.5, path_shortening_constant=1e-4meV,
                save_data_every=-1, log_every=800)
    @test band.converged
    @test (maximum(band.energies) - band.energies[1]) / meV ≈ 1.0 atol=1e-5
    @test maximum(abs.(vec(sum(gneb_images(band).^2; dims=1)) .- 1)) < 1e-12
end

@using_gpu()
test_functions("GNEB solvers", test_gneb_solvers;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("MMF saddle solvers", test_saddle_solvers;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("Spirit GNEB formulas", test_spirit_gneb_formulas;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("Transition tangent basis", test_tangent_basis_2N;
               platforms=["CPU", "CUDA"], precisions=[Float64])
test_functions("GNEB force enhancements", test_gneb_force_enhancements;
               platforms=["CPU", "CUDA"], precisions=[Float64])
