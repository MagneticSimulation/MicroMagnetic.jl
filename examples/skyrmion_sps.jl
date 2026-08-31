# @title Systematic Skyrmion Saddle-Point Search
# @description Discover skyrmion transitions using Hessian modes, MMF, and GNEB
# @tags atomistic skyrmion sps hessian mmf gneb gpu

# References:
# G. P. Müller et al., Phys. Rev. Lett. 121, 197202 (2018).
# https://doi.org/10.1103/PhysRevLett.121.197202
# P. F. Bessarab et al., Comput. Phys. Commun. 196, 335-347 (2015).
# https://doi.org/10.1016/j.cpc.2015.07.001

using MicroMagnetic
using CairoMakie
using Printf

@using_gpu()

const N = 40
const J = 1.0meV
const D = 0.45J
const MU_S = mu_B
const B_Z = 0.8D^2 / (MU_S * J)
const OUTPUT = joinpath(@__DIR__, "skyrmion_sps")

function build_sim()
    mesh = CubicMesh(; nx=N, ny=N, nz=1, pbc="open")
    sim = Sim(mesh; driver="SD", name="skyrmion_sps", save_data=false)
    set_mu_s(sim, MU_S)
    add_exch(sim, J)
    add_dmi(sim, D)
    add_zeeman(sim, (0, 0, B_Z))
    return sim
end

function initial_skyrmion(i, j, k, dx, dy, dz)
    x = i - (N + 1) / 2
    y = j - (N + 1) / 2
    r = hypot(x, y)
    theta = pi * clamp(1 - r / 8, 0, 1)
    r == 0 && return (0.0, 0.0, -1.0)
    return (sin(theta) * y / r, -sin(theta) * x / r, cos(theta))
end

mkpath(OUTPUT)
sim = build_sim()
init_m0(sim, initial_skyrmion)
relax(sim; max_steps=50_000, stopping_dmdt=1e-5,
      save_data_every=-1, save_m_every=-1)
minimum = Float64.(Array(sim.spin))

result = find_transitions(
    build_sim, minimum; n_modes=8, directions=(-1, 1),
    exploration_depth=1, n_transitions=3, images=21,
    krylov_dimension=64, output_folder=OUTPUT)

println(" mode dir      Q(left)     Q(right)  barrier(meV)  converged")
for path in result.transitions
    q_left = compute_skyrmion_number(path.images[:, 1], sim.mesh)
    q_right = compute_skyrmion_number(path.images[:, end], sim.mesh)
    @printf(" %4d %+3d  %+11.5f  %+11.5f  %12.6f  %s\n",
            path.mode_index, path.direction, q_left, q_right,
            path.forward_barrier / meV,
            path.ordinary_converged && path.climbing_converged)
end

isempty(result.transitions) || plot_transition_paths(
    result; output=joinpath(OUTPUT, "transition_paths.png"))
