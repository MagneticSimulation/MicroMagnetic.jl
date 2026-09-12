# Diamond domain formation on an adaptive composite grid, after
# García-Cervera & Roma, "Adaptive Mesh Refinement for Micromagnetics
# Simulations", IEEE Trans. Magn. 42(6), 1648 (2006), Section IV:
# a rectangular sample with strong out-of-plane anisotropy develops a
# diamond domain structure; the composite grid (base 128x32 plus three
# refinement levels, equivalent to a 1024x256 uniform grid) follows the
# walls, vortices and the sample boundary.
#
# NOTE on the physics: the paper's sample has infinite thickness and its
# stray field is the 2D log-kernel solution. This demo uses the package's
# thin-film (finite-thickness) demag kernel with a small dz instead, so the
# domain pattern is similar in spirit but not identical to Fig. 2 of the
# paper.

using MicroMagnetic

set_backend("cpu")
set_precision(Float64)

mesh = FDMesh(dx=7.8125e-9, dy=31.25e-9, dz=5e-9, nx=128, ny=32, nz=1)
# dt must stay at or below ~0.5 ps: the demag field is updated explicitly
# once per step, and at dt = 1 ps the self-consistency error pumps energy
# (~dt^2 per step) instead of relaxing.
amr = AMRSim(mesh; levels=4, Ms=8e5, A=1.3e-11, Ku=5e4, alpha=0.02,
             demag=true, dt=5e-13, remesh_interval=50,
             name="amr_diamond", save_data=true)

init_m0(amr, (1, 0.2, 0))
relax(amr; dt=5e-13, maxsteps=3000, stopping_dmdt=1.0, save_m_every=500,
      save_m_path="amr_diamond", verbose=true)
save_ovf(amr, "amr_diamond_final")
@info "done: steps=$(amr.nsteps) max_dmdt=$(amr.maxdmdt) patches=$(sum(length.(amr.patches)))"
