using MicroMagnetic
using CUDA
using CairoMakie

@assert CUDA.functional()

# Exercise the full user path so that every artifact needed at runtime
# (MKL, FFmpeg, fonts, ...) is present in the image: GPU simulation with
# demag, plotting, and movie export.
MicroMagnetic.set_backend("CUDA")
set_precision(Float32)

mesh = FDMesh(; dx=4e-9, dy=4e-9, dz=4e-9, nx=12, ny=12, nz=1)
sim = Sim(mesh; name="_smoke")
set_Ms(sim, 8e5)
init_m0(sim, (1, 0.1, 0.1))
add_exch(sim, 1.3e-11)
add_demag(sim)
add_zeeman(sim, (0, 0, 1e5))
relax(sim; max_steps=2, stopping_dmdt=0.0, save_m_every=1)

fig = MicroMagnetic.plot_m(sim; arrows=(3, 3))
save("/tmp/_smoke.png", fig)
MicroMagnetic.ovf2movie("_smoke_LLG"; output="/tmp/_smoke.mp4", framerate=2)

println("PROBE_OK png=$(filesize("/tmp/_smoke.png")) mp4=$(filesize("/tmp/_smoke.mp4"))")
