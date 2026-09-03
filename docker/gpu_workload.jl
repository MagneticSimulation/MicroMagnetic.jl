using MicroMagnetic
using CUDA

MicroMagnetic.set_backend("CUDA")
set_precision(Float32)

# Bake the hot GPU path: FDMesh field kernels (demag FFT, exchange,
# anisotropy, zeeman, STT) plus the LLG driver, so that CUDA kernels are
# compiled for the GPU architecture present at image build time.
#
# NOTE: keep this file identical to what was validated locally; extend it
# only after the added calls are verified to work on CUDA (e.g. the SD
# driver currently errors on CUDA, and the thermal/SOT GPU paths are not
# covered here yet).
mesh = FDMesh(; dx=4e-9, dy=4e-9, dz=4e-9, nx=16, ny=16, nz=2)
sim = Sim(mesh; name="_sysimage_wl", driver="LLG")
set_Ms(sim, 8e5)
init_m0(sim, (1, 0.25, 0.1))
add_exch(sim, 1.3e-11)
add_demag(sim)
add_anis(sim, 5e2, axis=(0, 1, 0))
add_zeeman(sim, (0, 0, 1e4))
add_stt(sim, model=:zhang_li, P=0.5, Ms=8e5, xi=0.05, J=(1e12, 0, 0))
relax(sim; max_steps=3, stopping_dmdt=0.0, save_m_every=-1)
run_sim(sim; steps=2, dt=1e-12, save_m_every=-1)

println("GPU_WORKLOAD_DONE")
