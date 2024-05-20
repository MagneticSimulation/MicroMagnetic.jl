using NuMag
using Printf
using Test
using MPI


MPI.Init()

NuMag.cuda_using_double(true)

NuMag.using_multiple_gpus() #comment out this line if you have multiple gpus

#julia test_neb_gpu_mpi.jl //for single MPI process
#mpiexec -np 2 julia test_neb_gpu_mpi.jl //for two MPI processes

function creat_sim()
    mesh =  FDMeshGPU(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh, name="neb", driver="none",save_data=false)
    set_Ms(sim, 8e5)

    init_m0(sim, (0.6, 0, -0.8))

    add_exch(sim, 1.3e-11)
    add_anis(sim, 5e4, axis=(1,0,0), name="Kx")
    add_anis(sim, -2e4, axis=(0,0,1), name="Kz")
    return sim
end

sim = creat_sim()

neb = NEB_MPI(sim, [ (1,0,0), (0,1,1), (-1,0,0)], [5, 5]; name="test")

relax(neb, maxsteps=200, stopping_dmdt=0.01, save_ovf_every=-1)
max_energy = maximum(neb.energy)
expected_energy_diff = 5e4*5e-9^3

all_max_energy = [0.0]

if neb.comm_rank == 0
    all_max_energy[1] = MPI.Reduce(max_energy, MPI.MAX, 0, MPI.COMM_WORLD)
else
    MPI.Reduce(max_energy, MPI.MAX, 0, MPI.COMM_WORLD)
end
MPI.Bcast!(all_max_energy, 0, MPI.COMM_WORLD)

if neb.comm_rank == 0
    expected_energy_diff = 5e4*5e-9^3
    final_energy_diff = all_max_energy[1] - neb.energy_cpu[1]
    println("expected_diff=$expected_energy_diff, got_energy_diff=$final_energy_diff")
    @test abs(final_energy_diff - expected_energy_diff)/expected_energy_diff < 1e-8
end
