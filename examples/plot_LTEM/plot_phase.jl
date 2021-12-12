using JuMag
using PyCall
using NPZ

# create a uniform slab 

function savem()
    mesh = FDMesh(nx=32,ny=64,nz=16)
    sim = Sim(mesh)
    set_Ms(sim, 1e5)
    init_m0(sim,(cos(5/3*pi),sin(5/3*pi),0))
    save_ovf(sim,"test_tools")
end


# generate phase map
function gen_magnetic_phase()
    ovf = read_ovf("test_tools")
    m = reshape(ovf.data, (3,32,64,16))
    A = compute_vector_potential(m, 62*2)
    phase = compute_phase(A, gamma=0.0)
    npzwrite("phase_gamma_0.npy", phase)

    phase = compute_phase(A, gamma=pi/3)
    npzwrite("phase_gamma_60.npy", phase)

    # Or
    phase = compute_magnetic_phase(m, gamma=pi/3)
    npzwrite("phase2_gamma_60.npy", phase)
end

# plot phase map
function show_phase(fname)
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    plt = pyimport("matplotlib.pyplot")

    phase = npzread(fname)
    plt.imshow(transpose(phase), origin="lower")
    plt.colorbar()
    plt.savefig(fname[1:end-4]*".png", dpi=300)    
    plt.close()
end

savem()

gen_magnetic_phase()

show_phase("phase_gamma_0.npy")
show_phase("phase_gamma_60.npy")
show_phase("phase2_gamma_60.npy")