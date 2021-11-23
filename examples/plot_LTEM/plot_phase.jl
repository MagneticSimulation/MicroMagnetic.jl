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
    get_magnetic_phase(m, alphas=[0.0, 45], betas=[0.0, 45], N=256)
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

show_phase("X_0.npy")