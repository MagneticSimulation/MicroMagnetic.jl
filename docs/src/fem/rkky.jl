
using MicroMagnetic
using CairoMakie

mesh = FEMesh("meshes/nanorod.mesh", unit_length=1e-9);

function m0_fun(x,y,z)
    if z <= 0
        return (0,-1,0)
    else
        return (0,1,0)
    end
end

function relax_system()
    sim = Sim(mesh; driver="SD", name="rkky")

    K, A, Ms = 1e4, 1.3e-11, 8e5
    j = 2
    h = 4
    J = sqrt(2*j*A*K)
    H = h*2*K/(mu_0*Ms)
    
    set_Ms(sim, Ms)  
    add_exch(sim, A)
    add_anis(sim, K, axis=(0, 1, 0))
    add_exch_int(sim, -J; k1=1, k2=3, name="rkky") 

    add_zeeman(sim, (H, 0, 0))

    init_m0(sim, m0_fun)
    relax(sim; stopping_dmdt=0.001)  # Relax the system

    save_vtk(sim, "rkky.vtu")
    return Array(sim.spin)
end

spin = relax_system()

zs = [i+eps() for i in 1:90]
points = [[0, 0, z] for z in zs]
m = interpolate_field(mesh, spin, points)
m = reshape(m, (3, length(zs)))
phi = asin.(m[2, :])/pi*180
println("max phi = ", phi[1])   

function plot_m(zs, phi)

    fig = Figure(; size=(400, 280))
    ax = Axis(fig[1, 1]; xlabel="z (nm)", ylabel=L"\phi")

    scatterlines!(ax, zs, phi; markersize=6, color=:blue, markercolor=:orange)
    return fig
end

fig = plot_m(zs, phi)
save("rkky.png", fig)





