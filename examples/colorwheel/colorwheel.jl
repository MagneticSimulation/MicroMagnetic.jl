using JuMag
using PyCall

function save_colorwheel()
    mesh = FDMeshGPU(nx=100,ny=100,nz=1)
    sim = Sim(mesh)
    set_Ms_cylindrical(sim,1e5)

    function m0_fun(i,j,k,dx,dy,dz)
        R = 50*dx
        xc,yc = 50*dx,50*dy
        x,y = i*dx,j*dy
        rx,ry = x-xc,y-yc
        r = sqrt(rx^2+ry^2)
        if r>=R
            return (0,0,0)
        end
        mx,my = rx/R,ry/R
        mz = sqrt(1-mx^2-my^2)
        return (mx,my,mz)
    end

    init_m0(sim,m0_fun)
    save_ovf(sim,"colorwheel")
    JuMag.show_mag("colorwheel.ovf", rgb=true)

    rm("colorwheel.ovf")
end

save_colorwheel()
