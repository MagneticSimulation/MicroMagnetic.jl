using JuMag
using Test
using Printf
using PyCall
using NPZ

function savem()
    sim = Sim(nx=32,ny=64,nz=16,Ms=1e5,GPU=false)
    init_m0(sim,(cos(5/3*pi),sin(5/3*pi),0))
    save_ovf(sim,"test_tools")
end


function is_very_close(a, b)
    if abs(a-b) < 1e-12
        return true
    else
        return false
    end
end

function test_sum()
    ovf = read_ovf("test_tools")
    mx, my, mz = sum_ovf(ovf,axis=ez)
    @test is_very_close(mx[1],16*cos(5/3*pi))
    @test is_very_close(my[1],16*sin(5/3*pi))
    @test is_very_close(mx[1],0)

    mx, my, mz = sum_ovf(ovf,axis=ex)
    @test is_very_close(mx[1],-32*sin(5/3*pi))
    @test is_very_close(my[1],0)
    @test is_very_close(mz[1],-32*cos(5/3*pi))

    mx, my, mz = sum_ovf(ovf,axis=ey)
    @test is_very_close(mx[1],64*cos(5/3*pi))
    @test is_very_close(my[1],0)
    @test is_very_close(mz[1],-64*sin(5/3*pi))
end

function test_projection()
    ovf = read_ovf("test_tools")
    m = reshape(ovf.data, (3,32,64,16))
    m_padded = JuMag.vector_padding(m, 128,128,128)
    m = JuMag.vector_field_projection(m_padded, 0.0, "alpha")
    mx, my, mz = m[1,:,:], m[2,:,:], m[3,:,:]
    @test abs(mx[48,64] - 0) < 1e-6
    @test abs(mx[49,64] - 16*cos(5/3*pi)) < 1e-6
    @test abs(my[48,64] - 0) < 1e-6
    @test abs(my[49,64] - 16*sin(5/3*pi)) < 1e-6
    @test abs(mz[48,64] - 0)< 1e-6
    @test abs(mz[49,64] - 0) < 1e-6
end




savem()
test_projection()

#OVF2LTEM("test_tools",N=128)
OVF2XRAY("test_tools")
OVF2MFM("test_tools")
rm("test_tools.ovf")


#TODO:phase test
#=function F0(x::Float64, y::Float64)
    eps = 1e-10
    if abs(x)<eps && abs(y)<eps
        return 0
    end
    return x*log(x^2+y^2) - 2*x + 2*y*atan(x/y)
end

function phim_uniformly_magnetized_slab(x, y, mx, my, Lx, Ly, Lz, Ms)
    mu0 = 4*pi*1e-7
    Phi0 = 2.067833e-15
    a = F0(x-Lx/2, y-Ly/2) - F0(x+Lx/2, y-Ly/2) - F0(x-Lx/2, y+Ly/2) + F0(x+Lx/2, y+Ly/2)
    b = F0(y-Ly/2, x-Lx/2) - F0(y+Ly/2, x-Lx/2) - F0(y-Ly/2, x+Lx/2) + F0(y+Ly/2, x+Lx/2)
    #println(a, " ", b, "  ", my*b-mx*a)
    return mu0*Ms*Lz/(4*Phi0)*(my*b-mx*a)*1e-18
end

function phase_theory()
    mx = cos(5/3*pi)
    my = sin(5/3*pi)
    Lx = 32
    Ly = 64
    Lz = 16
    Nx, Ny = 160, 160
    Ms = 1e5
    phi = zeros(Nx, Nx)
    for i = 1:Nx, j=1:Ny
        x = (i/Nx-0.5)*160
        y = (j/Ny-0.5)*160
        phi[i, j] = phim_uniformly_magnetized_slab(x, y, mx, my, Lx, Ly, Lz, Ms)
    end
    
    return phi
end

theory = phase_theory()
phase,intensity = OVF2LTEM("test_tools",Nx=160,Ny=160)

#@test isapprox(theory[])
@info("test tools passed!")=#