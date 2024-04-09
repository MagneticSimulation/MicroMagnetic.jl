
export init_m0_random, init_m0_skyrmion

"""
    init_m0_random(sim::MicroSim)

Set the initial magnetization with random direction.
"""
function init_m0_random(sim::AbstractSim)
    function random_m(i,j,k,dx,dy,dz)
        return 2*rand(3).-1
    end
    init_m0(sim, random_m)
end

function init_m0_skyrmion(sim::AbstractSim)
    mesh = sim.mesh
    nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
    dx,dy,dz = mesh.dx,mesh.dy,mesh.dz
    xc, yc = dx*(nx+1)/2, dy*(ny+1)/2
    R = 1/2 * min(nx*dx, ny*dy)
    init_m0_skyrmion(sim, (xc,yc), R)
end

"""
    init_m0_skyrmion(sim::AbstractSim, center::Tuple, R::Float64; ratio=0.7, p=-1, c=1, type="B")

Set the magnetization with skyrmions. Note that this function can be called mulitple times to add more skyrmons.

center :  the skyrmion center, should be a Tuple. For example, center = (50e-9,50e-9)

R : the skyrmion radius.

ratio : ratio=w/R where w is the width of domain wall. By default ratio = 0.7

p : polarity, +1 --> core up; -1 --> core down

c : chirality, +1 --> lefthand,for positive D; -1 --> righthand,for negative D

type : "B" or "N", representing Bloch or Neel skyrmions.

For example:
```julia
    init_m0_skyrmion(sim, (50e-9,50e-9), 2e-8, ratio=0.5, p=-1, c=1, type="B")
```
"""
function init_m0_skyrmion(sim::AbstractSim, center::Tuple, R::Float64; ratio=0.7, p=-1, c=1, type="B")
    if type!="N" && type!="B"
        @info("add_skyrmion: type should be \"N\" or \"B\" ")
        return nothing
    end
    x0,y0 = center[1],center[2]
    w = ratio*R
    mesh = sim.mesh
    nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
    dx,dy,dz = mesh.dx,mesh.dy,mesh.dz
    m = Array(sim.spin)
    b = reshape(m, (3,nx,ny,nz))
    Ms= Array(sim.Ms)

    eps = 1e-6
    for i=1:nx, j=1:ny
        x = (i+eps)*dx-x0
        y = (j+eps)*dy-y0
        r = sqrt(x^2+y^2)
        if x^2+y^2 <= 2*R^2
            theta = 2*atan(sinh(R/w)/sinh(r/w))
            if p == 1
                theta = theta + pi
            end
            if c == 1
                theta = -theta
            end
            if type == "N"
                mx = -sin(theta)*x/r
                my = -sin(theta)*y/r
                mz = cos(theta)
            elseif type == "B"
                mx = sin(theta)*y/r
                my = -sin(theta)*x/r
                mz = cos(theta)
            end

            for k=1:nz
                id = index(i,j,k,nx,ny,nz)
                if Ms[id] == 0
                    continue
                end

                b[1,i,j,k] = mx
                b[2,i,j,k] = my
                b[3,i,j,k] = mz
            end
        end
    end
    #normalise(sim.spin, sim.n_total)
    copyto!(sim.prespin, m)
    copyto!(sim.spin, m)
    return nothing
end

function rotation_2d(theta::Float64, v::Array{Float64,1})
    ct, st = cos(deg2rad(theta)), sin(deg2rad(theta))
    vnew = zeros(2)
    vnew[1] = ct*v[1] - st*v[2]
    vnew[2] = st*v[1] + ct*v[2]
    return vnew
end

"""
  init_m0_skyrmion_lattice(sim::AbstractSim, k::Tuple, lambda::Float64; p=-1, c=1, type="B")

Use 3-Q approximation to approach skyrmion lattice.

lambda: period

p: polarity. 1(up) or -1(down)
c: chirality. "r" (left-handed) or "l"(left-handed)
type: "B"(Bloch) or "N"(Neel)
"""

function init_m0_skyrmion_lattice(sim::JuMag.AbstractSim, lambda::Float64; p=-1, c="l", type="B")
    mesh = sim.mesh
    nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
    dx,dy,dz = mesh.dx,mesh.dy,mesh.dz
    spin = zeros(3*nx*ny*nz)
    b = reshape(spin,(3,nx,ny,nz))

    #normalize
    k=(1,0)
    kx, ky = k[1]/hypot(k[1], k[2]), k[2]/hypot(k[1], k[2])
    #wave vector
    Q1 = [2*pi*kx/lambda, 2*pi*ky/lambda]
    Q2 = rotation_2d(120., Q1)
    Q3 = rotation_2d(120., Q2)

    xc, yc = nx/2*dx, ny/2*dy
    for i=1:nx, j=1:ny
        x, y = i*dx-xc, j*dy-yc
        phase1 = Q1[1]*x + Q1[2]*y
        phase2 = Q2[1]*x + Q2[2]*y 
        phase3 = Q3[1]*x + Q3[2]*y 
        mx = sin(phase1)*cos(0) + sin(phase2)*cos(2*pi/3) + sin(phase3)*cos(4*pi/3)
        my = sin(phase1)*sin(0) + sin(phase2)*sin(2*pi/3) + sin(phase3)*sin(4*pi/3)
        if c=="l"
            mx,my = JuMag.rotation_2d(180., [mx, my])
        end

        if type == "B"
          mx , my = JuMag.rotation_2d(90., [mx, my])
          b[1,i,j,:] .= mx
          b[2,i,j,:] .= my
          b[3,i,j,:] .= cos(phase1) + cos(phase2) + cos(phase3)
        elseif type == "N"
          b[1,i,j,:] .= mx
          b[2,i,j,:] .= my
          b[3,i,j,:] .= cos(phase1) + cos(phase2) + cos(phase3)
        end


    end

    if p==-1
      spin .*= -1
    end
    init_m0(sim, spin)

end

"""
  init_hard_bubble(sim::AbstractSim; R::Float64, delta::Float64, center::Any="default", S::Int=10, Q::Int=-1, phi0::Float64=pi/2)

Initialize a bubble with arbitrary topological charge                                            
Ref: J.Appl.Phys.120,233901(2016)                                                                
Parameters: R:radius in nanometer; 
            delta: domain wall width= sqrt(A/k);
            center=(100e-9,100e-9);  
            S: winding number; 
            Q:polarity or "topological quantities"; 
            phi0= initial phase                                                                    

"""

function init_hard_bubble(sim::AbstractSim; R::Float64, delta::Float64, center::Any="default", S::Int=10, Q::Int=-1, phi0::Float64=pi/2)
    mesh = sim.mesh
    nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
    dx,dy,dz = mesh.dx,mesh.dy,mesh.dz
    Ms = sim.Ms

    if center=="default"
        center = ((nx+1)*dx/2, (ny+1)*dy/2)
    end
    x0,y0 = center[1],center[2]
    b = reshape(sim.spin,(3,nx,ny,nz))

    eps = 1e-6
    for i=1:nx, j=1:ny
        x = (i+eps)*dx-x0
        y = (j+eps)*dy-y0
        r = sqrt(x^2+y^2)
        if r>1.2*R
            continue
        end
        psi = atan(y,x)
        theta = 2*atan(exp(Q*(r-R)/delta))
        phi = S*psi + phi0
        mx = sin(theta)*cos(phi)
        my = sin(theta)*sin(phi)
        mz = cos(theta)

        for k=1:nz
            id = index(i,j,k,nx,ny,nz)
            if Ms[id] == 0
                continue
            end

            b[1,i,j,k] = mx
            b[2,i,j,k] = my
            b[3,i,j,k] = mz
        end
    end
end
