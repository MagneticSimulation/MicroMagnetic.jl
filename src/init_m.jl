
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

"""
    init_m0_skyrmion(sim::AbstractSim, center::Tuple, R::Float64; ratio=0.7, p=-1, c=1, type="B")

Set the magnetization with skyrmions. Note that this function can be called mulitple times to add more skyrmons.

center :  the skyrmion center, should be a Tuple. For example, center = (50e-9,50e-9)

R : the skyrmion radius.

ratio : ratio=w/R where w is the width of domain wall. By default ratio = 0.7

p : polarity, +1 --> core up; -1 --> core down

c : chirality, +1 --> clockwise; -1 --> counter-clockwise

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
    b = reshape(sim.spin,(3,nx,ny,nz))

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
                b[1,i,j,:] .= -sin(theta)*x/r
                b[2,i,j,:] .= -sin(theta)*y/r
                b[3,i,j,:] .= cos(theta)
            elseif type == "B"
                b[1,i,j,:] .= sin(theta)*y/r
                b[2,i,j,:] .= -sin(theta)*x/r
                b[3,i,j,:] .= cos(theta)
            end
      end
  end
  #normalise(sim.spin, sim.nxyz)
  copyto!(sim.prespin, sim.spin)
  return nothing
end
