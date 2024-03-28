using FFTW
using LinearAlgebra

mutable struct Demag
  nx_fft::Int64
  ny_fft::Int64
  nz_fft::Int64
  tensor_xx::Array{Float64, 3}
  tensor_yy::Array{Float64, 3}
  tensor_zz::Array{Float64, 3}
  tensor_xy::Array{Float64, 3}
  tensor_xz::Array{Float64, 3}
  tensor_yz::Array{Float64, 3}
  mx::Array{Float64, 3}
  my::Array{Float64, 3}
  mz::Array{Float64, 3}
  Mx::Array{Complex{Float64}, 3}
  My::Array{Complex{Float64}, 3}
  Mz::Array{Complex{Float64}, 3}
  Hx::Array{Complex{Float64}, 3}
  Hy::Array{Complex{Float64}, 3}
  Hz::Array{Complex{Float64}, 3}
  m_plan::Any
  h_plan::Any
  field::Array{Float64, 1}
  energy::Array{Float64, 1}
  name::String
end

function init_demag(sim::MicroSim, Nx::Int, Ny::Int, Nz::Int)
  mesh = sim.mesh
  max_size = max(mesh.dx, mesh.dy, mesh.dz)
  dx = mesh.dx/max_size
  dy = mesh.dy/max_size
  dz = mesh.dz/max_size

  nx = mesh.nx
  ny = mesh.ny
  nz = mesh.nz

  cn = 10
  nx_fft = mesh.nx > cn ? 2*mesh.nx : 2*mesh.nx - 1
  ny_fft = mesh.ny > cn ? 2*mesh.ny : 2*mesh.ny - 1
  nz_fft = mesh.nz > cn ? 2*mesh.nz : 2*mesh.nz - 1

  tensor = zeros(nx, ny, nz)
  long_tensor = zeros(nx_fft, ny_fft, nz_fft)

  compute_demag_tensors(tensor, demag_tensor_xx, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor)
  tensor_xx = real(FFTW.rfft(long_tensor)) #tensors xx, xy, .. should be pure real

  compute_demag_tensors(tensor, demag_tensor_yy, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor)
  tensor_yy = real(FFTW.rfft(long_tensor))

  compute_demag_tensors(tensor, demag_tensor_zz, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor)
  tensor_zz = real(FFTW.rfft(long_tensor))

  compute_demag_tensors(tensor, demag_tensor_xy, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor, tx=true, ty=true)
  tensor_xy = real(FFTW.rfft(long_tensor))

  compute_demag_tensors(tensor, demag_tensor_xz, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor, tx=true, tz=true)
  tensor_xz = real(FFTW.rfft(long_tensor))

  compute_demag_tensors(tensor, demag_tensor_yz, nx, ny, nz, Nx, Ny, Nz, dx, dy, dz)
  fill_demag_tensors(long_tensor, tensor, ty=true, tz=true)
  tensor_yz = real(FFTW.rfft(long_tensor))

  mx = zeros(nx_fft, ny_fft, nz_fft)
  my = zeros(nx_fft, ny_fft, nz_fft)
  mz = zeros(nx_fft, ny_fft, nz_fft)
  lenx = (nx_fft%2>0) ? nx : nx+1
  Mx = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  My = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Mz = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  m_plan = FFTW.plan_rfft(mx)

  Hx = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Hy = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Hz = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)

  h_plan = FFTW.plan_irfft(Hx, nx_fft)

  field = zeros(Float64, 3*sim.n_nodes)
  energy = zeros(Float64, sim.n_nodes)
  demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, mx, my, mz, Mx, My, Mz,
                Hx, Hy, Hz, m_plan, h_plan, field, energy, "Demag")
  return demag
end

function copy_spin_to_m(mx::Array{Float64, 3}, my::Array{Float64, 3}, mz::Array{Float64, 3},
	                   spin::Array{Float64, 1}, Ms::Array{Float64, 1}, nx::Int64, ny::Int64, nz::Int64)
  fill!(mx, 0.0)
  fill!(my, 0.0)
  fill!(mz, 0.0)
  for k=1:nz, j=1:ny, i=1:nx
	  p = (k-1) * nx*ny + (j-1) * nx + i
	  mx[i,j,k] = spin[3*p-2]*Ms[p]
	  my[i,j,k] = spin[3*p-1]*Ms[p]
	  mz[i,j,k] = spin[3*p]*Ms[p]
  end
end

function extract_effective_field(field::Array{Float64, 1}, fx::Array{Float64, 3}, fy::Array{Float64, 3},
	                          fz::Array{Float64, 3}, nx::Int64, ny::Int64, nz::Int64)
  for k=1:nz, j=1:ny, i=1:nx
    p = (k-1) * nx*ny + (j-1) * nx + i
	field[3*p-2] = -1.0*fx[i,j,k]
	field[3*p-1] = -1.0*fy[i,j,k]
	field[3*p] = -1.0*fz[i,j,k]
  end
end

function effective_field(demag::Demag, sim::MicroSim, spin::Array{Float64, 1}, t::Float64)
  nx, ny, nz = sim.mesh.nx, sim.mesh.ny, sim.mesh.nz
  copy_spin_to_m(demag.mx, demag.my, demag.mz, spin, sim.Ms, nx, ny, nz)

  mul!(demag.Mx, demag.m_plan, demag.mx)
  mul!(demag.My, demag.m_plan, demag.my)
  mul!(demag.Mz, demag.m_plan, demag.mz)

  demag.Hx .= demag.tensor_xx.*demag.Mx .+ demag.tensor_xy.*demag.My .+  demag.tensor_xz.*demag.Mz
  demag.Hy .= demag.tensor_xy.*demag.Mx .+ demag.tensor_yy.*demag.My .+  demag.tensor_yz.*demag.Mz
  demag.Hz .= demag.tensor_xz.*demag.Mx .+ demag.tensor_yz.*demag.My .+  demag.tensor_zz.*demag.Mz

  mul!(demag.mx, demag.h_plan, demag.Hx) #we use mx to store field
  mul!(demag.my, demag.h_plan, demag.Hy)
  mul!(demag.mz, demag.h_plan, demag.Hz)

  extract_effective_field(demag.field, demag.mx, demag.my, demag.mz, nx, ny, nz)

  mu0 = 4*pi*1e-7
  volume = sim.mesh.volume
  field = demag.field
  for i=1:sim.n_nodes
    j = 3*i
    demag.energy[i] = -0.5*mu0*volume*sim.Ms[i]*(field[j-2]*spin[j-2] + field[j-1]*spin[j-1] + field[j]*spin[j])
  end
end


function newell_f(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;
   R = sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   f = 1.0/6*(2*x2-y2-z2)*R

   if x2>0
    f -= x*y*z*atan(y*z/(x*R));
   end

   if x2+z2>0
     f += 0.5*y*(z2-x2)*asinh(y/(sqrt(x2+z2)))
   end

   if x2+y2>0
     f += 0.5*z*(y2-x2)*asinh(z/(sqrt(x2+y2)))
   end
  return f;
end


function newell_g(x::Float64, y::Float64, z::Float64)::Float64
   x2 = x*x;
   y2 = y*y;
   z2 = z*z;

   R = sqrt(x2+y2+z2);
   if R == 0.0
     return 0.0;
   end

   g = -1.0/3*x*y*R;

   if z2>0   g -= 1.0/6*z2*z*atan(x*y/(z*R)) end
   if y2>0   g -= 0.5*y2*z*atan(x*z/(y*R)) end
   if x2>0   g -= 0.5*x2*z*atan(y*z/(x*R)) end

   if x2+y2>0
     g += x*y*z*asinh(z/(sqrt(x2+y2)))
    end

   if y2+z2>0
     g += 1.0/6*y*(3*z2-y2)*asinh(x/(sqrt(y2+z2)))
   end

   if x2+z2>0
     g += 1.0/6*x*(3*z2-x2)*asinh(y/(sqrt(x2+z2)))
   end

  return g;
end

#Numerical Micromagnetics: Finite Difference Methods, Jacques E. Miltat1 and Michael J. Donahue. Page 14.

function demag_tensor_xx(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)

  tensor = 8.0*newell_f(x,y,z);

  tensor -= 4.0*newell_f(x+dx,y,z);
  tensor -= 4.0*newell_f(x-dx,y,z);
  tensor -= 4.0*newell_f(x,y-dy,z);
  tensor -= 4.0*newell_f(x,y+dy,z);
  tensor -= 4.0*newell_f(x,y,z-dz);
  tensor -= 4.0*newell_f(x,y,z+dz);


  tensor +=  2.0*newell_f(x+dx,y+dy,z);
  tensor +=  2.0*newell_f(x+dx,y-dy,z);
  tensor +=  2.0*newell_f(x-dx,y-dy,z);
  tensor +=  2.0*newell_f(x-dx,y+dy,z);
  tensor +=  2.0*newell_f(x+dx,y,z+dz);
  tensor +=  2.0*newell_f(x+dx,y,z-dz);
  tensor +=  2.0*newell_f(x-dx,y,z+dz);
  tensor +=  2.0*newell_f(x-dx,y,z-dz);
  tensor +=  2.0*newell_f(x,y-dy,z-dz);
  tensor +=  2.0*newell_f(x,y-dy,z+dz);
  tensor +=  2.0*newell_f(x,y+dy,z+dz);
  tensor +=  2.0*newell_f(x,y+dy,z-dz);

  tensor -= newell_f(x+dx,y+dy,z+dz);
  tensor -= newell_f(x+dx,y+dy,z-dz);
  tensor -= newell_f(x+dx,y-dy,z+dz);
  tensor -= newell_f(x+dx,y-dy,z-dz);
  tensor -= newell_f(x-dx,y+dy,z+dz);
  tensor -= newell_f(x-dx,y+dy,z-dz);
  tensor -= newell_f(x-dx,y-dy,z+dz);
  tensor -= newell_f(x-dx,y-dy,z-dz);

  return tensor/(4.0*pi*dx*dy*dz)

end

function demag_tensor_yy(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xx(y,x,z,dy,dx,dz);
end

function demag_tensor_zz(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xx(z,y,x,dz,dy,dx);
end

function demag_tensor_xy(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)

  tensor = 8.0*newell_g(x,y,z);

  tensor -= 4.0*newell_g(x+dx,y,z);
  tensor -= 4.0*newell_g(x-dx,y,z);
  tensor -= 4.0*newell_g(x,y-dy,z);
  tensor -= 4.0*newell_g(x,y+dy,z);
  tensor -= 4.0*newell_g(x,y,z-dz);
  tensor -= 4.0*newell_g(x,y,z+dz);


  tensor +=  2.0*newell_g(x+dx,y+dy,z);
  tensor +=  2.0*newell_g(x+dx,y-dy,z);
  tensor +=  2.0*newell_g(x-dx,y-dy,z);
  tensor +=  2.0*newell_g(x-dx,y+dy,z);
  tensor +=  2.0*newell_g(x+dx,y,z+dz);
  tensor +=  2.0*newell_g(x+dx,y,z-dz);
  tensor +=  2.0*newell_g(x-dx,y,z+dz);
  tensor +=  2.0*newell_g(x-dx,y,z-dz);
  tensor +=  2.0*newell_g(x,y-dy,z-dz);
  tensor +=  2.0*newell_g(x,y-dy,z+dz);
  tensor +=  2.0*newell_g(x,y+dy,z+dz);
  tensor +=  2.0*newell_g(x,y+dy,z-dz);

  tensor -= newell_g(x+dx,y+dy,z+dz);
  tensor -= newell_g(x+dx,y+dy,z-dz);
  tensor -= newell_g(x+dx,y-dy,z+dz);
  tensor -= newell_g(x+dx,y-dy,z-dz);
  tensor -= newell_g(x-dx,y+dy,z+dz);
  tensor -= newell_g(x-dx,y+dy,z-dz);
  tensor -= newell_g(x-dx,y-dy,z+dz);
  tensor -= newell_g(x-dx,y-dy,z-dz);

  return tensor/(4.0*pi*dx*dy*dz)

end

function demag_tensor_xz(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xy(x,z,y,dx,dz,dy);
end

function demag_tensor_yz(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64, dz::Float64)
    return demag_tensor_xy(y,z,x,dy,dz,dx)
end

function fill_demag_tensors(long_tensor::Array{T,3}, tensor::Array{T,3}; tx=false, ty=false, tz=false) where {T<:Float64}
    lnx, lny, lnz = size(long_tensor)
    nx, ny, nz = size(tensor)
    for i=1:lnx, j=1:lny, k=1:lnz
        if (lnx%2 == 0 && i == nx+1) || (lny%2 == 0 && j == ny+1) || (lnz%2 == 0 && k == nz+1)
          continue
        end
        x = (i<=nx) ? i : lnx - i + 2
        y = (j<=ny) ? j : lny - j + 2
        z = (k<=nz) ? k : lnz - k + 2
        sx = tx && (i>nx) ? -1 : 1
        sy = ty && (j>ny) ? -1 : 1
        sz = tz && (k>nz) ? -1 : 1
        long_tensor[i,j,k] = sx*sy*sz*tensor[x,y,z]
    end
end

function compute_demag_tensors(tensor::Array{T,3}, tensor_fun::Function,  nx::I, ny::I, nz::I, Nx::I, Ny::I, Nz::I, dx::T, dy::T, dz::T) where {T<:Float64, I<:Int}
    for i = 1:nx, j=1:ny, k=1:nz
        tensor[i,j,k] = 0
        for p = -Nx:Nx, q=-Ny:Ny, s= -Nz:Nz
            x = (i-1 + p*nx)*dx
            y = (j-1 + q*ny)*dy
            z = (k-1 + s*nz)*dz
            tensor[i,j,k] += tensor_fun(x,y,z,dx,dy,dz)
        end
    end
end
