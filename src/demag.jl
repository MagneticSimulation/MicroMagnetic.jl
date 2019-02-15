using FFTW
using LinearAlgebra

mutable struct Demag
	nx_fft::Int64
  ny_fft::Int64
  nz_fft::Int64
  tensor_xx::Array{Float64}
  tensor_yy::Array{Float64}
  tensor_zz::Array{Float64}
  tensor_xy::Array{Float64}
  tensor_xz::Array{Float64}
  tensor_yz::Array{Float64}
  m_field::Array{Float64}
  Mx::Array{Complex{Float64}}
  My::Array{Complex{Float64}}
  Mz::Array{Complex{Float64}}
	h_field::Array{Float64}
  H_field::Array{Complex{Float64}}
  m_plan::Any
  h_plan::Any
  field::Array{Float64}
  energy::Array{Float64}
  name::String
end

function init_demag(sim::SimData)
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

  tensor_xx = zeros(nx_fft, ny_fft, nz_fft)
  tensor_yy = zeros(nx_fft, ny_fft, nz_fft)
  tensor_zz = zeros(nx_fft, ny_fft, nz_fft)
  tensor_xy = zeros(nx_fft, ny_fft, nz_fft)
  tensor_xz = zeros(nx_fft, ny_fft, nz_fft)
  tensor_yz = zeros(nx_fft, ny_fft, nz_fft)
  for i = 1:nx_fft
    for j = 1:ny_fft
      for k = 1:nz_fft
        if (nx_fft%2 == 0 && i == nx+1) || (ny_fft%2 == 0 && j == ny+1) || (nz_fft%2 == 0 && k == nz+1)
          continue
        end
        x = (i<=nx) ? (i-1)*dx : (i-nx_fft-1)*dx
        y = (j<=ny) ? (j-1)*dy : (j-ny_fft-1)*dy
        z = (k<=nz) ? (k-1)*dz : (k-nz_fft-1)*dz
        tensor_xx[i,j,k] = demag_tensor_xx(x,y,z,dx,dy,dz)
        tensor_yy[i,j,k] = demag_tensor_xx(y,x,z,dy,dx,dz)
        tensor_zz[i,j,k] = demag_tensor_xx(z,y,x,dz,dy,dx)
        tensor_xy[i,j,k] = demag_tensor_xy(x,y,z,dx,dy,dz)
        tensor_xz[i,j,k] = demag_tensor_xy(x,z,y,dx,dz,dy)
        tensor_yz[i,j,k] = demag_tensor_xy(y,z,x,dy,dz,dx)
      end
    end
  end

  tensor_xx = real(FFTW.rfft(tensor_xx)) #tensors xx, xy, .. should be pure real
  tensor_yy = real(FFTW.rfft(tensor_yy))
  tensor_zz = real(FFTW.rfft(tensor_zz))
  tensor_xy = real(FFTW.rfft(tensor_xy))
  tensor_xz = real(FFTW.rfft(tensor_xz))
  tensor_yz = real(FFTW.rfft(tensor_yz))

	m_field = zeros(nx_fft, ny_fft, nz_fft)
  lenx = (nx_fft%2>0) ? nx : nx+1
  Mx = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  My = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  Mz = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)
  m_plan = FFTW.plan_rfft(m_field)

	h_field = zeros(nx_fft, ny_fft, nz_fft)
	H_field = zeros(Complex{Float64}, lenx, ny_fft, nz_fft)

  h_plan = FFTW.plan_irfft(H_field, nx_fft)

	field = zeros(Float64, 3*sim.nxyz)
	energy = zeros(Float64, sim.nxyz)
  demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz,
                tensor_xy, tensor_xz, tensor_yz, m_field, Mx, My, Mz, h_field,
                H_field, m_plan, h_plan, field, energy, "Demag")
  return demag
end

function copy_spin_to_m(m::Array{Float64}, spin::Array{Float64}, Ms::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, bias::Int64)
  for k=1:nz, j=1:ny, i=1:nx
        p = (k-1) * nx*ny + (j-1) * nx + i
        m[i,j,k] = spin[3*p+bias]*Ms[p]
  end
end

function copy_cfield_to_field(field::Array{Float64}, cfield::Array{Float64}, nx::Int64, ny::Int64, nz::Int64, bias::Int64)
  for k=1:nz, j=1:ny, i=1:nx
    p = (k-1) * nx*ny + (j-1) * nx + i
    field[3*p+bias] = -1.0*cfield[i,j,k]
  end
end

function effective_field(demag::Demag, sim::SimData, spin::Array{Float64}, t::Float64)
  nx, ny, nz = sim.mesh.nx, sim.mesh.ny, sim.mesh.nz
  copy_spin_to_m(demag.m_field, spin, sim.Ms, nx, ny, nz, -2)
  mul!(demag.Mx, demag.m_plan, demag.m_field)
	copy_spin_to_m(demag.m_field, spin, sim.Ms, nx, ny, nz, -1)
	mul!(demag.My, demag.m_plan, demag.m_field)
	copy_spin_to_m(demag.m_field, spin, sim.Ms, nx, ny, nz, 0)
	mul!(demag.Mz, demag.m_plan, demag.m_field)

  demag.H_field .= demag.tensor_xx.*demag.Mx .+ demag.tensor_xy.*demag.My .+  demag.tensor_xz.*demag.Mz
  mul!(demag.h_field, demag.h_plan, demag.H_field)
  copy_cfield_to_field(demag.field, demag.h_field, nx, ny, nz, -2)

	demag.H_field .= demag.tensor_xy.*demag.Mx .+ demag.tensor_yy.*demag.My .+  demag.tensor_yz.*demag.Mz
  mul!(demag.h_field, demag.h_plan, demag.H_field)
  copy_cfield_to_field(demag.field, demag.h_field, nx, ny, nz, -1)

	demag.H_field .= demag.tensor_xz.*demag.Mx .+ demag.tensor_yz.*demag.My .+  demag.tensor_zz.*demag.Mz
  mul!(demag.h_field, demag.h_plan, demag.H_field)
  copy_cfield_to_field(demag.field, demag.h_field, nx, ny, nz, 0)

  for i=1:sim.nxyz
    j = 3*i
    demag.energy[i] = -0.5*sim.Ms[i]*(demag.field[j-2]*spin[j-2] + demag.field[j-1]*spin[j-1] + demag.field[j]*spin[j])
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
