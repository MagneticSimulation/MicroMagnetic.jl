mutable struct NEB_SD <: NEBDriver
  gk::Array{Float64, 2}
  tau::Float64
  max_tau::Float64
  min_tau::Float64
  nsteps::Int64
end

function compute_tau(driver::NEB_SD, pre_images::Array{Float64, 2}, images::Array{Float64, 2}, h::Array{Float64, 2}, nxyz::Int64,N)

  if driver. == 0
    for n=1:N, i=0:nxyz-1  #compute gk[] for step 0
        j = 3*i+1
        fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
        gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
        driver.gk[j,n] = gx
        driver.gk[j+1,n] = gy
        driver.gk[j+2,n] = gz
      end
     driver.tau  = driver.min_tau
     return nothing
  end
  sum1, sum2, sum3 = zeros(N),zeros(N),zeros(N)
  for n=1:N
    for i=0:nxyz-1
        j = 3*i+1
        sx = images[j,n] - pre_images[j,n]
        sy = images[j+1,n] - pre_images[j+1,n]
        sz = images[j+2,n] - pre_images[j+2,n]
        fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
        gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
        fx = gx - driver.gk[j,n]
        fy = gy - driver.gk[j+1,n]
        fz = gz - driver.gk[j+2,n]
        driver.gk[j,n] = gx
        driver.gk[j+1,n] = gy
        driver.gk[j+2,n] = gz
        sum1[n] += sx*sx+sy*sy+sz*sz
        sum2[n] += sx*fx+sy*fy+sz*fz
        sum3[n] += fx*fx+fy*fy+fz*fz
     end
  end
     sum_1,sum_2,sum_3 = sum(sum1)/N,sum(sum2)/N,sum(sum3)/N
     tau1 = sum_2!=0.0 ? sum_1/sum_2 : driver.min_tau
     tau2 = sum_3!=0.0 ? sum_2/sum_3 : driver.min_tau
     driver.tau = driver.%2 == 0 ? abs(tau2) : abs(tau1)
     if driver.tau > driver.max_tau
        driver.tau = driver.max_tau
     end
     return nothing
end

function update_images(driver::NEB_SD,images::Array{Float64,2},h::Array{Float64,2},nxyz::Int64,N::Int64)
  tau = driver.tau
  new_images = zeros(3*nxyz,N)
  for n =1:N
    if (n == 1)||(n==N)
      for i=0:nxyz-1
        j = 3*i+1
        new_images[j,n] = images[j,n]
        new_images[j+1,n] = images[j+1,n]
        new_images[j+2,n] = images[j+2,n]
      end
    end
    for i=0:nxyz-1
      j = 3*i+1
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
      gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
      factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
      #driver.gk[j,n],driver.gk[j+1,n],driver.gk[j+2,n] = gx,gy,gz
      mx = (1-factor)*images[j,n] - tau*gx
      my = (1-factor)*images[j+1,n] - tau*gy
      mz = (1-factor)*images[j+2,n] - tau*gz
      new_images[j,n] = mx/(1+factor)
      new_images[j+1,n] = my/(1+factor)
      new_images[j+2,n] = mz/(1+factor)
    end
  end
  for n=1:N
    max_length_error = error_length_m(new_images[:,n], nxyz)
    if max_length_error > 1e-15
      new_images[:,n] = normalized(new_images[:,n], nxyz)                                 #fix me
    end
  end
  return new_images
end

function compute_norm_gk(gk::Array{Float64,2},N::Int)
  sum = 0.0
  for n =1:N
    if (n==0)||(n==N)
      continue
    end
    sum += model(gk[:,n])
  end
  return (sum/(N-2))
end

function compute_gk(images::Array{Float64,2},h::Array{Float64,2},nxyz::Int64,N::Int64)
  gk = zeros(3*nxyz,N)
  for n =1:N
    if (n == 1)||(n==N)
      continue
    end
    for i=0:nxyz-1
      j = 3*i+1
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
      gk[j,n],gk[j+1,n],gk[j+2,n] = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
    end
  end
  return gk
end

function search_tau(neb::NEB)
  sim = neb.sim
  nxyz = sim.nxyz
  N = neb.N
  driver = neb.driver
  images=zeros(N,3*nxyz)
  h0=zeros(N,3*nxyz)
  images[:] = neb.images[:]
  h0[:] = neb.field[:]
  for n=1:N, i=0:nxyz-1  #compute gk[] for step 0
    j = 3*i+1
    fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h0[j,n],h0[j+1,n],h0[j+2,n])
    gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
    driver.gk[j,n] = gx
    driver.gk[j+1,n] = gy
    driver.gk[j+2,n] = gz
  end
  tau1 = driver.min_tau
  dtau = 0.1*driver.min_tau
  tau = driver.min_tau
  tau2 = driver.min_tau
  driver.tau = tau1
  neb.images = update_images(driver,images,h0,nxyz,N)

  effective_field_NEB(neb, neb.spin)

  gk1 =compute_gk(neb.images,neb.field,nxyz,N)
  gk1_norm = compute_norm_gk(gk1,N)
  for i = 0:100
    tau2 = tau1 + dtau
    driver.tau = tau2
    neb.images = update_images(driver,images,h0,nxyz,N)
    effective_field_NEB(neb, neb.spin)
    gk2 =compute_gk(neb.images,neb.field,nxyz,N)
    gk2_norm = compute_norm_gk(gk2,N)
    if gk1_norm>=gk2_norm
      dtau = 2*dtau
      tau = tau1
      gk1_norm = gk2_norm
      tau1 = tau2
    elseif i==0
      dtau = -dtau
    else
      break
    end
  end
  a = min(tau,tau2)
  b = max(tau,tau2)
  for j=1:10
    tau1 = a + 1/3*(b-a)
    driver.tau = tau1
    neb.images = update_images(driver,images,h0,nxyz,N)
    effective_field_NEB(neb, neb.spin)
    gk1 =compute_gk(neb.images,neb.field,nxyz,N)
    gk1_norm = compute_norm_gk(gk1,N)

    tau2 = a + 2/3*(b-a)
    driver.tau = tau2
    neb.images = update_images(driver,images,h0,nxyz,N)
    effective_field_NEB(neb, neb.spin)
    gk2 =compute_gk(neb.images,neb.field,nxyz,N)
    gk2_norm = compute_norm_gk(gk2,N)
    if gk1_norm<gk2_norm
      a = tau1
    else
      b = tau2
    end
    if b-a<1e-14
      break
    end
  end
  driver.tau = a
  neb.field[:] = h0[:]
  neb.images[:] = images[:]
end


function run_step(neb::NEB)
  driver = neb.driver
  N = neb.N
  images = neb.images
  pre_images = neb.pre_images
  sim = neb.sim
  nxyz = sim.nxyz
  images=zeros(N,3*nxyz)
  effective_field_NEB(neb, neb.spin, 0.0)
  if driver.nsteps == 0
    search_tau(neb)
    #driver.tau = 1e-12
    #compute_tau(driver, pre_images, images, neb.field, nxyz, N)
  else
    compute_tau(driver, pre_images, images, neb.field, nxyz, N)
  end
  neb.pre_images[:] =  neb.images[:]
  h = neb.field
  gk = driver.gk
  tau = driver.tau
  for n =1:N
    if (n == 1)||(n==N)
      continue
    end
    for i=0:sim.nxyz-1
      j = 3*i+1
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
      factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
      mx = (1-factor)*images[j,n] - tau*gk[j,n]
      my = (1-factor)*images[j+1,n] - tau*gk[j+1,n]
      mz = (1-factor)*images[j+2,n] - tau*gk[j+2,n]
      neb.images[j,n] = mx/(1+factor)
      neb.images[j+1,n] = my/(1+factor)
      neb.images[j+2,n] = mz/(1+factor)
    end
  end
  for n=1:N
    max_length_error = error_length_m(neb.images[:,n], sim.nxyz)
    if max_length_error > 1e-15
      neb.images[:,n] = normalized(neb.images[:,n], sim.nxyz)
    end
  end
  driver.nsteps += 1
  return  nothing
end

function normalized(a::Array{T, 1}, N::Int64) where {T<:AbstractFloat}
  for i = 0:N-1
      j = 3*i+1
      length = sqrt(a[j]*a[j]+a[j+1]*a[j+1]+a[j+2]*a[j+2])
      if length > 0
        length = 1.0/length
        a[j] *= length
        a[j+1] *= length
        a[j+2] *= length
      end
    end
   return a
end

function relax_NEB(neb::NEB, maxsteps::Int64, stopping_torque::Float64, save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String,save_ovf_every::Int64,ovf_folder::String)
  N = neb.N
  sim = neb.sim
  gk_abs = zeros(Float64,3*sim.nxyz)
  driver = neb.driver
  maxtorque=zeros(N)
  if save_m_every>0
    compute_system_energy(neb)
    write_data(neb)
  end
  if save_vtk_every > 0
    save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, 0)))
  end
  if save_ovf_every > 0
    save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, 0)))
  end
  for i=1:maxsteps
      run_step(neb)
      #maxtorque = zeros(N)
      for n = 1:N
        if (n==1)||(n==N)
          continue
        end
        abs!(gk_abs, driver.gk[:,n])  #max_torque = maximum(abs.(driver.gk)) eats gpu memory???
        maxtorque[n] = maximum(gk_abs)
      end

      max_torque = maximum(maxtorque)
      @info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e", i, driver.tau, max_torque)

      if save_m_every>0 && i%save_m_every == 0
        compute_system_energy(neb)
        write_data(neb)
      end
      if save_vtk_every > 0 && i%save_vtk_every == 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
      end
      if save_ovf_every > 0 && i%save_ovf_every == 0
        save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
      end
      if max_torque < stopping_torque
        @info @sprintf("max_torque (mxmxH) is less than stopping_torque=%g, Done!", stopping_torque)
        if save_m_every>0
         compute_system_energy(neb)
          write_data(neb)
        end
        if save_vtk_every > 0
          save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
        end
        if save_ovf_every > 0
          save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
        end
      break
      end
  end
return nothing
end
