function uniform_random_sphere_kernel!(m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        j = 3*index - 2
        @inbounds phi = rnd[j]*2*pi;
        @inbounds ct = 2*rnd[j+1] - 1;
        st = CUDA.sqrt(1-ct*ct);
        @inbounds m[j] = st*CUDA.cos(phi);
        @inbounds m[j+1] = st*CUDA.sin(phi);
        @inbounds m[j+2] = ct;
    end
   return nothing
end

function uniform_random_circle_xy_kernel!(m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        j = 3*index - 2
        @inbounds phi = rnd[j]*2*pi;
        @inbounds m[j] = CUDA.cos(phi);
        @inbounds m[j+1] = CUDA.sin(phi);
        @inbounds m[j+2] = 0;
    end
   return nothing
end

function dE_zeeman_anisotropy_cubic_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              Hx::T, Hy::T, Hz::T, Ku::T, Kc::T,
                              ux::T, uy::T, uz::T,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #bias should be 0, 1 and 2 for cubic mesh
    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if (a+b+c)%3 != bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman

        @inbounds delta_E += Ku*CUDA.pow(m[i]*ux+m[i+1]*uy+m[i+2]*uz, 2) #Anisotropy
        @inbounds delta_E -= Ku*CUDA.pow(next_m[i]*ux+next_m[i+1]*uy+next_m[i+2]*uz, 2)
        @inbounds delta_E += Kc*(CUDA.pow(m[i+2],4) - CUDA.pow(next_m[i+2],4)) #Cubic Anisotropy for z
        @inbounds delta_E += Kc*(CUDA.pow(m[i+1],4) - CUDA.pow(next_m[i+1],4)) #Cubic Anisotropy for y
        @inbounds delta_E += Kc*(CUDA.pow(m[i],4) - CUDA.pow(next_m[i],4)) #Cubic Anisotropy for x
        @inbounds dE[index] = delta_E
   end
   return nothing
end

function total_E_zeeman_anisotropy_kernel!(m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1},
                              Hx::T, Hy::T, Hz::T, Ku::T, Kc::T,
                              ux::T, uy::T, uz::T,
                              nxyz::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #bias should be 0, 1 and 2 for cubic mesh
    if index <= nxyz
        @inbounds si = shape[index]
        if !si
            @inbounds dE[index] = 0
            return nothing
        end

        i = 3*index - 2
        @inbounds mx =  m[i]
        @inbounds my =  m[i+1]
        @inbounds mz = m[i+2]

        delta_E = -(mx*Hx+my*Hy+mz*Hz) #zeeman

        delta_E -= Ku*CUDA.pow(mx*ux+my*uy+mz*uz, 2) #Anisotropy

        delta_E += Kc*(CUDA.pow(mx,4)+CUDA.pow(my,4)+CUDA.pow(mz,4))

        @inbounds dE[index] = delta_E
   end
   return nothing
end

function dE_zeeman_anisotropy_triangular_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              Hx::T, Hy::T, Hz::T, Ku::T,
                              ux::T, uy::T, uz::T,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #bias should be 0, 1 and 2 for cubic mesh
    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if mod(a-b+c,3)!= bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman

        @inbounds delta_E += Ku*CUDA.pow(m[i]*ux+m[i+1]*uy+m[i+2]*uz, 2) #Anisotropy
        @inbounds delta_E -= Ku*CUDA.pow(next_m[i]*ux+next_m[i+1]*uy+next_m[i+2]*uz, 2)
        @inbounds dE[index] = delta_E
   end
   return nothing
end


function dE_zeeman_kagome_anisotropy_triangular_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              Hx::T, Hy::T, Hz::T, Ku::T,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #bias should be 0, 1 and 2 for cubic mesh
    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if mod(a-b+c,3)!= bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        ux,uy,uz = T(0),T(0),T(0)
        if a%2==1 && b%2==1
            ux,uy,uz = (-0.5,-sqrt(3)/2,0)
        elseif a%2==0 && b%2==1
            ux,uy,uz = (1,0,0)
        elseif a%2==1 && b%2==0
            ux,uy,uz = (-0.5,sqrt(3)/2,0)
        end


        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman

        @inbounds delta_E += Ku*CUDA.pow(m[i]*ux+m[i+1]*uy+m[i+2]*uz, 2) #Anisotropy
        @inbounds delta_E -= Ku*CUDA.pow(next_m[i]*ux+next_m[i+1]*uy+next_m[i+2]*uz, 2)
        @inbounds dE[index] = delta_E
   end
   return nothing
end

function dE_zeeman_kagome_anisotropy_6fold_triangular_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              Hx::T, Hy::T, Hz::T, K1::T, K2::T,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #bias should be 0, 1 and 2 for cubic mesh
    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if mod(a-b+c,3)!= bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        # ux,uy,uz = T(0),T(0),T(0)
        # if a%2==1 && b%2==1
        #     ux,uy,uz = (-0.5,-sqrt(3)/2,0)
        # elseif a%2==0 && b%2==1
        #     ux,uy,uz = (1,0,0)
        # elseif a%2==1 && b%2==0
        #     ux,uy,uz = (-0.5,sqrt(3)/2,0)
        # end
        u1x,u1y = T(0),T(0),T(0)
        u2x,u2y = T(0),T(0),T(0)
        u3x,u3y = T(0),T(0),T(0)
        u1x,u1y = (1,0)
        u2x,u2y = (-0.5,sqrt(3)/2)
        u3x,u3y = (-0.5,-sqrt(3)/2)
        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman

        # @inbounds delta_E += Ku*CUDA.pow(m[i]*ux+m[i+1]*uy+m[i+2]*uz, 2) #Anisotropy
        # @inbounds delta_E -= Ku*CUDA.pow(next_m[i]*ux+next_m[i+1]*uy+next_m[i+2]*uz, 2)
        @inbounds mu1s=CUDA.pow(m[i]*u1x+m[i+1]*u1y,2)
        @inbounds mu2s=CUDA.pow(m[i]*u2x+m[i+1]*u2y,2)
        @inbounds mu3s=CUDA.pow(m[i]*u3x+m[i+1]*u3y,2)
        delta_E += K1*(mu1s*mu2s+mu2s*mu3s+mu3s*mu1s)+K2*(mu1s*mu2s*mu3s)
        @inbounds nu1s=CUDA.pow(next_m[i]*u1x+next_m[i+1]*u1y,2)
        @inbounds nu2s=CUDA.pow(next_m[i]*u2x+next_m[i+1]*u2y,2)
        @inbounds nu3s=CUDA.pow(next_m[i]*u3x+next_m[i+1]*u3y,2)
        delta_E -= K1*(nu1s*nu2s+nu2s*nu3s+nu3s*nu1s)+K2*(nu1s*nu2s*nu3s)

        @inbounds dE[index] = delta_E
   end
   return nothing
end

#compute exchange and dmi energy for cubic mesh, including the next-nearest correction, can run in parallel
function add_dE_exch_dmi_cubic_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              nngbs::CuDeviceArray{Int32, 2}, n_ngbs::Int64, J0::CuDeviceArray{T, 1},
                              J1::CuDeviceArray{T, 1}, D0::CuDeviceArray{T, 2}, D1::CuDeviceArray{T, 2},
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if (a+b+c)%3 != bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = 0.0

        for j = 1:n_ngbs
            @inbounds id = ngbs[j, index]

            @inbounds J = J0[j]
            @inbounds Dx = D0[1, j]
            @inbounds Dy = D0[2, j]
            @inbounds Dz = D0[3, j]

            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) #exchange
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end


            @inbounds J = J1[j]
            @inbounds Dx = D1[1, j]
            @inbounds Dy = D1[2, j]
            @inbounds Dz = D1[3, j]

            @inbounds id = nngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) #exchange
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end

        end
        @inbounds dE[index] += delta_E
    end
   return nothing
end


function add_total_E_exch_dmi_cubic_mesh_kernel!(m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              nngbs::CuDeviceArray{Int32, 2}, n_ngbs::Int64, J0::CuDeviceArray{T, 1},
                              J1::CuDeviceArray{T, 1}, D0::CuDeviceArray{T, 2}, D1::CuDeviceArray{T, 2},
                              nxyz::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if index <= nxyz

        @inbounds si = shape[index]
        if !si
            @inbounds dE[index] = 0
            return nothing
        end

        i = 3*index - 2
        @inbounds mx = m[i]
        @inbounds my = m[i+1]
        @inbounds mz = m[i+2]

        delta_E = 0.0

        for j = 1:n_ngbs
            @inbounds id = ngbs[j, index]

            @inbounds J = J0[j]
            @inbounds Dx = D0[1, j]
            @inbounds Dy = D0[2, j]
            @inbounds Dz = D0[3, j]

            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(mx*sx + my*sy + mz*sz) #exchange
                delta_E += volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end


            @inbounds J = J1[j]
            @inbounds Dx = D1[1, j]
            @inbounds Dy = D1[2, j]
            @inbounds Dz = D1[3, j]

            @inbounds id = nngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(mx*sx + my*sy + mz*sz) #exchange
                delta_E += volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end

        end
        @inbounds dE[index] += 0.5*delta_E
    end
   return nothing
end


#compute nearest exchange and dmi energy for triangular mesh, can run in parallel
function add_dE_exch_dmi_triangular_mesh_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              n_ngbs::Int64, J0::CuDeviceArray{T, 1},
                              D0::CuDeviceArray{T, 2},
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        if mod(a-b+c,3)!= bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = 0.0

        for j = 1:n_ngbs
            @inbounds id = ngbs[j, index]

            @inbounds J = J0[j]
            @inbounds Dx = D0[1, j]
            @inbounds Dy = D0[2, j]
            @inbounds Dz = D0[3, j]

            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) #exchange
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end
        @inbounds dE[index] += delta_E
    end
   return nothing
end


function add_total_E_exch_dmi_triangular_mesh_kernel!(m::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              dE::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              n_ngbs::Int64, J0::CuDeviceArray{T, 1},
                              D0::CuDeviceArray{T, 2},
                              nxyz::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if index <= nxyz

        @inbounds si = shape[index]
        if !si
            @inbounds dE[index] = 0
            return nothing
        end

        i = 3*index - 2
        @inbounds mx = m[i]
        @inbounds my = m[i+1]
        @inbounds mz = m[i+2]

        delta_E = 0.0

        for j = 1:n_ngbs
            @inbounds id = ngbs[j, index]

            @inbounds J = J0[j]
            @inbounds Dx = D0[1, j]
            @inbounds Dy = D0[2, j]
            @inbounds Dz = D0[3, j]

            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(mx*sx + my*sy + mz*sz) #exchange
                delta_E += volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end
        @inbounds dE[index] += 0.5*delta_E
    end
   return nothing
end


function run_monte_carlo_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1},
                              shape::CuDeviceArray{Bool, 1},
                              energy::CuDeviceArray{T, 1}, temp::Float64,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64, cubic::Bool) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nx*ny*nz && shape[index]
        a, b, c = Tuple(CUDA.CartesianIndices((nx,ny,nz))[index])
        sign = cubic ? 1 : -1
        if mod(a+sign*b+c, 3) != bias
            return nothing
        end

        i = 3*index - 2
        update = false
        delta_E = T(0)

        @inbounds delta_E = energy[index]

        if delta_E < 0
            update = true
        else
            @inbounds rnd_i = rnd[i+2]
            if rnd_i < CUDA.exp(-delta_E/temp)
                update = true
            end
        end

        if update == true
          @inbounds m[i] = next_m[i];
          @inbounds m[i+1] = next_m[i+1];
          @inbounds m[i+2] = next_m[i+2]
        end

    end
   return nothing
end
