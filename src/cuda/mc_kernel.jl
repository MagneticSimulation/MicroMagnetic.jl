

function run_step_triangular_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1},
                              energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              J::T, lambda::T, D::T, Dz::T, bulk_dmi::Bool, Ku::T, Hx::T, Hy::T, Hz::T, temp::T,
                              nx::Int64, ny::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (T(1), T(1/2), T(-1/2), T(-1), T(-1/2), T(1/2))
    Ry = (T(0), T(sqrt(3)/2), T(sqrt(3)/2), T(0), T(-sqrt(3)/2), T(-sqrt(3)/2))

    delta_E = T(0.0)
    update = false
    #bias should be 0, 1 and 2.

    if index <= nx*ny
        a, b = Tuple(CuArrays.CartesianIndices((nx,ny))[index])
        if (a+b)%3 != bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman
        @inbounds delta_E += Ku*(m[i+2]*m[i+2] - next_m[i+2]*next_m[i+2]) #Anisotropy

        for j = 1:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) + lambda*dmz*sz #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                Dzt = j%2 == 0 ? Dz : -1*Dz
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dzt) #DMI
            end
        end
        if delta_E < 0
            update = true
        else
            @inbounds rnd_i = rnd[i+2]
            if rnd_i < CUDAnative.exp(-delta_E/temp)
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

function run_step_cubic_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1},
                              energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              nngbs::CuDeviceArray{Int32, 2}, J::T, Jz::T, J1::T, D::T, Dz::T, D1::T, bulk_dmi::Bool,
                              Ku::T, Kc::T, Hx::T, Hy::T, Hz::T, temp::T,
                              nx::Int64, ny::Int64, nz::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (1, -1, 0, 0, 0, 0)
    Ry = (0, 0, 1, -1, 0, 0)
    Rz = (0, 0, 0, 0, 1, -1)

    delta_E = T(0.0)
    update = false
    #bias should be 0, 1 and 2.

    if index <= nx*ny*nz
        a, b, c = Tuple(CuArrays.CartesianIndices((nx,ny,nz))[index])
        if (a+b+c)%3 != bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman
        @inbounds delta_E += Ku*(m[i+2]*m[i+2] - next_m[i+2]*next_m[i+2]) #Anisotropy
        @inbounds delta_E += Kc*(CUDAnative.pow(m[i+2],4) - CUDAnative.pow(next_m[i+2],4)) #Cubic Anisotropy for z
        @inbounds delta_E += Kc*(CUDAnative.pow(m[i+1],4) - CUDAnative.pow(next_m[i+1],4)) #Cubic Anisotropy for y
        @inbounds delta_E += Kc*(CUDAnative.pow(m[i],4) - CUDAnative.pow(next_m[i],4)) #Cubic Anisotropy for x

        for j = 1:4
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                Dz = bulk_dmi ? Dz*Rz[j] : T(0)
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end

        for j = 5:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= Jz*(dmx*sx + dmy*sy + dmz*sz) #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                Dz = bulk_dmi ? Dz*Rz[j] : T(0)
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end

        for j = 1:6
            id = nngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J1*(dmx*sx + dmy*sy + dmz*sz) #exchange

                Dx = bulk_dmi ? D1*Rx[j] : -D1*Ry[j]
                Dy = bulk_dmi ? D1*Ry[j] : D1*Rx[j]
                Dz = bulk_dmi ? D1*Rz[j] : T(0)
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end
        end
        if delta_E < 0
            update = true
        else
            @inbounds rnd_i = rnd[i+2]
            if rnd_i < CUDAnative.exp(-delta_E/temp)
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


function total_energy_triangular_kernel!(m::CuDeviceArray{T, 1}, energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              J::T, lambda::T, D::T, Dz::T, bulk_dmi::Bool, Ku::T, Hx::T, Hy::T, Hz::T,
                              nxy::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (T(1), T(1/2), T(-1/2), T(-1), T(-1/2), T(1/2))
    Ry = (T(0), T(sqrt(3)/2), T(sqrt(3)/2), T(0), T(-sqrt(3)/2), T(-sqrt(3)/2))

    E = T(0.0)

    if index <= nxy
        i = 3*index - 2
        @inbounds mx = m[i]
        @inbounds my = m[i+1]
        @inbounds mz = m[i+2]

        E = -(mx*Hx+my*Hy+mz*Hz) #zeeman
        E -= Ku*(mz*mz) #Anisotropy

        for j = 1:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                E -= 0.5*J*(mx*sx + my*sy + mz*sz) + 0.5*lambda*mz*sz #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]  #\hat{z} \times r_ij
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                Dzt = j%2 == 0 ? Dz : -1*Dz
                E += 0.5*volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dzt) #DMI
            end
        end

        @inbounds energy[index] = E

    end
   return nothing
end

function total_energy_cubic_kernel!(m::CuDeviceArray{T, 1}, energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              nngbs::CuDeviceArray{Int32, 2},
                              J::T, J1::T, D::T, D1::T, bulk_dmi::Bool, Ku::T, Hx::T, Hy::T, Hz::T,
                              nxyz::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (1, -1, 0, 0, 0, 0)
    Ry = (0, 0, 1, -1, 0, 0)
    Rz = (0, 0, 0, 0, 1,-1)

    E = T(0.0)

    if index <= nxyz
        i = 3*index - 2
        @inbounds mx = m[i]
        @inbounds my = m[i+1]
        @inbounds mz = m[i+2]

        E = -(mx*Hx+my*Hy+mz*Hz) #zeeman
        E -= Ku*(mz*mz) #Anisotropy

        for j = 1:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                E -= 0.5*J*(mx*sx + my*sy + mz*sz) #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                Dz = bulk_dmi ? D*Rz[j] : T(0)
                E += 0.5*volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end

            id = nngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                E -= J1*(mx*sx + my*sy + mz*sz) #exchange

                Dx = bulk_dmi ? D1*Rx[j] : -D1*Ry[j]
                Dy = bulk_dmi ? D1*Ry[j] : D1*Rx[j]
                Dz = bulk_dmi ? D1*Rz[j] : T(0)
                E += volume(mx, my, mz, sx, sy, sz, Dx, Dy, Dz) #DMI
            end

        end

        @inbounds energy[index] = E

    end
   return nothing
end
