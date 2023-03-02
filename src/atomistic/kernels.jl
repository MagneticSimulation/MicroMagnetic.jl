using CUDA

function __magnetoelectric_laser__kernel_001!(m, h, energy, lambda, Ex, Ey, Ez, Hx, Hy, Hz, mu_s, N)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0 < i <= N
        j = 3*(i-1)

        @inbounds ms_local = mu_s[i]
        if ms_local == 0.0
            @inbounds energy[i] = 0
            @inbounds h[j+1] = 0
            @inbounds h[j+2] = 0
            @inbounds h[j+3] = 0
            return nothing
        end
        ms_inv = lambda/ms_local
        @inbounds mx = m[j+1]
        @inbounds my = m[j+2]
        @inbounds mz = m[j+3]

        @inbounds h[j+1] = -(Ez*mx + Ex*mz)*ms_inv
        @inbounds h[j+2] = (Ez*my + Ey*mz)*ms_inv 
        @inbounds h[j+3] = (-Ex*mx + Ey*my)*ms_inv
        
        #compute the energy for the electric part
        @inbounds energy[i] = -0.5*ms_local*(mx*h[j+1] + my*h[j+2] + mz*h[j+3])

        #add the contribution of the magnetic part
        @inbounds h[j+1] += Hx
        @inbounds h[j+2] += Hy
        @inbounds h[j+3] += Hz

        @inbounds energy[i] -= ms_local*(mx*Hx + my*Hy + mz*Hz)
    end

    return nothing
end



function __magnetoelectric_laser__kernel_110!(m, h, energy, lambda, Ex, Ey, Ez, Hx, Hy, Hz, mu_s, N)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0 < i <= N
        j = 3*(i-1)

        @inbounds ms_local = mu_s[i]
        if ms_local == 0.0
            @inbounds energy[i] = 0
            @inbounds h[j+1] = 0
            @inbounds h[j+2] = 0
            @inbounds h[j+3] = 0
            return nothing
        end
        ms_inv = lambda/ms_local
        @inbounds mx = m[j+1]
        @inbounds my = m[j+2]
        @inbounds mz = m[j+3]

        @inbounds h[j+1] = -(Ey*mx + Ex*my)*ms_inv
        @inbounds h[j+2] = (-Ex*mx + Ez*mz)*ms_inv 
        @inbounds h[j+3] = (Ez*my + Ey*mz)*ms_inv
        
        #compute the energy for the electric part
        @inbounds energy[i] = -0.5*ms_local*(mx*h[j+1] + my*h[j+2] + mz*h[j+3])

        #add the contribution of the magnetic part
        @inbounds h[j+1] += Hx
        @inbounds h[j+2] += Hy
        @inbounds h[j+3] += Hz

        @inbounds energy[i] -= ms_local*(mx*Hx + my*Hy + mz*Hz)
    end

    return nothing
end

function __magnetoelectric_laser__kernel_111!(m, h, energy, lambda, Ex, Ey, Ez, Hx, Hy, Hz, mu_s, N)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if 0 < i <= N
        j = 3*(i-1)

        @inbounds ms_local = mu_s[i]
        if ms_local == 0.0
            @inbounds energy[i] = 0
            @inbounds h[j+1] = 0
            @inbounds h[j+2] = 0
            @inbounds h[j+3] = 0
            return nothing
        end
        ms_inv = -lambda/ms_local/CUDA.sqrt(3.0)
        @inbounds mx = m[j+1]
        @inbounds my = m[j+2]
        @inbounds mz = m[j+3]

        @inbounds h[j+1] = (CUDA.sqrt(2.0)*(Ey*mx + Ex*my)+Ez*mx+Ex*mz)*ms_inv
        @inbounds h[j+2] = (CUDA.sqrt(2.0)*(Ex*mx - Ey*my)+Ez*my+Ey*mz)*ms_inv 
        @inbounds h[j+3] = (Ex*mx + Ey*my-2*Ez*mz)*ms_inv
        
        #compute the energy for the electric part
        @inbounds energy[i] = -0.5*ms_local*(mx*h[j+1] + my*h[j+2] + mz*h[j+3])

        #add the contribution of the magnetic part
        @inbounds h[j+1] += Hx
        @inbounds h[j+2] += Hy
        @inbounds h[j+3] += Hz

        @inbounds energy[i] -= ms_local*(mx*Hx + my*Hy + mz*Hz)
    end

    return nothing
end