using KernelAbstractions

function newell_f(x::Float64, y::Float64, z::Float64)::Float64
    x2 = x * x
    y2 = y * y
    z2 = z * z
    R = sqrt(x2 + y2 + z2)
    if R == 0.0
        return 0.0
    end

    f = 1.0 / 6 * (2 * x2 - y2 - z2) * R
    if x2 > 0
        f -= x * y * z * atan(y * z / (x * R))
    end

    if x2 + z2 > 0
        f += 0.5 * y * (z2 - x2) * asinh(y / (sqrt(x2 + z2)))
    end

    if x2 + y2 > 0
        f += 0.5 * z * (y2 - x2) * asinh(z / (sqrt(x2 + y2)))
    end
    return f
end

function demag_tensor_xx(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64,
                         dz::Float64)
    tensor = 8.0 * newell_f(x, y, z)

    tensor -= 4.0 * newell_f(x + dx, y, z)
    tensor -= 4.0 * newell_f(x - dx, y, z)
    tensor -= 4.0 * newell_f(x, y - dy, z)
    tensor -= 4.0 * newell_f(x, y + dy, z)
    tensor -= 4.0 * newell_f(x, y, z - dz)
    tensor -= 4.0 * newell_f(x, y, z + dz)

    tensor += 2.0 * newell_f(x + dx, y + dy, z)
    tensor += 2.0 * newell_f(x + dx, y - dy, z)
    tensor += 2.0 * newell_f(x - dx, y - dy, z)
    tensor += 2.0 * newell_f(x - dx, y + dy, z)
    tensor += 2.0 * newell_f(x + dx, y, z + dz)
    tensor += 2.0 * newell_f(x + dx, y, z - dz)
    tensor += 2.0 * newell_f(x - dx, y, z + dz)
    tensor += 2.0 * newell_f(x - dx, y, z - dz)
    tensor += 2.0 * newell_f(x, y - dy, z - dz)
    tensor += 2.0 * newell_f(x, y - dy, z + dz)
    tensor += 2.0 * newell_f(x, y + dy, z + dz)
    tensor += 2.0 * newell_f(x, y + dy, z - dz)

    tensor -= newell_f(x + dx, y + dy, z + dz)
    tensor -= newell_f(x + dx, y + dy, z - dz)
    tensor -= newell_f(x + dx, y - dy, z + dz)
    tensor -= newell_f(x + dx, y - dy, z - dz)
    tensor -= newell_f(x - dx, y + dy, z + dz)
    tensor -= newell_f(x - dx, y + dy, z - dz)
    tensor -= newell_f(x - dx, y - dy, z + dz)
    tensor -= newell_f(x - dx, y - dy, z - dz)

    return tensor / (4.0 * pi * dx * dy * dz)
end

function naive_transpose_kernel!(a, b)
    i, j = @index(Global, NTuple)
    @inbounds b[i, j] = a[j, i]
end

@kernel function compute_tensors_kernel_xx!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum = sum + demag_tensor_xx(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

using CUDA
nx, ny, nz = 10, 5, 7
dx, dy, dz = 1.0, 1.0, 1.0
Nx, Ny, Nz = 3, 1, 1
tensor = CUDA.zeros(Float64, nx, ny, nz)

groupsize = 512
kernel! = compute_tensors_kernel_xx!(get_backend(tensor), groupsize)
kernel!(tensor, dx, dy, dz, Nx, Ny, Nz; ndrange=size(tensor))
KernelAbstractions.synchronize(get_backend(tensor))

print(Array(tensor))
