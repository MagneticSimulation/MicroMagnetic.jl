using NearestNeighbors
using Random
using StaticArrays
using PoissonDiskSampling

export voronoi, plot_voronoi

"""
    voronoi(mesh; min_dist=20, seed=123456, threshold=nothing, dims=3) -> grain_ids, gb_mask, points

Generate a Voronoi tessellation on a 2D/3D grid with grain boundaries detection using Poisson disk sampling.

# Arguments
- `mesh`: `FDMesh` object containing grid information
- `min_dist`: Minimum distance between Voronoi seeds in nanometers (default: 20)
- `seed`: Random number generator seed for reproducible results (default: 123456)
- `threshold`: Optional threshold for continuous boundary detection. If provided, uses distance-based boundary detection instead of discrete neighbor comparison.
- `dims`: Dimensionality of the tessellation (2 for 2D, 3 for 3D, default: 3)

# Returns
- `grain_ids`: Array of integers where each element represents the ID of the Voronoi cell it belongs to
- `gb_mask`: Boolean array marking grain boundary locations (true = boundary)
- `points`: Array of Voronoi seed points used for the tessellation

# Example
```julia
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=2)

# 3D Voronoi tessellation (default)
grain_ids, boundaries, seeds = voronoi(mesh; min_dist=20, seed=1000)

# 2D Voronoi tessellation
grain_ids, boundaries, seeds = voronoi(mesh; min_dist=20, seed=1000, dims=2)

# Continuous boundary detection with custom threshold
grain_ids, boundaries, seeds = voronoi(mesh; min_dist=20, seed=1000, threshold=0.2, dims=3)
```
"""
function voronoi(mesh; min_dist=20, seed=123456, threshold=nothing, dims=3)
    if dims == 2
        return voronoi2d(mesh; min_dist=min_dist, seed=seed, threshold=threshold)
    else
        return voronoi2d(mesh; min_dist=min_dist, seed=seed, threshold=threshold)
    end
end

function voronoi2d(mesh; min_dist=20, seed=123456, threshold=nothing)
    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx * 1e9, mesh.dy * 1e9  # Convert to nanometers
    Lx, Ly = nx * dx, ny * dy

    # Initialize random number generator
    rng = Xoshiro(seed)

    # Generate Poisson disk sampled seed points
    points = PoissonDiskSampling.generate(rng, min_dist, (0.0, Lx), (0.0, Ly))

    # Convert points to matrix format for KDTree
    points_matrix = hcat([collect(p) for p in points]...)

    # Build KDTree for efficient nearest neighbor queries
    kdtree = KDTree(points_matrix)

    # Create grid coordinates (cell centers)
    xs = [dx/2 + (i-1)*dx for i in 1:nx]
    ys = [dy/2 + (j-1)*dy for j in 1:ny]

    # Initialize output arrays
    grain_ids = zeros(Int, (nx, ny))
    gb_mask = falses((nx, ny))

    if threshold !== nothing
        # Distance-based approach
        boundary_strength = zeros(Float64, (nx, ny))
        D = sqrt(dx^2 + dy^2)
        for j in 1:ny, i in 1:nx
            pos = SVector(xs[i], ys[j])

            idxs, dists = NearestNeighbors.knn(kdtree, pos, 2)

            d1, d2 = dists[1], dists[2]
            grain_ids[i, j] = d1 < d2 ? idxs[1] : idxs[2]

            # Calculate boundary strength using exponential decay
            # The boundary is strongest when d1 ≈ d2
            boundary_strength[i, j] = exp(-abs(d2 - d1) / D)
            gb_mask[i, j] = boundary_strength[i, j] > threshold
        end
    else
        # Using neighbor comparison

        for j in 1:ny, i in 1:nx
            pos = SVector(xs[i], ys[j])
            idx, _ = NearestNeighbors.nn(kdtree, pos)
            grain_ids[i, j] = idx
        end

        for j in 1:ny, i in 1:nx
            if gb_mask[i, j]
                continue
            end

            if j<ny && grain_ids[i, j] != grain_ids[i, j+1]
                gb_mask[i, j] = true
            end

            if i<nx && grain_ids[i, j] != grain_ids[i+1, j]
                gb_mask[i, j] = true
            end

            if i<nx && j<ny && grain_ids[i, j] != grain_ids[i+1, j+1]
                gb_mask[i, j] = true
            end
        end
    end

    return grain_ids, gb_mask, points
end

function voronoi3d(mesh; min_dist=20, seed=123456, threshold=nothing)
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx * 1e9, mesh.dy * 1e9, mesh.dz * 1e9  # Convert to nanometers
    Lx, Ly, Lz = nx * dx, ny * dy, nz * dz

    # Initialize random number generator
    rng = Xoshiro(seed)

    # Generate Poisson disk sampled seed points
    points = PoissonDiskSampling.generate(rng, min_dist, (0.0, Lx), (0.0, Ly), (0.0, Lz))

    # Convert points to matrix format for KDTree
    points_matrix = hcat([collect(p) for p in points]...)

    # Build KDTree for efficient nearest neighbor queries
    kdtree = KDTree(points_matrix)

    # Create grid coordinates (cell centers)
    xs = [dx/2 + (i-1)*dx for i in 1:nx]
    ys = [dy/2 + (j-1)*dy for j in 1:ny]
    zs = [dz/2 + (k-1)*dz for k in 1:nz]

    # Initialize output arrays
    grain_ids = zeros(Int, (nx, ny, nz))
    gb_mask = falses((nx, ny, nz))

    if threshold !== nothing
        # Distance-based approach
        boundary_strength = zeros(Float64, (nx, ny, nz))
        D = sqrt(dx^2 + dy^2 + dz^2)
        for k in 1:nz, j in 1:ny, i in 1:nx
            pos = SVector(xs[i], ys[j], zs[k])

            idxs, dists = NearestNeighbors.knn(kdtree, pos, 2)

            d1, d2 = dists[1], dists[2]
            grain_ids[i, j, k] = d1 < d2 ? idxs[1] : idxs[2]

            # Calculate boundary strength using exponential decay
            # The boundary is strongest when d1 ≈ d2
            boundary_strength[i, j, k] = exp(-abs(d2 - d1) / D)
            gb_mask[i, j, k] = boundary_strength[i, j, k] > threshold
        end
    else
        # Using neighbor comparison
        for k in 1:nz, j in 1:ny, i in 1:nx
            pos = SVector(xs[i], ys[j], zs[k])
            idx, _ = NearestNeighbors.nn(kdtree, pos)
            grain_ids[i, j, k] = idx
        end

        for k in 1:nz, j in 1:ny, i in 1:nx
            if gb_mask[i, j, k]
                continue
            end

            if k<nz && grain_ids[i, j, k] != grain_ids[i, j, k+1]
                gb_mask[i, j, k] = true
            end

            if j<ny && grain_ids[i, j, k] != grain_ids[i, j+1, k]
                gb_mask[i, j, k] = true
            end

            if i<nx && grain_ids[i, j, k] != grain_ids[i+1, j, k]
                gb_mask[i, j, k] = true
            end

            if i<nx && j<ny && grain_ids[i, j, k] != grain_ids[i+1, j+1, k]
                gb_mask[i, j, k] = true
            end

            if i<nx && j<ny && k<nz && grain_ids[i, j, k] != grain_ids[i+1, j+1, k+1]
                gb_mask[i, j, k] = true
            end
        end
    end

    return grain_ids, gb_mask, points
end

function plot_voronoi end
