using NearestNeighbors
using Random
using StaticArrays 
using PoissonDiskSampling

export voronoi, plot_voronoi

"""
    voronoi(mesh; min_dist=20, seed=123456) -> grain_ids, gb_mask, points

Generate a Voronoi tessellation on a 2D grid with grain boundaries detection.

# Arguments
- `mesh`: `FDMesh`
- `min_dist`: Minimum distance between Voronoi seeds in nanometers (default: 20)
- `seed`: Random number generator seed for reproducible results (default: 123456)

# Returns
- `grain_ids`: Matrix of integers where each element represents the ID of the Voronoi cell it belongs to
- `gb_mask`: Boolean matrix marking grain boundary locations (true = boundary)
- `points`: Array of Voronoi seed points used for the tessellation

# Example
```julia
mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=2)
grain_ids, boundaries, seeds = voronoi(mesh; min_dist=40, seed=1000)
```
"""
function voronoi(mesh; min_dist = 20, seed=123456)
    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx*1e9, mesh.dy*1e9
    Lx, Ly = nx*dx, ny*dy

    rng = Xoshiro(seed);

    points = PoissonDiskSampling.generate(
        rng, min_dist, 
        (0, Lx), (0, Ly)
    );
    points_matrix = hcat([collect(p) for p in points]...)

    kdtree = KDTree(points_matrix)
    
    xs = [dx/2 + i*dx for i in 1:nx]
    ys = [dy/2 + j*dy for j in 1:ny]
     
    grain_ids = zeros(Int, (nx, ny))
    gb_mask = falses((nx, ny))

    for j in 1:ny, i in 1:nx
        pos = SVector(xs[i], ys[j])
        idx, _ = NearestNeighbors.nn(kdtree, pos)
        grain_ids[i,j] = idx
    end

    for j in 1:ny, i in 1:nx
        if gb_mask[i,j] 
            continue
        end

        if j<ny && grain_ids[i,j] != grain_ids[i, j+1] 
            gb_mask[i,j] = true
        end

        if i<nx && grain_ids[i,j] != grain_ids[i+1, j] 
            gb_mask[i,j] = true
        end

        if i<nx && j<ny && grain_ids[i,j] != grain_ids[i+1, j+1]
            gb_mask[i,j] = true
        end

    end

    return grain_ids, gb_mask, points
end


function plot_voronoi() end