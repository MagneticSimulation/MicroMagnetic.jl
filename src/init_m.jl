using Random

export init_m0_random, vortex, skyrmion, bubble2d, skyrmion_lattice

"""
    init_m0_random(sim::MicroSim; seed=nothing)

Set the initial magnetization with random direction.

# Keyword Arguments
- `seed::Union{Integer, Nothing}=nothing`: Random seed for reproducible results.
  If `nothing` (default), uses the global random number generator.
"""
function init_m0_random(sim::AbstractSim; seed=nothing)
    n_total = sim.n_total

    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    
    m_data = 2 .* rand(rng, 3 * n_total) .- 1
    
    init_m0(sim, m_data)
end

"""
    vortex(;center=(0,0), R=10e-9, p=1, c=1)

Create a vortex function that can be invoked from `init_m0`.

Parameters:
- center: Vortex center coordinates
- R: Vortex radius (default: 10e-9)
- p: Polarity (+1 for up, -1 for down)
- c: Chirality (+1 for counterclockwise, -1 for clockwise)

Example:
```julia
init_m0(sim, vortex(p=1, c=1))
```
"""
function vortex(;center=(0,0), R=10e-9, p=1, c=1)
    cx, cy = center[1], center[2] 
    function vortex_fun(x, y, z)
        dx = x - cx
        dy = y - cy
        r = sqrt(dx^2 + dy^2)
    
        if r < R
            return (0,0,p)
        end

        theta = atan(dy, dx)
        return (-c * sin(theta),c * cos(theta),0)
    end
    return vortex_fun
end


"""
    bubble2d(; center=(0,0), R=10e-9, w=5e-9, p=1, v=1, phi=0.0)

Generate a bubble profile function for topological magnetic textures (skyrmions, bubbles, etc.).

Mathematical formulation:
    m(r,ϕ) = (sinΘ(r) * cosΦ(ϕ), sinΘ(r) * sinΦ(ϕ), p * cosΘ(r))
    
    where:
    - Θ(r) = 2 * atan(exp((r - R) / w))  (polar angle profile)
    - Φ(ϕ) = v * ϕ + phi             (azimuthal angle profile)
    
    Topological charge: Q = p * v

Parameters:
    center      :: Tuple{Float64,Float64}  - (x,y) center position of the texture
    R           :: Float64                 - Bubble radius (distance where Θ = π/2)
    width       :: Float64                 - Domain wall width (controls transition sharpness)
    p           :: Int                     - Polarity (+1: up at center, -1: down at center)
    v           :: Int                     - Vorticity/winding number (S in literature)
    phi         :: Float64                 - Helicity angle (radians) - phase offset

Returns:
    A function f(x,y,z) -> (m_x, m_y, m_z) that returns normalized magnetization vector at position (x,y,z)

Examples:
```julia
    bubble = bubble2d(p=1, v=1, phi=0, w=5e-9, R=40e-9)
    init_m0(sim, bubble)
```
"""
function bubble2d(; center=(0.0,0.0), R=10e-9, w=5e-9, p=1, v=1, phi=0.0)
    
    # Extract center coordinates    
    cx, cy = center[1], center[2]
    
    function texture_function(x, y, z)
        # Calculate relative position
        dx = x - cx
        dy = y - cy
        r = sqrt(dx^2 + dy^2)      # Radial distance from center
        ϕ = atan(dy, dx)           # Azimuthal angle
        
        # Calculate polar angle Θ(r) using domain wall profile
        # This gives a smooth transition from Θ=0 at r=0 to Θ=π at r→∞
        if w > 0.0
            u = p * (r - R) / w
            Θ = 2.0 * atan(exp(u))
        else
            Θ = (r < R) ? 0.0 : π
        end
                
        # This determines the in-plane magnetization rotation
        Φ = v * ϕ + phi
        sinΘ = sin(Θ)
        cosΘ = cos(Θ)        
        m_x = sinΘ * cos(Φ)
        m_y = sinΘ * sin(Φ)
        m_z = cosΘ
        
        return (m_x, m_y, m_z)
    end
    
    return texture_function
end

"""
    skyrmion(;center=(0,0), R=10e-9, p=1, v=1, phi=0.0)

Create a skyrmion function that can be invoked from `init_m0`.

Parameters:
    center      :: Tuple{Float64,Float64}  - (x,y) center position of the texture
    R           :: Float64                 - Skyrmion radius (distance where Θ = π/2)
    p           :: Int                     - Polarity (+1: up at center, -1: down at center)
    v           :: Int                     - Vorticity/winding number (S in literature)
    phi         :: Float64                 - Helicity angle (radians) - phase offset

Example:
```julia
    # Néel skyrmion with upward center
    neel_skyrmion = skyrmion(p=1, v=1, phi=0.0)
    
    # Bloch skyrmion  
    bloch_skyrmion = skyrmion(p=1, v=1, phi=pi/2)
    
    # Antiskyrmion
    antiskyrmion = skyrmion(p=1, v=-1, phi=0.0)

    init_m0(sim, antiskyrmion)
```
"""
function skyrmion(center=(0.0,0.0), R=10e-9, p=1, v=1, phi=0.0)
    return bubble2d(center=center, R=R, w=R/2, p=p, v=v, phi=phi)
end


# Helper function for 2D rotation
function rotation_2d(angle_deg::Float64, vector::Vector{Float64})
    θ = deg2rad(angle_deg)
    cosθ, sinθ = cos(θ), sin(θ)
    x, y = vector[1], vector[2]
    return [cosθ*x - sinθ*y, sinθ*x + cosθ*y]
end


"""
    skyrmion_lattice(lambda::Float64; p=-1, c=+1, type=:bloch)

Create a function that returns magnetization for a skyrmion lattice.

Parameters:
    lambda    :: Float64     - Skyrmion lattice period (m)
    p         :: Int         - polarity: +1 (up) or -1 (down)
    c         :: Int         - chirality - +1 (clockwise) or -1 (counterclockwise)
    type      :: Symbol      - :bloch or :neel

Returns:
    A function f(x, y, z) -> (m_x, m_y, m_z) that returns normalized magnetization.

Examples:
```julia
    lattice = skyrmion_lattice(50e-9, p=-1, c=1, type=:bloch)
    init_m0(sim, lattice)
```
"""
function skyrmion_lattice(lambda::Float64; p=-1, c=1, type=:bloch)
    
    # Wave vectors for triangular lattice (120° apart)
    Q1 = [2π / lambda, 0.0]           # Primary direction
    Q2 = rotation_2d(120.0, Q1)       # Rotated 120°
    Q3 = rotation_2d(120.0, Q2)       # Rotated 240°
    
    function lattice_function(x, y, z)
        # Calculate phases for each Q-vector
        ϕ1 = Q1[1]*x + Q1[2]*y
        ϕ2 = Q2[1]*x + Q2[2]*y 
        ϕ3 = Q3[1]*x + Q3[2]*y 
        
        # Calculate in-plane components
        mx = sin(ϕ1) + sin(ϕ2 + 2π/3) + sin(ϕ3 + 4π/3)
        my = cos(ϕ1) + cos(ϕ2 + 2π/3) + cos(ϕ3 + 4π/3)
        
        # Apply chirality
        if c == -1
            mx, my = rotation_2d(180.0, [mx, my])
        end
        
        if type == :bloch
            # Rotate by 90° for Bloch type
            mx, my = rotation_2d(90.0, [mx, my])
        end
        
        mz = cos(ϕ1) + cos(ϕ2) + cos(ϕ3)
        
        norm = p*sqrt(mx^2 + my^2 + mz^2)
        return (mx/norm, my/norm, mz/norm)
    end
    
    return lattice_function
end