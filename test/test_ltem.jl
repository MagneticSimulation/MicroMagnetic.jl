using MicroMagnetic
using Test

# ---------------------------------------------------------------------------
# warp2d / rotate3d
# ---------------------------------------------------------------------------

function test_warp2d()
    # a 90-degree rotation must move a pixel-centered pulse to the rotated
    # position exactly (even sizes map centers onto centers)
    img = zeros(16, 16)
    img[10, 10] = 1.0
    w = MicroMagnetic.warp2d(img, pi / 2)
    @test w[7, 10] ≈ 1.0 atol = 1e-12
    @test sum(w) ≈ 1.0 atol = 1e-12

    # a constant field is invariant in the interior
    img = ones(16, 16)
    w = MicroMagnetic.warp2d(img, 0.3)
    @test all(w[6:11, 6:11] .≈ 1.0)
    @test w[1, 1] == 0.0   # corners sample outside the array
end

function test_rotate3d()
    m = zeros(3, 8, 8, 8)
    m[1, :, :, :] .= 1.0                    # constant field along x
    angles = [(0.3, 0.0, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 0.7), (0.2, 0.3, 0.1)]
    for (tx, ty, tz) in angles
        mr = MicroMagnetic.rotate3d(m, tx, ty, tz)
        R = MicroMagnetic.euler_matrix(tx, ty, tz)
        # interior voxel: m_new(r) = R * m_old(R^{-1} r) = R * (1,0,0)
        for c in 1:3
            @test mr[c, 4, 4, 4] ≈ R[c, 1] atol = 1e-10
        end
    end

    # single-axis rotations are exactly invertible: rotate + rotate back
    # reproduces a smooth field (up to interpolation error).  The field is
    # smooth, long-wavelength and compact so that the residual is dominated
    # by the bilinear interpolation, not by singularities or clipping.
    v = zeros(3, 16, 16, 16)
    for k in 1:16, j in 1:16, i in 1:16
        env = exp(-((i - 8.5)^2 + (j - 8.5)^2 + (k - 8.5)^2) / 18.0)
        v[1, i, j, k] = sin(2 * pi * i / 32) * cos(2 * pi * j / 48) * env
        v[2, i, j, k] = cos(2 * pi * i / 48) * sin(2 * pi * j / 32) * env
        v[3, i, j, k] = sin(2 * pi * k / 32) * env
    end
    for angles in [(0.3, 0.0, 0.0), (0.0, -0.25, 0.0), (0.0, 0.0, 0.4)]
        v1 = MicroMagnetic.rotate3d(v, angles...)
        v2 = MicroMagnetic.rotate3d(v1, (.-angles)...)
        @test maximum(abs.(v2 .- v)) / maximum(abs.(v)) < 5e-2
    end
end

# ---------------------------------------------------------------------------
# padding / grid size
# ---------------------------------------------------------------------------

function test_padding()
    @test MicroMagnetic.padding_range(5, 9) == 3:7
    @test MicroMagnetic.padding_range(5, 10) == 4:8
    @test MicroMagnetic.padding_range(12, 16) == 3:14

    v = rand(3, 5, 6, 7)
    p = MicroMagnetic.vector_padding(v, 16)
    @test size(p) == (3, 16, 16, 16)
    r1 = MicroMagnetic.padding_range(5, 16)
    r2 = MicroMagnetic.padding_range(6, 16)
    r3 = MicroMagnetic.padding_range(7, 16)
    @test p[:, r1, r2, r3] == v

    # grid size contains the rotated bounding box
    @test MicroMagnetic.rotation_grid_size(12, 12, 6, 0.0, 0.0, 0.0) == 16
    @test MicroMagnetic.rotation_grid_size(8, 8, 4, pi / 2, 0.0, 0.0) == 8
    @test MicroMagnetic.rotation_grid_size(12, 12, 6, 0.0, 0.4, 0.0) >= 16
end

# ---------------------------------------------------------------------------
# projection
# ---------------------------------------------------------------------------

function test_projection()
    # uniform out-of-plane magnetization, beam along z
    m = zeros(3, 8, 8, 4)
    m[3, :, :, :] .= 1.0
    p = MicroMagnetic.project3d(m)
    @test size(p) == (3, 8, 8)
    @test all(p[3, :, :] .≈ 4.0) && all(p[1, :, :] .≈ 0) && all(p[2, :, :] .≈ 0)

    # rotate by 90 degrees about x: mz -> -my (active rotation convention);
    # the box's 8-cell x dimension maps onto z, so the vertical chord is 8,
    # and the 4-cell z dimension maps onto y (footprint y in 3:6)
    vp = MicroMagnetic.vector_field_projection(m, 90.0, "x")
    @test size(vp) == (3, 8, 8)
    @test vp[2, 4, 4] ≈ -8.0 atol = 1e-12
    @test all(abs.(vp[2, :, 3:6] .+ 8.0) .< 1e-12)
    @test vp[1, 4, 4] == 0.0 && abs(vp[3, 4, 4]) < 1e-12
    @test vp[2, 1, 1] == 0.0

    # rotate by 90 degrees about y: mz -> mx; the padded box is centered in
    # the cube, so the rotated box (4, 8, 8) spans x in 3:6 with chord 8
    vp = MicroMagnetic.vector_field_projection(m, 90.0, "y")
    @test vp[1, 4, 4] ≈ 8.0 atol = 1e-12
    @test all(abs.(vp[1, 3:6, :] .- 8.0) .< 1e-12)
    @test vp[1, 1, 1] == 0.0
    @test maximum(abs.(vp[2, :, :])) < 1e-12

    # a 45-degree tilt: beam-integrated in-plane induction is nz (the longer
    # path exactly cancels the cos-projected component); no y component
    vp = MicroMagnetic.vector_field_projection(m, 45.0, "y")
    @test size(vp) == (3, 16, 16)
    @test maximum(vp[1, :, :]) ≈ 4.0 atol = 0.2
    @test maximum(abs.(vp[2, :, :])) < 1e-12
end

# ---------------------------------------------------------------------------
# magnetic phase: FFT kernel vs slow direct reference
# ---------------------------------------------------------------------------

function vortex_m(nx::Int, ny::Int, nz::Int)
    m = zeros(3, nx, ny, nz)
    for i in 1:nx, j in 1:ny
        phi = atan(j - (ny + 1) / 2, i - (nx + 1) / 2)
        env = exp(-(((i - (nx + 1) / 2)^2 + (j - (ny + 1) / 2)^2) / 25.0)^4)
        m[1, i, j, :] .= -sin(phi) * env
        m[2, i, j, :] .= cos(phi) * env
    end
    return m
end

function test_phase_direct_vs_fft()
    dx = dy = dz = 4e-9
    Ms = 8e5
    m = vortex_m(24, 24, 4)
    Lz = 4 * dz

    phi = compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz)
    @test size(phi) == (32, 32)                # nextpow2(24)

    # phase is odd under magnetization reversal (consistency check)
    phi_rev = compute_magnetic_phase(-m; Ms=Ms, dx=dx, dy=dy, dz=dz)
    scale = maximum(abs.(phi))
    @test maximum(abs.(phi .+ phi_rev)) < 1e-12 * scale

    # cross-check against the direct real-space kernel at interior points;
    # the 24-cell sample sits at padding_range(24, 32) = 5:28 on the output grid
    off = first(MicroMagnetic.padding_range(24, 32)) - 1
    mx2d, my2d = m[1, :, :, 1], m[2, :, :, 1]
    worst = 0.0
    for (i0, j0) in [(8, 12), (14, 9), (11, 16), (16, 14), (9, 15)]
        phi_d = MicroMagnetic.compute_magnetic_phase_direct(mx2d, my2d, dx, dy, Lz, Ms, i0, j0)
        worst = max(worst, abs(phi[i0+off, j0+off] - phi_d) / scale)
    end
    @test worst < 5e-3
end

# ---------------------------------------------------------------------------
# tilt: rigid-rotation covariance of the phase
# ---------------------------------------------------------------------------

function bilinear(img::AbstractMatrix, x::Float64, y::Float64)
    i1, j1 = floor(Int, x), floor(Int, y)
    fx, fy = x - i1, y - j1
    n = size(img, 1)
    (1 <= i1 < n && 1 <= j1 < n) || return NaN
    return (1 - fx) * (1 - fy) * img[i1, j1] + fx * (1 - fy) * img[i1+1, j1] +
           (1 - fx) * fy * img[i1, j1+1] + fx * fy * img[i1+1, j1+1]
end

function test_phase_tilt_covariance()
    # Rotating the (sample + beam) system rigidly must leave the phase
    # invariant as a scalar field:  phi_alpha(r) = phi_0(R^{-1} r).
    # The residual is dominated by the O(alpha^2) bilinear resampling error
    # of the tilted volume (even in alpha, verified separately), so the
    # tolerance is loose; the sharp physical discrimination of the tilt
    # convention is done in test_phase_tilt_perpendicular.
    dx = dy = dz = 4e-9
    Ms = 8e5
    m = vortex_m(24, 24, 4)
    N = 32
    alpha = 0.4
    phi0 = compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz, N=N)
    phiA = compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz, ty=alpha, N=N)

    c = (N + 1) / 2
    R = MicroMagnetic.euler_matrix(0.0, -alpha, 0.0)   # back-rotation
    scale = maximum(abs.(phi0))
    worst = 0.0
    for (i, j) in [(14, 16), (17, 15), (15, 18), (12, 14), (18, 17)]
        p = R * [i - c, j - c, 0.0]
        ref = bilinear(phi0, p[1] + c, p[2] + c)
        worst = max(worst, abs(phiA[i, j] - ref) / scale)
    end
    @test worst < 0.15
end

function test_phase_tilt_perpendicular()
    # Tilted perpendicularly magnetized film: the beam sees the in-plane
    # component Ms*sin(tilt) integrated along the longer path t/cos(tilt),
    # so the phase must scale like tan(alpha) and must be NONZERO.
    # (Rotating the projected vector back into the sample frame — as
    # MagRecon/maglab does — cancels these two factors exactly and would
    # give zero contrast at every tilt angle, contradicting experiment.)
    dx = dy = dz = 4e-9
    Ms = 8e5
    m = zeros(3, 32, 32, 3)
    m[3, :, :, :] .= 1.0

    phi0 = compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz, ty=0.0, N=48)
    @test maximum(abs.(phi0)) == 0.0        # untilted: no in-plane induction

    p1 = maximum(abs.(compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz, ty=0.2, N=48)))
    p2 = maximum(abs.(compute_magnetic_phase(m; Ms=Ms, dx=dx, dy=dy, dz=dz, ty=0.4, N=48)))
    @test p1 > 0.01                          # nonzero contrast at tilt
    @test p2 / p1 ≈ tan(0.4) / tan(0.2) rtol = 0.05

    # absolute scale: interior phase ramp ~ (e/hbar)*mu0*Ms*t*tan(alpha)*L,
    # reduced by edge/demag effects; a 2pi- or mu0-level coefficient bug
    # would push this out of the bracket
    t = 3 * dz
    ramp = (pi / MicroMagnetic.magnetic_flux_quantum) * mu_0 * Ms * t * tan(0.2) * (32 * dx)
    @test 0.15 < p1 / ramp < 0.6
end

# ---------------------------------------------------------------------------
# LTEM imaging
# ---------------------------------------------------------------------------

function test_ltem()
    dx = dy = dz = 5e-9
    Ms = 8e5
    m = vortex_m(24, 24, 4)

    # at zero defocus the image of a pure phase object is flat: |exp(i phi)|^2 = 1
    phi, img = LTEM(m; Ms=Ms, dx=dx, dy=dy, dz=dz, df=0, V0=0)
    @test size(img) == size(phi) == (32, 32)
    @test img ≈ ones(32, 32) atol = 1e-10

    # defocus produces contrast; phase is real and finite
    phi, img = LTEM(m; Ms=Ms, dx=dx, dy=dy, dz=dz, df=200, V0=0)
    @test maximum(img) - minimum(img) > 1e-3
    @test all(0 .<= img .<= 4)
    @test all(isfinite, phi)

    # tilted sample still yields a valid image
    phi_t, img_t = LTEM(m; Ms=Ms, dx=dx, dy=dy, dz=dz, df=200, V0=0, ty=deg2rad(30))
    @test size(img_t) == (32, 32)
    @test maximum(img_t) - minimum(img_t) > 1e-3
end

# ---------------------------------------------------------------------------
# ovf interface + electron optics helpers
# ---------------------------------------------------------------------------

function test_ovf_path()
    m = vortex_m(16, 16, 4)
    ovf = MicroMagnetic.OVF2{Float64}()
    ovf.xnodes, ovf.ynodes, ovf.znodes = 16, 16, 4
    ovf.xstepsize, ovf.ystepsize, ovf.zstepsize = 4e-9, 4e-9, 4e-9
    ovf.data = vec(m)

    phi_ovf = compute_magnetic_phase(ovf; Ms=8e5)
    phi_arr = compute_magnetic_phase(m; Ms=8e5, dx=4e-9, dy=4e-9, dz=4e-9)
    @test phi_ovf == phi_arr

    # tilting via the positional (theta, axis) wrapper agrees with the kwarg
    phi_w = compute_magnetic_phase(m, deg2rad(20), "y"; Ms=8e5, dx=4e-9, dy=4e-9, dz=4e-9)
    phi_k = compute_magnetic_phase(m; Ms=8e5, dx=4e-9, dy=4e-9, dz=4e-9, ty=deg2rad(20))
    @test phi_w ≈ phi_k atol = 1e-12

    # electron optics helpers
    lambda = MicroMagnetic.electron_wavelength(300e3)
    @test isapprox(lambda, 1.9687e-12, rtol=1e-3)          # 300 kV
    @test isapprox(MicroMagnetic.interaction_constant(300e3), 6.53e6, rtol=1e-2)
    l2, phi_E = MicroMagnetic.compute_electric_phase(300e3, 20.0, 20e-9, 0.0)
    @test l2 == lambda
    @test isapprox(phi_E, 6.53e6 * 20 * 20e-9, rtol=1e-2)
end

@testset "tools/warp" begin
    test_warp2d()
end

@testset "tools/rotate" begin
    test_rotate3d()
end

@testset "tools/padding" begin
    test_padding()
end

@testset "tools/projection" begin
    test_projection()
end

@testset "tools/magnetic phase" begin
    test_phase_direct_vs_fft()
end

@testset "tools/tilt covariance" begin
    test_phase_tilt_covariance()
    test_phase_tilt_perpendicular()
end

@testset "tools/LTEM" begin
    test_ltem()
end

@testset "tools/ovf path" begin
    test_ovf_path()
end
