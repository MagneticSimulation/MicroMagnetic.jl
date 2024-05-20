using Test
using NuMag

warp = NuMag.warp
radon = NuMag.radon

function test_warp()
    a=zeros(4,4)
    a[2,2] = 1
    a[2,3] = 2
    a[3,2] = 3
    a[3,3] = 4

    b = warp(a, pi/2)
    @test isapprox(b[2,2], 2)
    @test isapprox(b[2,3], 4)
    @test isapprox(b[3,2], 1)
    @test isapprox(b[3,3], 3)
end

function test_radon()
    a = zeros(4,4)
    a[2,2] = 1
    a[3,2] = 10
    a[2,3] = 100
    a[3,3] = 1000

    sinogram = radon(a, 0.0)
    @test isapprox(sinogram, [0, 11, 1100, 0.])

    sinogram = radon(a, pi/2)
    @test isapprox(sinogram, [0., 1010., 101., 0.])

    sinogram = radon(a, pi*1)
    @test isapprox(sinogram, [0., 1100., 11., 0.])

    sinogram = radon(a, -pi/2)
    @test isapprox(sinogram,[0., 101., 1010., 0.])
end

function test_radon_3d()
    a = zeros(4,4,4)
    a[2,2,2] = 1e0
    a[3,2,2] = 1e1
    a[2,3,2] = 1e2
    a[3,3,2] = 1e3

    a[2,2,3] = 1e4
    a[3,2,3] = 1e5
    a[2,3,3] = 1e6
    a[3,3,3] = 1e7

    projection = NuMag.radon_x(a, 0.0)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.0001e4 1.0001e6 0.;
            0. 1.0001e5 1.0001e7 0.;
            0. 0. 0. 0.]
    )


    projection = NuMag.radon_x(a, pi/2)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.01e6 1.01e2 0.;
            0. 1.01e7 1.01e3 0.;
            0. 0. 0. 0.]
    )

    projection = NuMag.radon_x(a, -pi/2)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.01e2 1.01e6 0.;
            0. 1.01e3 1.01e7 0.;
            0. 0. 0. 0.]
    )

    projection = NuMag.radon_y(a, 0.0)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.0001e4 1.0001e6 0.;
            0. 1.0001e5 1.0001e7 0.;
            0. 0. 0. 0.]
    )


    projection = NuMag.radon_y(a, pi/2)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.1e1 1.1e3 0.;
            0. 1.1e5 1.1e7 0.;
            0. 0. 0. 0.]
    )

    projection = NuMag.radon_y(a, -pi/2)
    @test isapprox(
        projection, 
        [   0. 0. 0. 0.;
            0. 1.1e5 1.1e7 0.;
            0. 1.1e1 1.1e3 0.;
            0. 0. 0. 0.]
    )
end

function test_vector_field_projection()
    a = zeros(3,4,4,4)
    for i = 2:3, j= 2:3, k=2:3
        a[1,i,j,k] = 1
        a[2,i,j,k] = 2
        a[3,i,j,k] = 3
    end

    gamma = pi/3
    projection = NuMag.radon_x(a, gamma)
    vy = NuMag.radon_x(a[2,:,:,:], gamma)
    vz = NuMag.radon_x(a[3,:,:,:], gamma)
    @test isapprox(
        projection, -1 .* (cos(gamma).*vz + sin(gamma).*vy)
    )

    beta = pi/4
    projection = NuMag.radon_y(a, beta)
    vx = NuMag.radon_y(a[1,:,:,:], beta)
    vz = NuMag.radon_y(a[3,:,:,:], beta)
    @test isapprox(
        projection, -1 .* (cos(beta).*vz - sin(beta).*vx)
    )
end

test_warp()
test_radon()
test_radon_3d()
test_vector_field_projection()
