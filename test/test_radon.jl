using Test
using JuMag

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

    sinogram = radon(a, 0)
    @test isapprox(sinogram[2,1], 101)
    @test isapprox(sinogram[3,1], 1010)

    sinogram = radon(a, 90)
    @test isapprox(sinogram[2,1], 11)
    @test isapprox(sinogram[3,1], 1100)

    sinogram = radon(a, 180)
    @test isapprox(sinogram[2,1], 1010)
    @test isapprox(sinogram[3,1], 101)

    sinogram = radon(a, -90)
    @test isapprox(sinogram[2,1], 1100)
    @test isapprox(sinogram[3,1], 11)
end

function test_radon_3d()
    a = zeros(4,4,4)
    a[2,2,2] = 1
    a[3,2,2] = 10
    a[2,3,2] = 100
    a[3,3,2] = 1000

    a[2,2,3] = 10000
    a[3,2,3] = 100000
    a[2,3,3] = 1000000
    a[3,3,3] = 10000000

    sinogram = radon_3d(a, [0, 90, -90], "x")
    @test isapprox(sinogram[2,2,1], 10001)
    @test isapprox(sinogram[3,2,1], 100010)
    @test isapprox(sinogram[2,3,1], 1000100)
    @test isapprox(sinogram[3,3,1], 10001000)

    @test isapprox(sinogram[2,2,2], 101)
    @test isapprox(sinogram[3,2,2], 1010)
    @test isapprox(sinogram[2,3,2], 1010000)
    @test isapprox(sinogram[3,3,2], 10100000)

    @test isapprox(sinogram[2,2,3], 1010000)
    @test isapprox(sinogram[3,2,3], 10100000)
    @test isapprox(sinogram[2,3,3], 101)
    @test isapprox(sinogram[3,3,3], 1010)

    sinogram = radon_3d(a, [0, 90, -90], "y")
    @test isapprox(sinogram[2,2,1], 10001)
    @test isapprox(sinogram[3,2,1], 100010)
    @test isapprox(sinogram[2,3,1], 1000100)
    @test isapprox(sinogram[3,3,1], 10001000)

    @test isapprox(sinogram[2,2,2], 110000)
    @test isapprox(sinogram[3,2,2], 11)
    @test isapprox(sinogram[2,3,2], 11000000)
    @test isapprox(sinogram[3,3,2], 1100)

    @test isapprox(sinogram[2,2,3], 11)
    @test isapprox(sinogram[3,2,3], 110000)
    @test isapprox(sinogram[2,3,3], 1100)
    @test isapprox(sinogram[3,3,3], 11000000)
end

function test_vector_field_projection()
    a = zeros(3,4,4,4)
    for i = 2:3, j= 2:3, k=2:3
        a[1,i,j,k] = 1
        a[2,i,j,k] = 2
        a[3,i,j,k] = 3
    end

    p1 = vector_field_projection(a, 0, "x")
    @test isapprox(p1[2,2,1], -6)
    @test isapprox(p1[2,3,1], -6)
    @test isapprox(p1[3,2,1], -6)
    @test isapprox(p1[3,3,1], -6)

    p1 = vector_field_projection(a, 90, "x")
    @test isapprox(p1[2,2,1], 4)
    @test isapprox(p1[2,3,1], 4)
    @test isapprox(p1[3,2,1], 4)
    @test isapprox(p1[3,3,1], 4)

    p1 = vector_field_projection(a, 90, "y")
    @test isapprox(p1[2,2,1], -2)
    @test isapprox(p1[2,3,1], -2)
    @test isapprox(p1[3,2,1], -2)
    @test isapprox(p1[3,3,1], -2)
end

function test_vector_field_projectionx()
    a = zeros(3,20,20,20)
    vx, vy, vz = 2, 3, 5
    for i = 6:15,j=6:15,k=6:15
        a[1,i,j,k] = vx * i
        a[2,i,j,k] = vy * i
        a[3,i,j,k] = vz * i
    end
    lx, ly = 10, sqrt(2) * 10
    xc, yc = 21/2, 21/2
    x1, y1 = ceil(Int, xc-lx/2), ceil(Int, yc-ly/2)
    x2, y2 = floor(Int, xc+lx/2), floor(Int, yc+ly/2)

    p1 = vector_field_projection(a, 45, "x")
    rad = pi/4
    dens = -1 * vz * cos(rad) + vy * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = ly - 2 * abs(j-yc)
        @test isapprox(p1[i,j,1], thickness*dens* i, atol=1* i)
    end

    rad = -pi/4
    p1 = vector_field_projection(a, -45, "x")
    dens = -1 * vz * cos(rad) + vy * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = ly - 2 * abs(j-yc)
        @test isapprox(p1[i,j,1], thickness*dens* i, atol=1* i)
    end

    rad = 3*pi/4
    p1 = vector_field_projection(a, 135, "x")
    dens = -1 * vz * cos(rad) + vy * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = ly - 2 * abs(j-yc)
        @test isapprox(p1[i,j,1], thickness*dens* i, atol=1* i)
    end

    rad = 5*pi/4
    p1 = vector_field_projection(a, 225, "x")
    dens = -1 * vz * cos(rad) + vy * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = ly - 2 * abs(j-yc)
        @test isapprox(p1[i,j,1], thickness*dens* i, atol=1* i)
    end
end

function test_vector_field_projectiony()
    a = zeros(3,20,20,20)
    for i = 6:15,j=6:15,k=6:15
        a[1,i,j,k] = 2 * j
        a[2,i,j,k] = 3 * j
        a[3,i,j,k] = 5 * j
    end
    lx, ly = sqrt(2) * 10, 10
    xc, yc = 21/2, 21/2
    x1, y1 = ceil(Int, xc-lx/2), ceil(Int, yc-ly/2)
    x2, y2 = floor(Int, xc+lx/2), floor(Int, yc+ly/2)

    p1 = vector_field_projection(a, 45, "y")
    rad = pi/4
    dens = -5 * cos(rad) - 2 * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = lx - 2 * abs(i-xc)
        @test isapprox(p1[i,j,1], thickness*dens * j, atol=1*j)
    end

    rad = -pi/4
    p1 = vector_field_projection(a, -45, "y")
    dens = -5 * cos(rad) - 2 * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = lx - 2 * abs(i-xc)
        @test isapprox(p1[i,j,1], thickness*dens* j, atol=1*j)
    end

    rad = 3*pi/4
    p1 = vector_field_projection(a, 135, "y")
    dens = -5 * cos(rad) - 2 * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = lx - 2 * abs(i-xc)
        @test isapprox(p1[i,j,1], thickness*dens* j, atol=1*j)
    end

    rad = 5*pi/4
    p1 = vector_field_projection(a, 225, "y")
    dens = -5 * cos(rad) - 2 * sin(rad)
    for i=x1:x2, j=y1:y2
        thickness = lx - 2 * abs(i-xc)
        @test isapprox(p1[i,j,1], thickness*dens* j, atol=1*j)
    end
end

test_warp()
test_radon()
test_radon_3d()
test_vector_field_projection()
test_vector_field_projectionx()
test_vector_field_projectiony()