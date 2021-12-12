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

function test_Euler()
    b = [1, 0, 0]
    a = Euler(0.0, 0.0, pi/2)
    b1 = a * b
    @test isapprox(b1[1], 1, atol=1e-6)
    @test isapprox(b1[2], 0.0, atol=1e-6)
    @test isapprox(b1[3], 0.0, atol=1e-6)

    a = Euler(0.0, pi/2, 0.0)
    b1 = a * b
    @test isapprox(b1[1], 0.0, atol=1e-6)
    @test isapprox(b1[2], 0.0, atol=1e-6)
    @test isapprox(b1[3], -1.0, atol=1e-6)

    a = Euler(0.0, -pi/2, 0.0)
    b1 = a * b
    @test isapprox(b1[1], 0.0, atol=1e-6)
    @test isapprox(b1[2], 0.0, atol=1e-6)
    @test isapprox(b1[3], 1.0, atol=1e-6)

    a = Euler(pi/2, 0.0, 0.0)
    b1 = a * b
    @test isapprox(b1[1], 0.0, atol=1e-6)
    @test isapprox(b1[2], 1.0, atol=1e-6)
    @test isapprox(b1[3], 0.0, atol=1e-6)

    a = Euler(-pi/2, 0.0, 0.0)
    b1 = a * b
    @test isapprox(b1[1], 0.0, atol=1e-6)
    @test isapprox(b1[2], -1.0, atol=1e-6)
    @test isapprox(b1[3], 0.0, atol=1e-6)

end

function test_tilt()
    a = zeros(4,4,4)
    a[2,2,2] = 1
    a[3,2,2] = 10
    a[2,3,2] = 100
    a[3,3,2] = 1000

    a[2,2,3] = 10000
    a[3,2,3] = 100000
    a[2,3,3] = 1000000
    a[3,3,3] = 10000000

    p = tilt(a, 0.0, pi/2, 0.0)

    @test isapprox(p[2,2,2], 10)
    @test isapprox(p[3,2,2], 1e5)
    @test isapprox(p[2,3,2], 1e3)  
    @test isapprox(p[3,3,2], 1e7)  

    @test isapprox(p[2,2,3], 1)
    @test isapprox(p[3,2,3], 1e4)
    @test isapprox(p[2,3,3], 1e2)  
    @test isapprox(p[3,3,3], 1e6)  

    p = tilt(a, 0.0, 0.0, pi/2)

    @test isapprox(p[2,2,2], 1e4)
    @test isapprox(p[3,2,2], 1e5)
    @test isapprox(p[2,3,2], 1)  
    @test isapprox(p[3,3,2], 1e1)  

    @test isapprox(p[2,2,3], 1e6)
    @test isapprox(p[3,2,3], 1e7)
    @test isapprox(p[2,3,3], 1e2)  
    @test isapprox(p[3,3,3], 1e3)
end

function test_radon()
    a = zeros(4,4)
    a[2,2] = 1
    a[3,2] = 10
    a[2,3] = 100
    a[3,3] = 1000

    p = radon(a, 0)
    @test isapprox(p[2], 101)
    @test isapprox(p[3], 1010)

    p = radon(a, pi/2)
    @test isapprox(p[2], 11)
    @test isapprox(p[3], 1100)

    p = radon(a, pi)
    @test isapprox(p[2], 1010)
    @test isapprox(p[3], 101)

    p = radon(a, -pi/2)
    @test isapprox(p[2], 1100)
    @test isapprox(p[3], 11)
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

    p = radon_3d_object(a, 0.0, 0.0, 0.0)
    @test isapprox(p[2,2], 10001)
    @test isapprox(p[3,2], 100010)
    @test isapprox(p[2,3], 1000100)
    @test isapprox(p[3,3], 10001000)

    p = radon_3d_object(a, 0.0, 0.0, -pi/2)
    @test isapprox(p[2,2], 101)
    @test isapprox(p[3,2], 1010)
    @test isapprox(p[2,3], 1010000)
    @test isapprox(p[3,3], 10100000)

    p = radon_3d_object(a, 0.0, 0.0, pi/2)
    @test isapprox(p[2,2], 1.01e6)
    @test isapprox(p[3,2], 1.01e7)
    @test isapprox(p[2,3], 101)
    @test isapprox(p[3,3], 1010)

    p = radon_3d_object(a, 0.0, -pi/2, 0.0)
    @test isapprox(p[2,2], 110000)
    @test isapprox(p[3,2], 11)
    @test isapprox(p[2,3], 11000000)
    @test isapprox(p[3,3], 1100)

    p = radon_3d_object(a, 0.0, pi/2, 0.0)
    @test isapprox(p[2,2], 11)
    @test isapprox(p[3,2], 110000)
    @test isapprox(p[2,3], 1100)
    @test isapprox(p[3,3], 11000000)
end

function test_radon_vecfld()
    a = zeros(3,4,4,4)
    a[1,2,2,2] = 1
    a[1,3,2,2] = 2
    a[1,2,3,2] = 3
    a[1,3,3,2] = 4

    a[1,2,2,3] = 5
    a[1,3,2,3] = 6
    a[1,2,3,3] = 7
    a[1,3,3,3] = 8
    for i = 2:3, j= 2:3, k=2:3
        a[2,i,j,k] = a[1,i,j,k] * 10
        a[3,i,j,k] = a[1,i,j,k] * 100
    end

    p1 = radon_vecfld(a, 0.0, 0.0, 0.0)
    @test isapprox(p1[2,2], 600)
    @test isapprox(p1[3,2], 800)
    @test isapprox(p1[2,3], 1000)
    @test isapprox(p1[3,3], 1200)

    p1 = radon_vecfld(a, 0.0, 0.0, -pi/2)
    @test isapprox(p1[2,2], -40)
    @test isapprox(p1[3,2], -60)
    @test isapprox(p1[2,3], -120)
    @test isapprox(p1[3,3], -140)

    p1 = radon_vecfld(a, 0.0, 0.0, pi/2)
    @test isapprox(p1[2,2,1], 120)
    @test isapprox(p1[3,2,1], 140)
    @test isapprox(p1[2,3,1], 40)
    @test isapprox(p1[3,3,1], 60)


    p1 = radon_vecfld(a, 0.0, -pi/2, 0.0)
    @test isapprox(p1[2,2,1], 11)
    @test isapprox(p1[3,2,1], 3)
    @test isapprox(p1[2,3,1], 15)
    @test isapprox(p1[3,3,1], 7)

    p1 = radon_vecfld(a, 0.0, pi/2, 0.0)
    @test isapprox(p1[2,2,1], -3)
    @test isapprox(p1[3,2,1], -11)
    @test isapprox(p1[2,3,1], -7)
    @test isapprox(p1[3,3,1], -15)
end


test_warp()
test_Euler()
test_tilt()
test_radon_3d()
test_radon_vecfld()