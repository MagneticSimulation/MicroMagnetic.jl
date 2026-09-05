using MicroMagnetic
using MicroMagnetic: SFlat, FlatTerm, mkterm, coefof, collect_terms
using Test
using LinearAlgebra
using Random

"""
SFlat{CAP} fixed-capacity symbolic dual numbers (SFLAT_PLAN.md):

  1. arithmetic vs FlatTerm reference (same first-order semantics)
  2. set_precision(SFlat{16}) plumbing (bare `SFlat` alias, buffer eltype)
  3. matrix-free mul! in SFlat mode vs AbstractFloat mode — local terms only
  4. matrix-free mul! with demag (exercises the Float64 parameter path through
     the n==0 conditional convert inside the demag kernels)
  5. dense build_matrix guard in SFlat mode

The AbstractFloat (Epsilon/FlatTerm) mode is the unbounded reference; SFlat
must reproduce it to bitwise (local arithmetic) / ≤1e-12 (full mul!) level.
"""

@testset "SFlat" begin
    @testset "arithmetic vs FlatTerm" begin
        RNG = MersenneTwister(20260905)
        for CAP in (1, 4, 16)
            S = SFlat{CAP}
            for _ in 1:50
                c1, c2, r = randn(RNG), randn(RNG), randn(RNG)
                nterms = rand(RNG, 1:min(3, CAP))
                ids = rand(RNG, 1:8, nterms)
                vs = 0.5 .* randn(RNG, nterms)
                # same op sequence on both representations:
                #   x = (c1 + Σ ε) ; y = (c2 + Σ ε) ; x+y, x*y, r*x
                ft = FlatTerm(c1, Dict{Int,Float64}(ids[1] => vs[1]))
                sf = MicroMagnetic.mkterm(S, c1, ids[1], vs[1])
                for t in 2:nterms
                    step_ft = FlatTerm(0.0, Dict{Int,Float64}(ids[t] => vs[t]))
                    step_sf = MicroMagnetic.mkterm(S, 0.0, ids[t], vs[t])
                    ft = ft + step_ft
                    sf = sf + step_sf
                end
                ft2 = FlatTerm(c2, Dict{Int,Float64}(ids[1] => -vs[1]))
                sf2 = MicroMagnetic.mkterm(S, c2, ids[1], -vs[1])

                for (a_ft, a_sf) in ((ft, sf), (ft2, sf2), (ft + ft2, sf + sf2),
                                     (ft * ft2, sf * sf2), (r * ft, r * sf),
                                     (-ft, -sf), (ft / 2.0, sf / 2.0))
                    tf = collect_terms(a_ft)
                    ts = collect_terms(a_sf)
                    @test sort(collect(keys(tf))) == sort(collect(keys(ts)))
                    for k in keys(tf)
                        @test isapprox(tf[k], ts[k]; rtol = 1e-12, atol = 1e-14)
                    end
                    @test MicroMagnetic.cpart(a_ft) ≈ MicroMagnetic.cpart(a_sf) atol = 1e-14
                end
            end
        end
    end

    @testset "coefof / capacity" begin
        S = SFlat{4}
        x = MicroMagnetic.mkterm(S, 1.5, 3, 0.25)
        @test MicroMagnetic.coefof(x, 3) == 0.25
        @test MicroMagnetic.coefof(x, 1) == 0.0
        @test MicroMagnetic.coefof(2.0, 1) == 0.0
        # overflow beyond CAP throws loudly with guidance
        acc = zero(S)
        for id in 1:4
            acc = acc + mkterm(S, 0.0, id, 1.0)
        end
        @test acc.n == 4
        @test_throws OverflowError acc + mkterm(S, 0.0, 5, 1.0)
        # tagged → Float64 refuses, constant-only converts
        @test_throws ErrorException Float64(x)
        @test Float64(zero(S)) === 0.0
        @test MicroMagnetic.SFlat{16} <: AbstractFloat
    end

    @testset "set_precision plumbing" begin
        MicroMagnetic.set_precision(SFlat{16})
        @test MicroMagnetic.Float[] == SFlat{16}
        @test eltype(MicroMagnetic.create_zeros(4)) == SFlat{16}
        MicroMagnetic.set_precision(SFlat)          # bare alias → default CAP
        @test MicroMagnetic.Float[] == SFlat{16}
        MicroMagnetic.set_precision(Float64)
    end

    # ---- matrix-free mul! cross-mode -----------------------------------------
    function _op_pair(nx, ny, nz; demag::Bool, alpha::Float64, CAP::Int=1)
        m0fun = (x, y, z) -> (sin(0.3x + 0.1y), cos(0.2y - 0.1z), sin(0.5z + 0.4x))

        set_precision(AbstractFloat)
        meshF = FDMesh(nx = nx, ny = ny, nz = nz, dx = 5e-9, dy = 5e-9, dz = 5e-9)
        simF = Sim(meshF; driver = "LLG")
        set_Ms(simF, 8e5)
        init_m0(simF, m0fun)
        add_exch(simF, 1.3e-11)
        add_anis(simF, 5e4; axis = (0, 0, 1))
        add_zeeman(simF, (0, 0, 1000))
        demag && add_demag(simF)
        opF = build_matrix(simF; matrixfree = true, alpha = alpha)

        set_precision(SFlat{CAP})
        meshS = FDMesh(nx = nx, ny = ny, nz = nz, dx = 5e-9, dy = 5e-9, dz = 5e-9)
        simS = Sim(meshS; driver = "LLG")
        set_Ms(simS, 8e5)
        init_m0(simS, m0fun)
        add_exch(simS, 1.3e-11)
        add_anis(simS, 5e4; axis = (0, 0, 1))
        add_zeeman(simS, (0, 0, 1000))
        demag && add_demag(simS)
        opS = build_matrix(simS; matrixfree = true, alpha = alpha)

        set_precision(Float64)
        return opF, opS
    end

    @testset "matrix-free mul!: SFlat vs AbstractFloat (local)" begin
        # CAP=1 is the recommended matrix-free capacity (the mul! regime carries
        # exactly one ε tag — see SFLAT_PLAN.md §1.2); CAP=16 pins the general case.
        for CAP in (1, 16)
            opF, opS = _op_pair(4, 3, 2; demag = false, alpha = 0.01, CAP = CAP)
            RNG = MersenneTwister(7)
            N2 = size(opF, 1)
            for _ in 1:5
                x = randn(RNG, N2)
                outF = zeros(N2); outS = zeros(N2)
                mul!(outF, opF, x)
                mul!(outS, opS, x)
                @test maximum(abs.(outF .- outS)) <=
                      1e-12 * max(1.0, maximum(abs, outF))
            end
        end
    end

    @testset "matrix-free mul!: SFlat vs AbstractFloat (with demag)" begin
        opF, opS = _op_pair(4, 4, 3; demag = true, alpha = 0.01, CAP = 1)
        RNG = MersenneTwister(11)
        N2 = size(opF, 1)
        for _ in 1:5
            x = randn(RNG, N2)
            outF = zeros(N2); outS = zeros(N2)
            mul!(outF, opF, x)
            mul!(outS, opS, x)
            @test maximum(abs.(outF .- outS)) <=
                  1e-12 * max(1.0, maximum(abs, outF))
        end
        # undamped variant exercises the α=0 branch of the linearisation
        opF0, opS0 = _op_pair(3, 3, 2; demag = true, alpha = 0.0, CAP = 1)
        x = randn(MersenneTwister(13), size(opF0, 1))
        outF0 = zeros(length(x)); outS0 = zeros(length(x))
        mul!(outF0, opF0, x)
        mul!(outS0, opS0, x)
        @test maximum(abs.(outF0 .- outS0)) <= 1e-12 * max(1.0, maximum(abs, outF0))
    end

    @testset "dense guard in SFlat mode" begin
        set_precision(SFlat{16})
        mesh = FDMesh(nx = 2, ny = 2, nz = 1, dx = 5e-9, dy = 5e-9, dz = 5e-9)
        sim = Sim(mesh; driver = "LLG")
        set_Ms(sim, 8e5)
        init_m0(sim, (1.0, 0.0, 0.0))
        add_exch(sim, 1.3e-11)
        @test_throws ErrorException build_matrix(sim)
        set_precision(Float64)
    end
end
