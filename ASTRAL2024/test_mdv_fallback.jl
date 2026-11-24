using Test
using Random

include(joinpath(@__DIR__, "recompute_mdv_fallback.jl"))

@testset "MDV fallback wind" begin
    t0 = DateTime(2024, 5, 20)
    sonde = (times=[t0 - Hour(1), t0 + Hour(1)],
             u=[8.0 12.0; NaN NaN], v=[4.0 6.0; NaN NaN])
    ship = (times=[t0 - Hour(1), t0 + Hour(1)], u=[2.0, 2.0], v=[1.0, 1.0])
    daily = (times=[t0, t0 + Minute(10)],
             ur=[4.0 9.0; 6.0 NaN], vr=[2.0 3.0; 4.0 NaN])
    context = (sonde=sonde, ship=ship, daily=Dict(Date(t0) => daily))
    ur, vr = chunk_winds(context, [t0, t0 + Minute(10)], 2)
    @test ur[1, 1] ≈ 6.0
    @test vr[1, 1] ≈ 3.0
    @test ur[1, 2] == 9.0
    @test isnan(ur[2, 2])
    @test isnan(interpolate_at(ship.times, ship.u, t0 - Hour(2)))
    missing_daily = (sonde=sonde, ship=ship, daily=Dict{Date, Any}(Date(t0) => nothing))
    ur_sonde, vr_sonde = chunk_winds(missing_daily, [t0], 2)
    @test ur_sonde[1, 1] ≈ 8.0
    @test vr_sonde[1, 1] ≈ 4.0
    @test isnan(ur_sonde[1, 2])
end

@testset "Structure-function fit" begin
    Random.seed!(2)
    n = 220
    w = hcat(cumsum(randn(n)), cumsum(randn(n)))
    angles = zeros(n)
    ur = hcat(fill(6.0, n), fill(10.0, n))
    vr = zeros(n, 2)
    fit = fit_epsilon(w, angles, angles, ur, vr)
    @test all(>(0), fit.epsi)
    for iz in 1:2
        rho = filter(isfinite, collect(skipmissing(fit.rho[:, iz])))
        distance = rho .^ (3 / 2)
        @test !isempty(distance)
        @test all(48.0 - 1e-8 .<= distance .<= 125.0 + 1e-8)
        @test fit.nbins[iz] >= 3
    end
    @test fit.nrho[1] != fit.nrho[2]
    ur[:, 2] .= NaN
    fit_without_wind = fit_epsilon(w, angles, angles, ur, vr)
    @test fit_without_wind.epsi[2] == -4.0

    ci1, ci2, _, _, it1, iz1, it2, iz2 = lidarindices(60, 1, 1)
    pitch = range(0.0, 0.2; length=60)
    result = small_displacements(ci1, ci2, it1, iz1, it2, iz2, 4.0, 0.0,
                                 pitch, zeros(60), w[1:60, 1:1], 0.035, 0.035;
                                 Lmin=48.0, Lmax=125.0)
    @test length(result[4]) == count(dt -> 48 <= 4 * dt <= 125, it2 .- it1)
end

@testset "Paired comparison" begin
    a = [1.0 NaN; 3.0 1.0; 100.0 3.0]
    b = [2.0 4.0; 4.0 2.0; NaN 4.0]
    va, vb, count_paired = paired_variance(a, b)
    @test va == [2.0, 2.0]
    @test vb == [2.0, 2.0]
    @test count_paired == [2, 2]
end
