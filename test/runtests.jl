using InsarLOS
using Test


@testset "InsarLOS.jl" begin

    dbfile = joinpath(@__DIR__, "S1A_IW_SLC__1SDV_20150316T005120_20150316T005148_005050_006570_8AB8.SAFE.db.1")
    lat = 31.0
    lon = -103.0
    expected_xyz = [0.727; 0.603; -0.3265];
    @test isapprox(InsarLOS.calculate_los_xyz(lat, lon, dbfile=dbfile), expected_xyz, atol=1e-3)

    lat = 31.0
    lon = -102.0
    expected_xyz = [0.78237; 0.55250; -0.28744]
    @test isapprox(InsarLOS.calculate_los_xyz(lat, lon, dbfile=dbfile), expected_xyz, atol=1e-3)


    lat = 31.0
    lon = -103.0
    xyz = [0.727; 0.603; -0.3265];
    expected_enu = [0.70140; 0.19368; -0.68559]
    @test isapprox(InsarLOS.rotate_xyz_to_enu(xyz, lat, lon), expected_enu, atol=1e-3)

    lat = 31.0
    lon = -102.0
    xyz = [0.78237; 0.55250; -0.28744]
    expected_enu = [0.75734; 0.20000; -0.62163]
    @test isapprox(InsarLOS.rotate_xyz_to_enu(xyz, lat, lon), expected_enu, atol=1e-3)

    lat = 31.0
    lon = -103.0
    expected_enu = [0.70140; 0.19368; -0.68559]
    @test isapprox(InsarLOS.get_los_enu([lat, lon], dbfile=dbfile), expected_enu, atol=1e-3)


#     enumat = hcat([0.75734; 0.20000; -0.62163], [0.70140; 0.19368; -0.68559])
#     llmat = hcat([31.; -103.], [31.; -102.])
#     @show enumat, llmat
#     @test isapprox(InsarLOS.get_los_enu(llmat, dbfile=dbfile), enumat, atol=1e-3)
end
