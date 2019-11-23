using InsarLOS
using Test


@testset "InsarLOS.jl" begin

    dbfile = joinpath(@__DIR__, "S1A_IW_SLC__1SDV_20150316T005120_20150316T005148_005050_006570_8AB8.SAFE.db.1")
    lat = 31.0
    lon = -103.0
    xyz = [0.72732057375709702; 0.60360651252993913; -0.32659449018789627]
    @test isapprox(InsarLOS.calculate_los_xyz(lat, lon, dbfile=dbfile), xyz, atol=1e-3)
    enu = [0.572897473168904; 0.107232543160571; -0.812582098574662]
    @test isapprox(InsarLOS.rotate_xyz_to_enu(xyz, lat, lon), enu, atol=1e-3)

    # Combo function
    @test isapprox(InsarLOS.get_los_enu([lat, lon], dbfile=dbfile), enu, atol=1e-3)

    lat = 31.0
    lon = -102.0
    xyz = [0.78237765815518601; 0.55250630461773964; -0.28744039969575685]
    @test isapprox(InsarLOS.calculate_los_xyz(lat, lon, dbfile=dbfile), xyz, atol=1e-3)
    enu = [0.650408309211683; 0.115737820558811; -0.750715517490279]
    @test isapprox(InsarLOS.rotate_xyz_to_enu(xyz, lat, lon), enu, atol=1e-3)



    # InsarLOS._compute_xyz_ground(31., -102.)
    # xyzground = [-1.1378e6; -5.3529e6; 3.26624e6];


#     enumat = hcat([0.75734; 0.20000; -0.62163], [0.70140; 0.19368; -0.68559])
#     llmat = hcat([31.; -103.], [31.; -102.])
#     @show enumat, llmat
#     @test isapprox(InsarLOS.get_los_enu(llmat, dbfile=dbfile), enumat, atol=1e-3)
end
