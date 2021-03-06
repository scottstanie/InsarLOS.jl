__precompile__(true)

module InsarLOS


include("./projections.jl")

using Sario
import MapImages

import Glob
using SQLite
using LinearAlgebra
using HDF5

const EARTH_SMA = 6378137.0
const EARTH_E2 = 0.0066943799901499996

const SOL = 299792458

const MAP_FILENAME = "los_map.h5"


"""
    get_los_enu(lat_lons::AbstractArray{<:AbstractFloat, 2}; xyz_los_vecs=nothing, dbfile=nothing, param_dict=nothing)

Calculate the LOS vector in ENU from a lat/lon Array
Can also take in the pre-computed xyz los vectors

`lat_lons` is `2 x N` for N different points
"""
function get_los_enu(lat_lons::AbstractArray{<:AbstractFloat, 2}; xyz_los_vecs=nothing, 
                     dbfile=nothing, param_dict=nothing)
    if isnothing(xyz_los_vecs)
        xyz_los_vecs = calculate_los_xyz(lat_lons, dbfile=dbfile, param_dict=param_dict)
    end
    # In projections file:
    return convert_xyz_latlon_to_enu(reshape(lat_lons, 2, :),
                                     reshape(xyz_los_vecs, 3, :))
end

function get_los_enu(lat_lons::AbstractArray{<:AbstractFloat, 1}; xyz_los_vecs=nothing,
                     dbfile=nothing, param_dict=nothing)
    return get_los_enu(reshape(lat_lons, :, 1), xyz_los_vecs=xyz_los_vecs, dbfile=dbfile, param_dict=param_dict)
end


filepath(full_path) = joinpath(splitpath(full_path)[1:end-1]...)
#
# Start of fortran converted functions:

"""Calculate the line of sight vector (from satellite to ground) in X,Y,Z coordinates at one lat/lon point"""
function calculate_los_xyz(lat::T, lon::T; dbfile::Union{String, Nothing}=nothing, param_dict=nothing) where {T<:AbstractFloat}
    if isnothing(param_dict)
        if isnothing(dbfile)
            dbfile_list = Glob.glob("*.db*")
            dbfile = dbfile_list[1]
        end
        param_dict = load_all_params(dbfile)
    end
    xyz_ground = _compute_xyz_ground(lat, lon)

    orbinfo_filename = param_dict["orbinfo"]  # The .db file doesn't save path
    dbpath = filepath(dbfile)
    orbinfo_file = joinpath(dbpath, orbinfo_filename)

    timeorbit, xx, vv = read_orbit_vector(orbinfo_file)

    for idx=1:param_dict["azimuthBursts"]
        idx = 9
        vecr = compute_burst_vec(xyz_ground, param_dict, idx, timeorbit, xx, vv)
        if !isnothing(vecr)
            return vecr
        end
    end
    return nothing
end


function calculate_los_xyz(lat::T, lon::T, dem, demrsc, param_dict, timeorbit, xx, vv) where {T<:AbstractFloat}
    xyz_ground = _compute_xyz_ground(lat, lon, dem, demrsc)
    for idx=1:param_dict["azimuthBursts"]
        idx = 9
        vecr = compute_burst_vec(xyz_ground, param_dict, idx, timeorbit, xx, vv)
        if !isnothing(vecr)
            return vecr
        end
    end
    return nothing
end

function calculate_los_xyz(lat_lon::AbstractArray{<:AbstractFloat, 1}; dbfile=nothing, param_dict=nothing)
    lat, lon = lat_lon
    return calculate_los_xyz(lat, lon, dbfile=dbfile, param_dict=param_dict)
end

function calculate_los_xyz(lat_lon_vecs::AbstractArray{<:AbstractFloat, 2}; dbfile=nothing, param_dict=nothing)
    num_vecs = size(lat_lon_vecs, 2)
    xyz_vecs = Array{eltype(lat_lon_vecs)}(undef, (3, ))
    for i=1:num_vecs
        lat, lon = lat_lon_vecs[:, i]
        xyz_vecs[:, i] .= calculate_los_xyz(lat, lon, dbfile=dbfile, param_dict=param_dict)
    end
    return xyz_vecs
end

function _compute_xyz_ground(lat, lon)
    # TODO: do i wanna get rid of this "params" file?
    # dem_file, dem_rsc_file = readlines("params")

    # TODO: what directory should this be??
    dem_rsc_file = Sario.find_rsc_file(directory="../")
    demrsc = Sario.load(dem_rsc_file)

    row, col = MapImages.latlon_to_rowcol(demrsc, lat, lon)
    # println("($lat, $lon) is at ($row, $col)")

    dem_file = replace(dem_rsc_file, ".rsc" => "")
    # Load just the 1 value from the DEM
    dem_height = Sario.load(dem_file, (row, col))

    llh = [ deg2rad(lat), deg2rad(lon), dem_height ]
    xyz_ground = llh_to_xyz(llh)
    return xyz_ground
end

function _compute_xyz_ground(lat, lon, dem, demrsc)
    row, col = MapImages.nearest_pixel(demrsc, lat, lon)

    # Load just the 1 value from the DEM
    dem_height = dem[row, col]

    llh = [ deg2rad(lat), deg2rad(lon), dem_height ]
    xyz_ground = llh_to_xyz(llh)
    return xyz_ground
end


function compute_burst_vec(xyz_ground, param_dict, idx, timeorbit, xx, vv)
    dtaz = 1.0 / param_dict["prf"]  # Nazlooks / prf
    tstart = param_dict["azimuthTimeSeconds$idx"]
    tend  = tstart + (param_dict["linesPerBurst"] - 1) * dtaz
    tmid = (tstart + tend) / 2.0

    # println("Burst $idx, Start, stop Acquisition time: $tstart,$tend")
  
    rngstart = param_dict["slantRangeTime"] * SOL / 2.0
    dmrg = SOL / 2.0 / param_dict["rangeSamplingRate"] # Nnum rnglooks * drho
    rngend = rngstart + (param_dict["samplesPerBurst"]-1)*dmrg
    rngmid = 0.50*(rngstart+rngend)

    # TODO: Does it matter if the lat/lon is within this specific burst?
    # If we interpolate to the middle it seems to always be the same
    # 
    # latlons = bounds(tstart, tend, rngstart, rngend, timeorbit, xx, vv)
    # if lat > latlons[1] || lat < latlons[2]
    #     print("$lat not within bounds $latlons")
    #     return nothing
    # end


    xyz_mid, vel_mid = intp_orbit(timeorbit, xx, vv, tmid)

    # println("Satellite midpoint time,position,velocity: $tmid $xyz_mid $vel_mid")

    tline, range = orbitrangetime(xyz_ground, timeorbit, xx, vv, tmid, xyz_mid, vel_mid)
    if isnothing(tline) || isnothing(range)
        println("Failed on burst $idx")
        return nothing
    end


    satx, satv = intp_orbit(timeorbit, xx, vv, tline)
    # TESTING: see how much the East and North change at beginning/end of orbit
    # satx = xx[:, end]
    # satv = vv[:, end]

    # @show satx, satv
    # Note: pointing AWAY from satellite (hence -satx)
    dr = xyz_ground - satx
    vecr = dr / range

    return vecr

end

function load_table_params(dbfile, param_list, param_types=nothing; tablename="file")
    db = SQLite.DB(dbfile)

    # TODO: How do you really do string interpolation for this library :\
    param_list_str = join(param_list, "\",\"")
    querystr = """SELECT * FROM $tablename WHERE name IN ("$param_list_str") """

    table = SQLite.Query(db, querystr)
    if isnothing(param_types)
        return Dict(row.name => row.value for row in table)
    else
        return Dict(row.name => parse(param_types[row.name], row.value)
                    for row in table)
    end
end

function load_all_params(dbfile; tablename="file")
    param_dict = Dict{String, Any}()
    param_list = [
        "azimuthBursts",
        "linesPerBurst",
        "samplesPerBurst",
        "prf",
        "slantRangeTime",
        "rangeSamplingRate",
    ]
    # To parse into the correct types, since all are strings
    param_types = Dict(
        "azimuthBursts" => Int,
        "linesPerBurst" => Int,
        "samplesPerBurst" => Int,
        "prf" => Float64,
        "rangeSamplingRate" => Float64,
        "slantRangeTime" => Float64,
    )
    d = load_table_params(dbfile, param_list, param_types, tablename=tablename)
    merge!(param_dict, d)

    # For this, one of each of these per burst up to "azimuthBursts"
    for idx = 1:param_dict["azimuthBursts"]
        burst_params = [
            "azimuthTimeSeconds$idx",
        ]
        types = Dict(
            "azimuthTimeSeconds$idx" => Float64,
        )
        d = load_table_params(dbfile, burst_params, types, tablename=tablename)
        merge!(param_dict, d)
    end

    # Finally, the orbinfo file is a string
    merge!(param_dict, 
           load_table_params(dbfile, ["orbinfo"], tablename=tablename))
    return param_dict
end


function read_orbit_vector(orbtiming_file)
    lines = readlines(orbtiming_file)
    # timefirst, timeend = map(line -> parse(Float64, line), lines[1:2])
    nlines, num_state_vec = map(line -> parse(Int, line), lines[3:4])

    timeorbit = zeros(num_state_vec)
    xx = Array{Float64, 2}(undef, (3, num_state_vec))
    vv = similar(xx)
    # Note: we don't use the accel
    for i = 1:num_state_vec

        floats = map(num -> parse(Float64, num), split(lines[4 + i]))
        timeorbit[i] = floats[1]
        xx[:, i] = floats[2:4]
        vv[:, i] = floats[5:7]
    end
    return timeorbit, xx, vv
end


function llh_to_xyz(llh::AbstractArray{<:AbstractFloat, 1})
    lat, lon, h = llh

    rad_earth = EARTH_SMA/sqrt(1.0 - EARTH_E2*sin(lat)^2)

    xyz = similar(llh)
    xyz[1] = (rad_earth + h)*cos(lat)*cos(lon)
    xyz[2] = (rad_earth + h)*cos(lat)*sin(lon)
    xyz[3] = (rad_earth*(1.0 - EARTH_E2) + h) * sin(lat)
    return xyz
end

function orbitrangetime(xyz, timeorbit, xx, vv, tline0, satx0, satv0, tol=1e-6)
    # starting state
    tline = tline0
    satx = satx0
    satv = satv0
    
    idx, max_iter = 1, 100
    tprev = tline + 1  # Need starting guess
    while abs(tline - tprev) > tol && idx < max_iter
        tprev = tline

        dr = xyz - satx
        range = norm(dr, 2)

        fn = dr' * satv
        fnprime = -norm(satv, 2)^2

        tline = tline - fn / fnprime

        satx, satv = intp_orbit(timeorbit, xx, vv, tline)
        if isnothing(satx) || isnothing(satv)
            return nothing, nothing
        end
        idx += 1
    end
    if abs(tline - tprev) > tol
        println("Warning: orbitrangetime didn't converge within $tol (residual $(abs(tline - tprev))")
    end

    dr = xyz - satx
    range = norm(dr, 2)

    return tline, range
end


function intp_orbit(timeorbit, xx, vv, time)
    if isinf(time) || isnan(time)
        return nothing, nothing
    end
    ilocation = (time - timeorbit[1]) / (timeorbit[2]-timeorbit[1])
    # @show time
    # @show ilocation
    ilocation = round(Int, clamp(ilocation, 2, length(timeorbit) - 2))

    xyz_mid, vel_mid = orbithermite(xx[:,ilocation-1:end], vv[:,ilocation-1:end],
                                    timeorbit[ilocation-1:end], time)
    # xyz_mid, vel_mid = orbithermite(xx, vv,timeorbit, time)

    return xyz_mid, vel_mid
end


"""orbithermite - hermite polynomial interpolation of orbits"""
function orbithermite(x,v,t,time)
# inputs
#  x - 3x4 matrix of positions at four times
#  v - 3x4 matrix of velocities
#  t - 4-vector of times for each of the above data points
#  time - time to interpolate orbit to
# 
# outputs
#  xx - position at time time
#  vv - velocity at time time
    n = 4

    h, hdot = zeros(n), zeros(n)
    f0, f1 = zeros(n), zeros(n)
    g0, g1 = zeros(n), zeros(n)
    sum, product = 0, 0


    for i=1:n
        f1[i]=time-t[i]
        sum=0.
        for j=1:n
            if j != i
                sum += 1.0/(t[i]-t[j])
            end
        end
        f0[i]=1.0 - 2.0 * (time-t[i])*sum
    end

    for i=1:n
        product=1.0
        for k=1:n
            if k != i
                product *= (time-t[k])/(t[i]-t[k])
            end
        end
        h[i]=product

        sum=0.0
        for j=1:n
            product=1.0
            for k=1:n
               if k != i && k != j
                   product *= (time-t[k])/(t[i]-t[k])
               end
            end
            if j != i
                sum += 1.0 / (t[i]-t[j])*product
            end
        end
        hdot[i]=sum
    end

    for i=1:n
        g1[i]=h[i] + 2.0 * (time-t[i]) * hdot[i]
        sum=0.
        for j=1:n
            if i != j
                sum += 1.0 / (t[i]-t[j])
            end
        end
        g0[i]=2.0 *(f0[i]*hdot[i]-h[i]*sum)
    end
    
    xout = zeros(3)
    vout = zeros(3)
    for k=1:3
        sum=0.0
        for i=1:n
            sum += (x[k,i]*f0[i]+v[k,i]*f1[i])*h[i]*h[i]
        end
        xout[k]=sum

        sum=0.
        for i=1:n
            sum += (x[k,i]*g0[i]+v[k,i]*g1[i])*h[i]  #*h[i] extra in pdf
        end
        vout[k] = sum
    end

    return xout, vout
end


function _find_db_path(geo_path)
    extra_path = joinpath(geo_path, "extra_files")
    if isdir(extra_path) && length(Glob.glob(extra_path * "/*.db*")) > 0
        return extra_path
    else
        return geo_path
    end
end



"""Create a map with 3 layers (E, N, U) of the ENU vectors for every pixel within a demrsc"""
function create_los_map(;directory=".", demrsc=nothing, outfile="los_map.h5", 
                        full_elevation_dem=nothing, dbpath=nothing, dbfile=nothing)
    if isnothing(dbfile)
        isnothing(dbpath) && error("need dbfile or dbpath")
        @show dbpath
        dbfile = Glob.glob("*.db*", dbpath)[1]
    end
    @show dbfile

    if isnothing(demrsc)
        demrsc = Sario.load(Sario.find_rsc_file(directory))
    end

    param_dict = load_all_params(dbfile)
    orbinfo_filename = param_dict["orbinfo"]  # The .db file doesn't save path
    dbpath = filepath(dbfile)
    orbinfo_file = joinpath(dbpath, orbinfo_filename)
    timeorbit, xorbit, vorbit = read_orbit_vector(orbinfo_file)


    if isnothing(full_elevation_dem)
        # TODO: this is brittle
        dem_file = replace(Sario.find_rsc_file(directory=".."), ".rsc" => "")
        @show dem_file
        full_elevation_dem = Sario.load(dem_file)
    end

    xx, yy = MapImages.grid(demrsc, sparse=true)
    enu_out = Array{Float64, 3}(undef, (length(yy), length(xx), 3))
    # lats = Array{Float64, 2}(undef, (length(yy), length(xx)))
    # lons = similar(lats)

    Threads.@threads for j in 1:length(xx)
        for i in 1:length(yy)
    # for (j, x) in enumerate(xx)  # Doesn't seem to work with thread?
        # for (i, y) in enumerate(yy)
            # @show i, j, y, x
            y = yy[i]
            x = xx[j]
            xyz_los_vecs = calculate_los_xyz(y, x, full_elevation_dem, demrsc, param_dict, timeorbit, xorbit, vorbit)
            # println("$y $x is at ", MapImages.latlon_to_rowcol(demrsc, y, x))
            enu_out[i, j, :] = get_los_enu([y, x], xyz_los_vecs=xyz_los_vecs)
        end
    end

    if !isnothing(outfile)
        println("Writing to $outfile dset 'stack'")
        h5open(outfile, "w") do f
            write(f, "stack", permutedims(enu_out, (2, 1, 3)))
            write(f, "lats", collect(yy))
            write(f, "lons", collect(xx))
        end
        Sario.save_dem_to_h5(outfile, demrsc)
    end

    return enu_out
end


_check_oob(arr::AbstractArray, num) = @assert (num > minimum(arr) && num < maximum(arr)) "$num out of bounds for $(extrema(arr))"

findnearest(A::AbstractArray, t) = findmin(abs.(A .- t))[2]


function read_los_map(lat_lons::AbstractArray{<:AbstractFloat, 2}, los_map_file::String=MAP_FILENAME)
    los_map = permutedims(h5read(los_map_file, "stack"), (2, 1, 3))
    gridlats = h5read(los_map_file, "lats")
    gridlons = h5read(los_map_file, "lons")
    return read_los_map(lat_lons, los_map, gridlats, gridlons)
end

# """Given an array of lat/lon points, read the E,N,U coeffs from the LOS map"""
function read_los_map(lat_lons::AbstractArray{<:AbstractFloat, 2}, los_map::AbstractArray{<:AbstractFloat, 3}, gridlats, gridlons)
    num_points = size(lat_lons, 2)
    out = Array{eltype(lat_lons), 2}(undef, (3, num_points))
    for i = 1:num_points
        lat, lon = lat_lons[:, i]
        _check_oob(gridlats, lat)
        row = findnearest(gridlats, lat)

        _check_oob(gridlons, lon)
        col = findnearest(gridlons, lon)

        out[:, i] = los_map[row, col, :]
    end
    return out
end

read_los_map(lat_lons::AbstractArray{<:AbstractFloat, 1}, los_map_file::String=MAP_FILENAME) = read_los_map(reshape(lat_lons, :, 1), los_map_file)

end # module
