using HDF5
using MapImages

function solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
    # @show size(asc_los_map), size(desc_los_map), size(asc_img), size(desc_img)

    east = similar(asc_img)
    up = similar(asc_img)

    @show size(asc_los_map), size(asc_img), size(desc_img)
    for jj = 1:size(asc_los_map, 2)
        for ii = 1:size(asc_los_map, 1)
            asc_eu = asc_los_map[ii, jj, [1, 3]]
            desc_eu = desc_los_map[ii, jj, [1, 3]]

            A = hcat(asc_eu, desc_eu)'
            b = [asc_img[ii, jj] ; desc_img[ii, jj]]

            x = A \ b
            east[ii, jj] = x[1]
            up[ii, jj] = x[2]
        end
    end
    return east, up
end

function solve_east_up(asc_path::AbstractString, desc_path::AbstractString,
                       asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, 
                       dset="velos/1")
    asc_img = permutedims(h5read(joinpath(asc_path, asc_fname), dset))
    desc_img = permutedims(h5read(joinpath(desc_path, desc_fname), dset))

    asc_img, desc_img = _mask_asc_desc(asc_img, desc_img)

    asc_los_map = permutedims(h5read(joinpath(asc_path, "los_map.h5"), "stack"), (2, 1, 3))
    desc_los_map = permutedims(h5read(joinpath(desc_path, "los_map.h5"), "stack"), (2, 1, 3))
    return solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
end

function plot_eu(east, up; cmap="seismic_wide", vm=20, east_scale=1.0, title="")
    fig, axes = plt.subplots(1, 2, sharex=true, sharey=true)
    axim1 = axes[1].imshow(up, cmap=cmap, vmin=-vm, vmax=vm)
    fig.colorbar(axim1, ax=axes[1])

    vmeast = east_scale * vm
    axim2 = axes[2].imshow(east, cmap=cmap, vmin=-vmeast, vmax=vmeast)
    fig.colorbar(axim2, ax=axes[2])

    axes[1].set_title("up")
    axes[2].set_title("east")
    fig.suptitle(title)
    return fig, axes
end

function demo_east_up(fn="velocities_prune_l1.h5"; dset="velos/1", full=false, east_scale=1.0)
    if full
        asc_path, desc_path = ("/data1/scott/pecos/path78-bbox2/igrams_looked/", "/data4/scott/path85/stitched/igrams_looked/")
        asc_fname, desc_fname = map(x -> joinpath(x, fn), (asc_path, desc_path))
        asc_img, desc_img = find_overlap_patches(asc_fname, desc_fname, dset)

        asc_demrsc = Sario.load_dem_from_h5(asc_fname)
        desc_demrsc = Sario.load_dem_from_h5(desc_fname)
        asc_los_map = permutedims(h5read(joinpath(asc_path, "los_map.h5"), "stack"), (2, 1, 3))
        desc_los_map = permutedims(h5read(joinpath(desc_path, "los_map.h5"), "stack"), (2, 1, 3))
        asc_idxs, desc_idxs = find_overlap_idxs(asc_los_map, desc_los_map, asc_demrsc, desc_demrsc)

        east, up = solve_east_up(asc_img, desc_img, asc_los_map[asc_idxs..., :], desc_los_map[desc_idxs..., :])
    else
        asc_path, desc_path = ("/data3/scott/pecos/zoom_pecos_full_78/igrams_looked/", "/data3/scott/pecos/zoom_pecos_full_85/igrams_looked/", )
        east, up = solve_east_up(asc_path, desc_path, fn, fn, dset)
    end
    fig, axes = plot_eu(east, up; cmap="seismic_wide", vm=20, title="$fn: $dset", east_scale=east_scale)
    return east, up, fig, axes
end


function find_overlap_idxs(asc_img, desc_img, asc_demrsc, desc_demrsc)
    left, right, bottom, top = intersection_corners(asc_demrsc, desc_demrsc)
    println(left, right, bottom, top)

    # asc_patch = asc_velo[32.3:30.71, -104.1:-102.31]
    # desc_patch = desc_velo[32.3:30.71, -104.1:-102.31]

    row1, col1 = nearest_pixel(asc_demrsc, top, left)
    row2, col2 = nearest_pixel(asc_demrsc, bottom, right)
    asc_idxs = (row1:row2, col1:col2)

    # asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc_demrsc, top, left)
    row2, col2 = nearest_pixel(desc_demrsc, bottom, right)
    desc_idxs = (row1:row2, col1:col2)
    # desc_patch = desc_img[row1:row2, col1:col2]

    return asc_idxs, desc_idxs
    # return asc_patch, desc_patch
end

function find_overlaps(asc_img::MapImage, desc_img::MapImage)
    left, right, bottom, top = intersection_corners(asc_demrsc, desc_demrsc)
    println(left, right, bottom, top)

    # asc_patch = asc_velo[32.3:30.71, -104.1:-102.31]
    # desc_patch = desc_velo[32.3:30.71, -104.1:-102.31]

    row1, col1 = nearest_pixel(asc_demrsc, top, left)
    row2, col2 = nearest_pixel(asc_demrsc, bottom, right)
    asc_idxs = (row1:row2, col1:col2)

    # asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc_demrsc, top, left)
    row2, col2 = nearest_pixel(desc_demrsc, bottom, right)
    desc_idxs = (row1:row2, col1:col2)
    # desc_patch = desc_img[row1:row2, col1:col2]

    return asc_idxs, desc_idxs
    # return asc_patch, desc_patch
end

function _mask_asc_desc(a, d)
    m1 = a .== 0;
    m2 = d .== 0;
    mask = m1 .| m2
    a[mask] .= 0;
    d[mask] .= 0;
    return a, d
end

function find_overlap_patches(asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, dset="velos/1")
    asc_img = permutedims(h5read(asc_fname, dset))
    desc_img = permutedims(h5read(desc_fname, dset))

    asc_demrsc = Sario.load_dem_from_h5(asc_fname)
    desc_demrsc = Sario.load_dem_from_h5(desc_fname)
    asc_idxs, desc_idxs = find_overlap_idxs(asc_img, desc_img, asc_demrsc, desc_demrsc)

    a, d = asc_img[asc_idxs...], desc_img[desc_idxs...]
    a, d = _mask_asc_desc(a, d)
    return a, d
end

# function _check_bounds(idx_arr, bound)
#     int_idxs = Int.(round.(idx_arr))
#     bad_idxs = int_idxs .< 0 .| int_idxs .>= bound
# 
#     if any(bad_idxs)
#         # Need to check for single numbers, shape ()
#         if int_idxs.shape:
#             # Replaces locations of bad_idxs with none
#             int_idxs = findall(bad_idxs, None, int_idxs)
#         else:
#             int_idxs = None
#         end
#     end
# 
#     return int_idxs
# end

"""Find the nearest row, col to a given lat and/or lon"""
function nearest_pixel(demrsc::DemRsc, lat::AbstractFloat, lon::AbstractFloat)
    @show lon - demrsc.x_first
    @show demrsc.x_step
    col_idx = 1 + (lon - demrsc.x_first) / demrsc.x_step
    # out_row_col[2] = _check_bounds(col_idx_arr, ncols)
    row_idx = 1 + (lat - demrsc.y_first) / demrsc.y_step
    # out_row_col[1] = _check_bounds(row_idx_arr, nrows)

    return Int.(round.((row_idx, col_idx)))
end
nearest_pixel(demrsc, lats::AbstractArray{AbstractFloat}, 
              lons::AbstractArray{AbstractFloat}) = [nearest_pixel(demrsc, lat, lon)
                                                     for (lat, lon) in zip(lats, lons)]

_max_min(a, b) = max(minimum(a), minimum(b))
_least_common(a, b) = min(maximum(a), maximum(b))

# TODO: switch all dems to symbols??
_symdict(d) = Dict(Symbol(k) => v for (k, v) in d)
function intersection_corners(dem1, dem2)
    """
    Returns:
        tuple[float]: the boundaries of the intersection box of the 2 areas in order:
        (lon_left,lon_right,lat_bottom,lat_top)
    """
    dem1 = _symdict(dem1)
    dem2 = _symdict(dem2)
    corners1 = MapImages.grid_corners(;dem1...)
    corners2 = MapImages.grid_corners(;dem2...)
    lons1, lats1 = zip(corners1...)
    lons2, lats2 = zip(corners2...)
    left = _max_min(lons1, lons2)
    right = _least_common(lons1, lons2)
    bottom = _max_min(lats1, lats2)
    top = _least_common(lats1, lats2)
    return left, right, bottom, top
end

