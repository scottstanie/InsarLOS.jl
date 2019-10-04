using HDF5

function solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
    # @show size(asc_los_map), size(desc_los_map), size(asc_img), size(desc_img)

    east = similar(asc_img)
    up = similar(asc_img)

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

function solve_east_up(asc_path, desc_path, dset, asc_fname, desc_fname)
    asc_img = permutedims(h5read(joinpath(asc_path, asc_fname), dset))
    desc_img = permutedims(h5read(joinpath(desc_path, desc_fname), dset))

    asc_los_map = permutedims(h5read(joinpath(asc_path, "los_map.h5"), "stack"), (2, 1, 3))
    desc_los_map = permutedims(h5read(joinpath(desc_path, "los_map.h5"), "stack"), (2, 1, 3))
    return solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
end

function plot_eu(east, up; cmap="seismic_wide", vm=20)
    fig, axes = plt.subplots(1, 2, sharex=true, sharey=true)
    axim1 = axes[1].imshow(up, cmap=cmap, vmin=-vm, vmax=vm)
    fig.colorbar(axim1, ax=axes[1])

    vmeast = 0.5 * vm
    axim2 = axes[2].imshow(east, cmap=cmap, vmin=-vmeast, vmax=vmeast)
    fig.colorbar(axim2, ax=axes[2])

    axes[1].set_title("up")
    axes[2].set_title("east")
    fig.suptitle(title)
    return fig, axes
end


function find_overlap_patches(asc_img, desc_img, asc_dem_rsc, desc_dem_rsc)
    left, right, bottom, top = intersection_corners(asc_dem_rsc, desc_dem_rsc)
    println(left, right, bottom, top)

    # asc_patch = asc_velo[32.3:30.71, -104.1:-102.31]
    # desc_patch = desc_velo[32.3:30.71, -104.1:-102.31]

    row1, col1 = nearest_pixel(asc_dem_rsc, lon=left, lat=top)
    row2, col2 = nearest_pixel(asc_dem_rsc, lon=right, lat=bottom)
    asc_patch = asc_img[row1:row2, col1:col2]

    row1, col1 = nearest_pixel(desc_dem_rsc, lon=left, lat=top)
    row2, col2 = nearest_pixel(desc_dem_rsc, lon=right, lat=bottom)
    desc_patch = desc_img[row1:row2, col1:col2]

    return asc_patch, desc_patch
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
function nearest_pixel(dem_rsc, lat::AbstractFloat, lon::AbstractFloat)
    @show lon - dem_rsc["x_first"]
    @show dem_rsc["x_step"]
    col_idx = 1 + (lon - dem_rsc["x_first"]) / dem_rsc["x_step"]
    # out_row_col[2] = _check_bounds(col_idx_arr, ncols)
    row_idx = 1 + (lat - dem_rsc["y_first"]) / dem_rsc["y_step"]
    # out_row_col[1] = _check_bounds(row_idx_arr, nrows)

    return Int.(round.((row_idx, col_idx)))
end
nearest_pixel(dem_rsc, lats::AbstractArray{AbstractFloat}, 
              lons::AbstractArray{AbstractFloat}) = [nearest_pixel(dem_rsc, lat, lon)
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
    corners1 = grid_corners(;dem1...)
    corners2 = grid_corners(;dem2...)
    lons1, lats1 = zip(corners1...)
    lons2, lats2 = zip(corners2...)
    left = _max_min(lons1, lons2)
    right = _least_common(lons1, lons2)
    bottom = _max_min(lats1, lats2)
    top = _least_common(lats1, lats2)
    return left, right, bottom, top
end


function grid_corners(; kwargs...)
    """Takes sizes and spacing from .rsc info, finds corner points in (x, y) form

    Returns:
        list[tuple[float]]: the corners of the latlon grid in order:
        (top right, top left, bottom left, bottom right)
    """
    left, right, bot, top = grid_extent(;kwargs...)
    return [(right, top), (left, top), (left, bot), (right, bot)]
end


"""
Returns:
    tuple[float]: the boundaries of the latlon grid in order:
    (lon_left,lon_right,lat_bottom,lat_top)
"""
function grid_extent(;rows=nothing,
                     cols=nothing,
                     y_step=nothing,
                     x_step=nothing,
                     y_first=nothing,
                     x_first=nothing,
                     file_length=nothing,
                     width=nothing,
                     kwargs...)
    rows = isnothing(rows) ? file_length : rows
    cols = isnothing(cols) ? width : cols
    @show cols
    return (x_first, x_first .+ x_step * (cols - 1), y_first + y_step * (rows - 1), y_first)
end

