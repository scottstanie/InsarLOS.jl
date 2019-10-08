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
        asc_img, desc_img = find_overlaps(asc_fname, desc_fname, dset)

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



function find_overlap_idxs(asc_img::MapImage, desc_img::MapImage)
    left, right, bottom, top = intersection_corners(asc.demrsc, desc.demrsc)
    println(left, right, bottom, top)

    row1, col1 = MapImages.nearest_pixel(asc_img, top, left)
    row2, col2 = MapImages.nearest_pixel(asc_img, bottom, right)
    asc_idxs = (row1:row2, col1:col2)

    row1, col1 = MapImages.nearest_pixel(desc_img, top, left)
    row2, col2 = MapImages.nearest_pixel(desc_img, bottom, right)
    desc_idxs = (row1:row2, col1:col2)

    return asc_idxs, desc_idxs
end

function find_overlaps(asc_img::MapImage{T, 2}, desc_img::MapImage{T, 2}) where {T}
    asc_idxs, desc_idxs = find_overlap_idxs(asx_img, desc_img)
    a, d = asc_img[asc_idxs...], desc_img[desc_idxs...]
    a, d = _mask_asc_desc(a, d)
    return a, d
end

# TODO: these should really be one... figure out
function find_overlaps(asc_img::MapImage{T, 3}, desc_img::MapImage{T, 3}) where {T}
    asc_idxs, desc_idxs = find_overlap_idxs(asx_img, desc_img)
    a, d = asc_img[asc_idxs..., :], desc_img[desc_idxs..., :]
    a, d = _mask_asc_desc(a, d)
    return a, d
end

function find_overlaps(asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, dset="velos/1")
    asc_img = MapImage(asc_fname, dset)
    desc_img = MapImage(desc_fname, dset)
    return find_overlaps(asc_img, desc_img)
end

function _mask_asc_desc(a, d)
    m1 = a .== 0;
    m2 = d .== 0;
    mask = m1 .| m2
    a[mask] .= 0;
    d[mask] .= 0;
    return a, d
end

