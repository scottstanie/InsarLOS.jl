using HDF5
using MapImages

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

function solve_east_up(asc_path::AbstractString, desc_path::AbstractString,
                       asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, 
                       dset="velos/1")
    asc_img = permutedims(h5read(joinpath(asc_path, asc_fname), dset))
    desc_img = permutedims(h5read(joinpath(desc_path, desc_fname), dset))

    asc_img, desc_img = MapImages._mask_asc_desc(asc_img, desc_img)

    asc_los_map = permutedims(h5read(joinpath(asc_path, "los_map.h5"), "stack"), (2, 1, 3))
    desc_los_map = permutedims(h5read(joinpath(desc_path, "los_map.h5"), "stack"), (2, 1, 3))
    return solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
end

function plot_eu(east, up; cmap="seismic_wide", vm=20, east_scale=1.0, title="")
    fig, axes = plt.subplots(1, 2, sharex=true, sharey=true)
    axim1 = axes[1].imshow(up, cmap=cmap, vmin=-vm, vmax=vm)
    # fig.colorbar(axim1, ax=axes[1])

    vmeast = east_scale * vm
    axim2 = axes[2].imshow(east, cmap=cmap, vmin=-vmeast, vmax=vmeast)
    # fig.colorbar(axim2, ax=axes[2])

    axes[1].set_title("up")
    axes[2].set_title("east")
    fig.suptitle(title)
    return fig, axes
end

function demo_east_up(fn="velocities_prune_l1.h5"; dset="velos/1", full=false, vm=20, east_scale=1.0, shifta=0.0, shiftd=0.0)
    if full
        asc_path, desc_path = ("/data1/scott/pecos/path78-bbox2/igrams_looked/", "/data4/scott/path85/stitched/igrams_looked/")
        asc_fname, desc_fname = map(x -> joinpath(x, fn), (asc_path, desc_path))
        asc_img, desc_img = MapImages.find_overlaps(asc_fname, desc_fname, dset)

        asc_los_map = MapImage(joinpath(asc_path, "los_map.h5"), dset_name="stack")
        desc_los_map = MapImage(joinpath(desc_path, "los_map.h5"), dset_name="stack")
        #
        asc_idxs, desc_idxs = MapImages.find_overlap_idxs(asc_los_map, desc_los_map)
        # @show size(asc_los_map[asc_idxs..., :]), size(desc_los_map[desc_idxs..., :])

        east, up = solve_east_up(asc_img .+ shifta, desc_img .+ shiftd, asc_los_map[asc_idxs..., :], desc_los_map[desc_idxs..., :])
    else
        asc_path, desc_path = ("/data3/scott/pecos/zoom_pecos_full_78/igrams_looked/", "/data3/scott/pecos/zoom_pecos_full_85/igrams_looked/", )
        east, up = solve_east_up(asc_path, desc_path, fn, fn, dset)
    end
    fig, axes = plot_eu(east, up; cmap="seismic_wide", vm=vm, title="$fn: $dset", east_scale=east_scale)
    # return east, up, fig, axes
    return east, up
end


function eastups(stations, east, up)
    easts, ups = [], []
    for s in stations
        rc = gps.station_rowcol(s, Dict(up.demrsc))
        push!(easts, east[rc...])
        push!(ups, up[rc...])
    end
    return easts, ups
end
