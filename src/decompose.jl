using HDF5
using MapImages
using MAT

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

function plot_eu(east, up; cmap="seismic_wide_y", vm=20, east_scale=1.0, title="", show_cbar=false)
    fig, axes = plt.subplots(1, 2, sharex=true, sharey=true)
    axim1 = axes[1].imshow(up, cmap=cmap, vmin=-vm, vmax=vm)
    # fig.colorbar(axim1, ax=axes[1])

    vmeast = east_scale * vm
    axim2 = axes[2].imshow(east, cmap=cmap, vmin=-vmeast, vmax=vmeast)
    show_cbar && fig.colorbar(axim2, ax=axes[2])

    axes[1].set_title("up")
    axes[2].set_title("east")
    fig.suptitle(title)
    return fig, axes
end

function demo_east_up(fn="velocities_prune_l1.h5"; dset="velos/1", full=false, vm=20, east_scale=1.0, shifta=0.0, shiftd=0.0, cmap="seismic_wide_y", show=true)
    if full
        asc_path, desc_path = ("/data1/scott/pecos/path78-bbox2/igrams_looked_18/", "/data4/scott/path85/stitched/igrams_looked_18/")
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
    if show
        fig, axes = plot_eu(east, up; cmap=cmap, vm=vm, title="$fn: $dset", east_scale=east_scale)
    end
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

function slope_mm_yr(dts, data)
    # Convert to "days since start" for line fitting
    gps_poly = fit_line(dts, data)
    slope = length(gps_poly) == 2 ? Polynomials.coeffs(gps_poly)[2] : Polynomials.coeffs(gps_poly)[1]
    return 365 * 10 * slope
end

function eastup_gps(stations, start_date=Date(2014,11,1), end_date=Date(2019,1,1))

    easts, ups = [], []
    for s in stations
        dts, east, north, up = get_gps_enu(s)
        push!(easts, slope_mm_yr(dts, east))
        push!(ups, slope_mm_yr(dts, up))
    end
    return easts, ups
end

function save_east_up_decomp(fname="velocities_max700_sigma3_noshrink_avgpix.h5", shifta=-0.8, shiftd=-0.5, vm=12)
# lats = read(ff, "lats");
    # lons = read(ff, "lons");
    east18, up18 = demo_east_up( ;full=true, dset="velos/1",  shifta=shifta, shiftd=shiftd, vm=vm)
    # TODO: use hunjoos
    latrange = (30.901666673336, 31.90000000028)
    lonrange = (-104.0, -103.001666673056)
    eastcut = east18[latrange, lonrange];
    upcut = up18[latrange, lonrange];
    gx, gy = MapImages.grid(east18.demrsc);
    matwrite("zoom_pecos_vertical_east.mat", Dict("lats"=> gy, "lons" => gx, "up" => upcut ./ 10 ./ 365 .* 1147, "east" => eastcut ./ 10 ./ 365 .* 1147));
end

station_overlap = ["TXMH", "TXFS", "TXAD", "TXS3", "NMHB"]
eeups(fname, shifta, shiftd, dset="velos/1", full=true, show=true) = extrema.(eastups(station_overlap, 
                                                                           demo_east_up(fname ;full=full, dset=dset, shifta=shiftd, shiftd=shifta, show=show)...))
