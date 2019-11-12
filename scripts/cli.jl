using ArgParse

include("../src/InsarLOS.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--outfile", "-o"
            default = "los_map.h5"
            help = "Name output .h5 file to save solution stack"
        "--full-elevation-dem"
            help = "File with full elevation.dem used to make .geo SLCs"
        "--dbfile"
            help = "Name output .db file to use for satllite params"
        "--dbpath"
            help = "Path to directory contining some relevant .db files"
        "--demrsc"
            default = "dem.rsc"
            help = "Name of .rsc file for igrams- the output LOS file will match this DEM area"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        @show arg, val
    end

    # Assign the appropriate dataset name based on the order
    InsarLOS.create_los_map(outfile=parsed_args["outfile"],
                            demrsc=parsed_args["demrsc"],
                            dbpath=parsed_args["dbpath"],
                            dbfile=parsed_args["dbfile"],
                            full_elevation_dem=parsed_args["full-elevation-dem"])
end


