using Documenter, InsarLOS

makedocs(;
    modules=[InsarLOS],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/scottstanie/InsarLOS.jl/blob/{commit}{path}#L{line}",
    sitename="InsarLOS.jl",
    authors="scott <scott.stanie@gmail.com>",
    assets=String[],
)

deploydocs(
    repo = "github.com/scottstanie/InsarLOS.jl.git",
)
