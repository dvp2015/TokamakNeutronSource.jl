using TokamakNeutronSource
using Documenter

DocMeta.setdocmeta!(TokamakNeutronSource, :DocTestSetup, :(using TokamakNeutronSource); recursive=true)

makedocs(;
    modules=[TokamakNeutronSource],
    authors="dvp2015 <dmitri_portnov@yahoo.com>",
    repo="https://github.com/dvp2015/TokamakNeutronSource.jl/blob/{commit}{path}#{line}",
    sitename="TokamakNeutronSource.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dvp2015.github.io/TokamakNeutronSource.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dvp2015/TokamakNeutronSource.jl",
    devbranch="master",
)
