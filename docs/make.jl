using MolecularDynamics
using Documenter

DocMeta.setdocmeta!(MolecularDynamics, :DocTestSetup, :(using MolecularDynamics); recursive=true)

makedocs(;
    modules=[MolecularDynamics],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/singularitti/MolecularDynamics.jl/blob/{commit}{path}#{line}",
    sitename="MolecularDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/MolecularDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/MolecularDynamics.jl",
    devbranch="main",
)
