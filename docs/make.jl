using MolecularDynamics
using Documenter

DocMeta.setdocmeta!(MolecularDynamics, :DocTestSetup, :(using MolecularDynamics); recursive=true)

makedocs(;
    modules=[MolecularDynamics],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="MolecularDynamics.jl",
    format=Documenter.HTML(;
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
