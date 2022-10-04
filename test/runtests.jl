using MolecularDynamics
using Test

@testset "MolecularDynamics.jl" begin
    include("system.jl")
    include("mechanics.jl")
end
