using MolecularDynamics
using Test

@testset "MolecularDynamics.jl" begin
    include("particles.jl")
    include("mechanics.jl")
end
