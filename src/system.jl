using StaticArrays: FieldVector, Size
using StructEquality: @struct_hash_equal_isequal_isapprox
using ToggleableAsserts: @toggled_assert

import StaticArrays: similar_type

export Coordinates,
    Velocity,
    Particle,
    CuboidCell,
    distance,
    generate_neighbors,
    cellsize,
    cellvolume,
    number_density,
    getcoordinates

struct Coordinates{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct Velocity{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

similar_type(::Type{<:Coordinates}, ::Type{T}, s::Size{(3,)}) where {T} = Coordinates{T}
similar_type(::Type{<:Velocity}, ::Type{T}, s::Size{(3,)}) where {T} = Velocity{T}

@struct_hash_equal_isequal_isapprox mutable struct Particle{M,C,V}
    mass::M
    coordinates::Coordinates{C}
    velocity::Velocity{V}
end
Particle(mass, coordinates::AbstractArray, velocity::AbstractArray) =
    Particle(mass, Coordinates(coordinates), Velocity(velocity))

const Particles{M,C,V} = AbstractVector{Particle{M,C,V}}

abstract type Cell end
struct CuboidCell{T} <: Cell
    dimensions::NTuple{3,T}
end
CuboidCell(number::Integer, density) = CuboidCell(cbrt(number / density))

distance(𝐫, 𝐫′) = sqrt(sum(abs2, 𝐫 .- 𝐫′))  # Much faster than `norm`
distance(a::Particle, b::Particle) = distance(a.coordinates, b.coordinates)

"""
    generate_neighbor(a::Particle, b::Particle, cell::CubicCell)

Return a particle based on `b`, as if it were in the nearest image to `a` under periodic boundary conditions (PBCs).

This implementation ensures that interactions between particles consider the minimum image
convention (MIC), effectively simulating an infinite system using a finite cell.
"""
function generate_neighbor(a::Particle, b::Particle, cell::CuboidCell)
    @toggled_assert b in cell "the particle is not in the simulation cell!"  # Ensures b's coordinates are wrapped into the primary simulation cell, addressing cases where b might have moved beyond the cell boundaries.
    𝐫′ = map(b.coordinates, b.coordinates - a.coordinates, cellsize(cell)) do rᵢ, Δrᵢ, Lᵢ  # Adjust coordinates for nearest image, ensuring MIC is followed.
        if Δrᵢ > Lᵢ / 2
            rᵢ - Lᵢ
        elseif Δrᵢ < -Lᵢ / 2
            rᵢ + Lᵢ
        else  # abs(Δrᵢ) <= Lᵢ / 2
            rᵢ  # Do not shift, already nearest
        end
    end
    return Particle(b.mass, 𝐫′, b.velocity)
end

function generate_neighbors(a::Particle, particles, cell::Cell)
    @toggled_assert a in particles
    return map(filter(!=(a), particles)) do b
        generate_neighbor(a, b, cell)
    end
end

cellsize(cell::CuboidCell) = cell.dimensions

cellvolume(cell::CuboidCell) = reduce(*, cellsize(cell))

number_density(particles, cell::CuboidCell) = length(particles) / cellvolume(cell)

function Base.in(particle::Particle, cell::CuboidCell)
    sizes = cellsize(cell)
    return all(@. zero(sizes) <= particle.coordinates <= sizes)
end

movein(particle::Particle, cell::CuboidCell) = Particle(
    particle.mass,
    map(Base.Fix2(mod, cell.dimensions), particle.coordinates),
    particle.velocity,
)

function movein!(particle::Particle, cell::CuboidCell)
    particle.coordinates = map(Base.Fix2(mod, cell.dimensions), particle.coordinates)
    return particle
end

function getcoordinates(particles)
    allcoordinates = map(particles) do particle
        particle.coordinates
    end
    return function (; x=true, y=true, z=true)
        results = ntuple(_ -> Float64[], 3)
        map(allcoordinates) do coordinates
            if x
                push!(results[1], coordinates.x)
            end
            if y
                push!(results[2], coordinates.y)
            end
            if z
                push!(results[3], coordinates.z)
            end
        end
        return results
    end
end

function Base.show(io::IO, ::MIME"text/plain", particle::Particle)
    summary(io, particle)
    println(io)
    for name in propertynames(particle)
        println(io, " ", name, ": ", getfield(particle, name))
    end
    return nothing
end
