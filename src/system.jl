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
    dimensions,
    volume,
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
    𝐫′ = map(b.coordinates, b.coordinates - a.coordinates, dimensions(cell)) do rᵢ, Δrᵢ, Lᵢ  # Adjust coordinates for nearest image, ensuring MIC is followed.
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

dimensions(cell::CuboidCell) = cell.dimensions

volume(cell::CuboidCell) = reduce(*, dimensions(cell))

number_density(particles, cell::CuboidCell) = length(particles) / volume(cell)

function Base.in(particle::Particle, cell::CuboidCell)
    dimensions = dimensions(cell)
    return all(@. zero(dimensions) <= particle.coordinates <= dimensions)
end

function movein(particle::Particle, cell::CuboidCell)
    coordinates = Coordinates(
        mod(coordinate, dimension) for
        (coordinate, dimension) in zip(particle.coordinates, dimensions(cell))
    )
    return Particle(particle.mass, coordinates, particle.velocity)
end

function movein!(particle::Particle, cell::CuboidCell)
    particle.coordinates = Coordinates(
        mod(coordinate, dimension) for
        (coordinate, dimension) in zip(particle.coordinates, dimensions(cell))
    )
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
