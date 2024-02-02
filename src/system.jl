using StaticArrays: FieldVector, Size
using StructEquality: @struct_hash_equal_isequal_isapprox

import StaticArrays: similar_type

export Coordinates, Velocity, Force, Acceleration, Particle, CubicCell
export distance,
    find_neighbors, boxsize, boxvolume, number_density, getcoordinates, getvelocities

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

struct Force{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct Acceleration{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

similar_type(::Type{<:Coordinates}, ::Type{T}, s::Size{(3,)}) where {T} = Coordinates{T}
similar_type(::Type{<:Velocity}, ::Type{T}, s::Size{(3,)}) where {T} = Velocity{T}
similar_type(::Type{<:Force}, ::Type{T}, s::Size{(3,)}) where {T} = Force{T}
similar_type(::Type{<:Acceleration}, ::Type{T}, s::Size{(3,)}) where {T} = Acceleration{T}

@struct_hash_equal_isequal_isapprox mutable struct Particle{M,C,V}
    mass::M
    coordinates::Coordinates{C}
    velocity::Velocity{V}
end
Particle(mass, coordinates::AbstractArray, velocity::AbstractArray) =
    Particle(mass, Coordinates(coordinates), Velocity(velocity))

abstract type Cell end
struct CubicCell{T} <: Cell
    side_length::T
end
CubicCell(number::Integer, density) = CubicCell(cbrt(number / density))

distance(𝐫, 𝐫′) = sqrt(sum(abs2, 𝐫 .- 𝐫′))  # Much faster than `norm`
distance(a::Particle, b::Particle) = distance(a.coordinates, b.coordinates)

"""
    find_nearest_image(b::Particle, cell::CubicCell)

Return a function that, when given a particle `a`, computes the position of `b` as if it
were in the nearest image to `a` under periodic boundary conditions (PBCs).

This implementation ensures that interactions between particles consider the minimum image
convention (MIC), effectively simulating an infinite system using a finite cell.

# Arguments
- `b::Particle`: The particle for which we want to find the nearest image relative to another particle `a`.
- `cell::CubicCell`: The simulation cell which defines the boundaries for PBCs.
"""
function find_nearest_image(b::Particle, cell::CubicCell)
    L = cell.side_length
    𝐫 = map(Base.Fix2(mod, L), b.coordinates)  # Ensures b's coordinates are wrapped into the primary simulation cell, addressing cases where b might have moved beyond the cell boundaries.
    return function (a::Particle)
        Δ𝐫 = 𝐫 - a.coordinates  # Compute displacement
        𝐫′ = map(𝐫, Δ𝐫) do rᵢ, Δrᵢ  # Adjust coordinates for nearest image, ensuring MIC is followed.
            if Δrᵢ > L / 2
                rᵢ - L
            elseif Δrᵢ < -L / 2
                rᵢ + L
            else  # abs(Δrᵢ) <= L / 2
                rᵢ  # Do not shift, already nearest
            end
        end
        return Particle(b.mass, 𝐫′, b.velocity)
    end
end

function find_neighbors(i::Integer, particles, cell::Cell)
    return map(filter(!=(i), eachindex(particles))) do j
        find_nearest_image(particles[j], cell)(particles[i])
    end
end
function find_neighbors(a::Particle, particles, cell::Cell)
    @assert a in particles
    return map(filter(!=(a), particles)) do b
        find_nearest_image(b, cell)(a)
    end
end

boxsize(cell::CubicCell) = ntuple(_ -> cell.side_length, 3)

boxvolume(cell::CubicCell) = reduce(*, boxsize(cell))

number_density(particles, cell::CubicCell) = length(particles) / boxvolume(cell)

Base.in(particle::Particle, cell::CubicCell) = all(
    zero(cell.side_length) <= coordinate <= cell.side_length for
    coordinate in particle.coordinates
)

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

function getvelocities(particles)
    velocities = map(particles) do particle
        particle.velocity
    end
    return function (; x=true, y=true, z=true)
        results = ntuple(_ -> Float64[], 3)
        map(velocities) do velocity
            if x
                push!(results[1], velocity.x)
            end
            if y
                push!(results[2], velocity.y)
            end
            if z
                push!(results[3], velocity.z)
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
