using StaticArrays: FieldVector, Size
using StructEquality: @struct_hash_equal_isequal_isapprox

import StaticArrays: similar_type

export Coordinates, Velocity, Force, Acceleration, Particle, CubicBox
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

abstract type Box end
struct CubicBox{T} <: Box
    side_length::T
end
CubicBox(number::Integer, density::Real) = CubicBox(cbrt(number / density))

distance(𝐫, 𝐫′) = sqrt(sum(abs2, 𝐫 .- 𝐫′))  # Much faster than `norm`
distance(a::Particle, b::Particle) = distance(a.coordinates, b.coordinates)

function find_nearest_image(b::Particle, box::Box)
    L = box.side_length
    𝐫 = map(Base.Fix2(mod, L), b.coordinates)
    return function (a::Particle)
        Δ𝐫 = 𝐫 - a.coordinates
        𝐫′ = map(𝐫, Δ𝐫) do rᵢ, Δrᵢ
            if Δrᵢ > L / 2
                rᵢ - L
            elseif Δrᵢ < -L / 2
                rᵢ + L
            else  # abs(Δrᵢ) <= L / 2
                rᵢ  # Do not shift
            end
        end
        return Particle(b.mass, 𝐫′, b.velocity)
    end
end

function find_neighbors(i::Integer, particles, box::Box)
    return map(filter(!=(i), eachindex(particles))) do j
        find_nearest_image(particles[j], box)(particles[i])
    end
end
function find_neighbors(a::Particle, particles, box::Box)
    @assert a in particles
    return map(filter(!=(a), particles)) do b
        find_nearest_image(b, box)(a)
    end
end

boxsize(box::CubicBox) = ntuple(_ -> box.side_length, 3)

boxvolume(box::CubicBox) = reduce(*, boxsize(box))

number_density(particles, box::CubicBox) = length(particles) / boxvolume(box)

Base.in(particle::Particle, box::CubicBox) = all(
    zero(box.side_length) <= coordinate <= box.side_length for
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
