using LinearAlgebra: norm
using StaticArrays: FieldVector
using StructEquality: @struct_hash_equal_isequal_isapprox

export Coordinates, Velocity, Force, Particle, CubicBox, shift
export distance,
    find_neighbors,
    boxsize,
    boxvolume,
    number_density,
    init_coordinates!,
    init_velocities!,
    init!,
    damp!,
    getcoordinates,
    getvelocities

mutable struct Coordinates <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Velocity <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Force <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

@struct_hash_equal_isequal_isapprox mutable struct Particle
    coordinates::Coordinates
    velocity::Velocity
end
Particle(particle::Particle, velocity) = Particle(particle.coordinates, velocity)
Particle(coordinates, particle::Particle) = Particle(coordinates, particle.velocity)

abstract type Box end
struct CubicBox <: Box
    side_length::Float64
    function CubicBox(side_length)
        if side_length <= zero(side_length)
            throw(ArgumentError("the box's side length must be larger than zero!"))
        else
            return new(side_length)
        end
    end
end
CubicBox(number::Integer, density::Real) = CubicBox(cbrt(number / density))

distance(𝐫, 𝐫′) = norm(𝐫 .- 𝐫′)
distance(a::Particle, b::Particle) = distance(a.coordinates, b.coordinates)

function shift(𝐫::Coordinates, box::Box)
    L = box.side_length
    return map(𝐫) do rᵢ
        if rᵢ > L / 2
            rᵢ - L
        elseif rᵢ < -L / 2
            rᵢ + L
        else  # abs(rᵢ) <= L / 2
            rᵢ  # Do not shift
        end
    end
end

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
        return Particle(𝐫′, b.velocity)
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

function Base.in(particle::Particle, box::CubicBox)
    return all(particle.coordinates) do x
        0 <= x <= box.side_length
    end
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
