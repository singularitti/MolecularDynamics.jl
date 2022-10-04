using LinearAlgebra: norm
using StaticArrays: MVector, FieldVector
using StructEquality: @struct_hash_equal_isequal_isapprox

export Position, Velocity, Acceleration, Particle, CubicBox
export distance,
    find_neighbors,
    boxsize,
    boxvolume,
    number_density,
    init_positions!,
    init_velocities!,
    init!,
    damp!,
    getcoordinates,
    getvelocities

mutable struct Position <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Velocity <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Acceleration <: FieldVector{3,Float64}
    x::Float64
    y::Float64
    z::Float64
end

@struct_hash_equal_isequal_isapprox mutable struct Particle
    position::Position
    velocity::Velocity
end
Particle(particle::Particle, velocity::Velocity) = Particle(particle.position, velocity)
Particle(position::Position, particle::Particle) = Particle(position, particle.velocity)

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

distance(𝐫, 𝐫′) = norm(𝐫 .- 𝐫′)
distance(a::Particle, b::Particle) = distance(a.position, b.position)

function find_nearest_image(b::Particle, box::Box)
    L = box.side_length
    𝐫 = map(Base.Fix2(mod, L), b.position)
    return function (a::Particle)
        Δ𝐫 = 𝐫 - a.position
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
function find_neighbors(i::Integer, new_position, particles, box::Box)
    return map(filter(!=(i), eachindex(particles))) do j
        find_nearest_image(particles[j], box)(Particle(new_position, particles[i].velocity))
    end
end
function find_neighbors(a::Particle, new_position, particles, box::Box)
    @assert a in particles
    return map(filter(!=(a), particles)) do b
        find_nearest_image(b, box)(Particle(new_position, a.velocity))
    end
end

function init_positions!(particles, box::Box)
    # for particle in particles
    #     particle.position = boxsize(box) .* rand(3)
    # end
    for (particle, r) in zip(particles, vec(collect(Iterators.product(1:10, 1:10, 1:10))))
        particle.position = collect(r) * 1.1
    end
    @assert unique(particles) == particles
    return particles
end

function init_velocities!(particles)
    for particle in particles
        particle.velocity = zeros(Velocity)
    end
    return particles
end

function init!(particles, box::Box)
    init_positions!(particles, box)
    init_velocities!(particles)
    return particles
end

function damp!(particles, box, n, Δt)
    take_n_steps!(particles, box, n, Δt, VelocityVerlet())
    init_velocities!(particles)
    return particles
end

boxsize(box::CubicBox) = ntuple(_ -> box.side_length, 3)

boxvolume(box::CubicBox) = reduce(*, boxsize(box))

number_density(particles, box::CubicBox) = length(particles) / boxvolume(box)

function Base.in(particle::Particle, box::CubicBox)
    return all(particle.position) do x
        0 <= x <= box.side_length
    end
end

function getcoordinates(particles)
    positions = map(particles) do particle
        particle.position
    end
    return function (; x=true, y=true, z=true)
        results = ntuple(_ -> Float64[], 3)
        map(positions) do position
            if x
                push!(results[1], position.x)
            end
            if y
                push!(results[2], position.y)
            end
            if z
                push!(results[3], position.z)
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
