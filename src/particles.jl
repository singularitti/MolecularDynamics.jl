using LinearAlgebra: norm
using StaticArrays: MVector, FieldVector
using StructEquality: @struct_hash_equal_isequal_isapprox

export Position, Velocity, Acceleration, Particle
export distance,
    find_neighbors,
    boxsize,
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

struct CubicBox
    side_length::Float64
end

distance(𝐫, 𝐫′) = norm(𝐫 .- 𝐫′)
distance(a::Particle, b::Particle) = distance(a.position, b.position)

function find_nearest_image(cell::Cell, b::Particle)
    L = boxlength(cell)
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

function find_neighbors(cell::Cell, i::Integer)
    a = cell.particles[i]
    return map(filter(!=(i), eachindex(cell.particles))) do j
        b = cell.particles[j]
        find_nearest_image(cell, b)(a)
    end
end
function find_neighbors(cell::Cell, a::Particle)
    @assert a in cell
    return map(filter(!=(a), cell.particles)) do b
        find_nearest_image(cell, b)(a)
    end
end
function find_neighbors(cell::Cell, i::Integer, new_position)
    a = cell.particles[i]
    a′ = Particle(new_position, a.velocity)
    return map(filter(!=(i), eachindex(cell.particles))) do j
        b = cell.particles[j]
        find_nearest_image(cell, b)(a′)
    end
end
function find_neighbors(cell::Cell, a::Particle, new_position)
    @assert a in cell
    a′ = Particle(new_position, a.velocity)
    return map(filter(!=(a), cell.particles)) do b
        find_nearest_image(cell, b)(a′)
    end
end

function init_positions!(cell::Cell)
    L = boxlength(cell)
    for particle in eachparticle(cell)
        particle.position = L * rand(3)
    end
    @assert unique(cell.particles) == cell.particles
    return cell
end

function init_velocities!(cell::Cell)
    for particle in eachparticle(cell)
        particle.velocity = zeros(Velocity)
    end
    return cell
end

function init!(cell)
    init_positions!(cell)
    init_velocities!(cell)
    return cell
end

function damp!(cell, n, Δt)
    take_n_steps!(cell, n, Δt, VelocityVerlet())
    init_velocities!(cell)
    return cell
end

boxsize(box::CubicBox) = box.side_length^3

number_density(particles, box::CubicBox) = length(particles) / boxsize(box)

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
