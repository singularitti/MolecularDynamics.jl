using StaticArrays: MVector, FieldVector

export Position, Velocity, Acceleration, Particle, Cell
export distance,
    list_neighbors,
    eachparticle,
    boxsize,
    boxlength,
    init_positions!,
    init_velocities!,
    init!,
    damp!

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

mutable struct Particle
    position::Position
    velocity::Velocity
    Particle(position, velocity) = new(position, velocity)
    Particle(position) = new(position)
    Particle() = new()  # Incomplete initialization
end

struct Cell
    particles::Vector{Particle}
    density::Float64
end

distance(𝐫, 𝐫′) = sqrt(sum(abs2, 𝐫 .- 𝐫′))
distance(a::Particle, b::Particle) = distance(a.position, b.position)

function list_images(particle::Particle, L)
    return map(Iterators.product(Iterators.repeated((-L, 0, L), 3)...)) do shift
        Particle(particle.position .+ shift)
    end
end

function find_nearest(a::Particle, L)
    function (b::Particle)
        images = list_images(b, L)
        distances = map(images) do b′
            distance(a, b′)
        end
        index = argmin(distances)
        return images[index]
    end
end

function list_neighbors(cell::Cell, a::Particle)
    L = boxlength(cell)
    f = find_nearest(a, L)
    return map(f, filter(!=(a), cell.particles))
end
function list_neighbors(cell::Cell)
    return map(cell.particles) do particle
        list_neighbors(cell, particle)
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

boxsize(cell::Cell) = particlenumber(cell) / cell.density

boxlength(cell::Cell) = cbrt(boxsize(cell))

particlenumber(cell::Cell) = length(cell.particles)

struct EachParticle
    cell::Cell
end

eachparticle(cell::Cell) = EachParticle(cell)

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::EachParticle, state=1)
    if state > length(iter)
        return nothing
    else
        return iter.cell.particles[state], state + 1
    end
end

Base.eltype(::EachParticle) = Particle

Base.length(iter::EachParticle) = length(iter.cell.particles)

Base.IteratorSize(::Type{<:EachParticle}) = Base.HasLength()
