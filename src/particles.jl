using StaticArrays: MVector

export Particle, Cell
export distance,
    list_neighbors,
    eachparticle,
    boxsize,
    boxlength,
    init_positions!,
    init_velocities!,
    init!,
    damp!

mutable struct Particle
    position::MVector{3,Float64}
    velocity::MVector{3,Float64}
    Particle(position, velocity) = new(position, velocity)
    Particle(position) = new(position)
    Particle() = new()  # Incomplete initialization
end

struct Cell
    particles::Vector{Particle}
    density::Float64
end

function distance(particle::Particle, particle′::Particle)
    return sqrt(sum(abs2, particle.position - particle′.position))
end

function list_neighbors(cell::Cell, i)
    return map(filter(!=(i), eachindex(cell.particles))) do j
        particleᵢ, particleⱼ = cell.particles[[i, j]]
        𝐫ᵢⱼ, L = particleⱼ.position - particleᵢ.position, boxlength(cell)
        position = map(𝐫ᵢⱼ) do xᵢⱼ  # Finid the nearest image of particle `j`
            if xᵢⱼ > L / 2
                xᵢⱼ - L
            elseif xᵢⱼ < -L / 2
                xᵢⱼ + L
            else  # abs(x) < L / 2
                xᵢⱼ
            end
        end
        Particle(position, particleⱼ.velocity)
    end
end
function list_neighbors(cell::Cell)
    return map(Base.Fix1(list_neighbors, cell), eachindex(cell.particles))
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
        particle.velocity = zeros(MVector{3,Float64})
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
