using StaticArrays: MVector

export Particle, SimulationCell
export distance, list_interacting_particles, eachparticle, boxsize, boxlength

mutable struct Particle
    position::MVector{3,Float64}
    velocity::MVector{3,Float64}
    Particle(position, velocity) = new(position, velocity)
    Particle(position) = new(position)
    Particle() = new()  # Incomplete initialization
end

struct SimulationCell
    particles::Vector{Particle}
    density::Float64
end

function distance(particle::Particle, particle′::Particle)
    return sqrt(sum(abs2, particle.position - particle′.position))
end

function apply_pbc(x, L)
    return if x > L / 2
        x + L
    elseif x < -L / 2
        x - L
    else  # abs(x) < L / 2
        x
    end
end

function list_interacting_particles(cell::SimulationCell, i)
    return map(filter(!=(i), eachindex(cell.particles))) do j
        particleᵢ, particleⱼ = cell.particles[[i, j]]
        Δ𝐫 = particleᵢ.position - particleⱼ.position
        coordinates = map(Base.Fix2(apply_pbc, boxlength(cell)), Δ𝐫)
        Particle(coordinates)
    end
end
function list_interacting_particles(cell::SimulationCell)
    return map(Base.Fix1(list_interacting_particles, cell), eachindex(cell.particles))
end

boxsize(cell::SimulationCell) = length(cell.particles) / cell.density

boxlength(cell::SimulationCell) = cbrt(boxsize(cell))

struct EachParticle
    cell::SimulationCell
end

eachparticle(cell::SimulationCell) = EachParticle(cell)

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
