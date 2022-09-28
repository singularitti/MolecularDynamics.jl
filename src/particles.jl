using StaticArrays: MVector

export Particle, SimulationCell
export distance,
    list_interacting_particles,
    eachparticle,
    boxsize,
    boxlength,
    init_positions!,
    init_velocities!,
    init!

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

function list_interacting_particles(cell::SimulationCell, i)
    return map(filter(!=(i), eachindex(cell.particles))) do j
        particleᵢ, particleⱼ = cell.particles[[i, j]]
        Δ𝐫, L = particleⱼ.position - particleᵢ.position, boxlength(cell)
        coordinates = map(Δ𝐫) do Δr  # Finid the nearest image of particle `j`
            if Δr > L / 2
                Δr - L
            elseif Δr < -L / 2
                Δr + L
            else  # abs(x) < L / 2
                Δr
            end
        end
        Particle(coordinates)
    end
end
function list_interacting_particles(cell::SimulationCell)
    return map(Base.Fix1(list_interacting_particles, cell), eachindex(cell.particles))
end

function init_positions!(cell::SimulationCell)
    L = boxlength(cell)
    for particle in eachparticle(cell)
        particle.position = L * rand(3)
    end
    @assert unique(cell.particles) == cell.particles
    return cell
end

function init_velocities!(cell::SimulationCell)
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

boxsize(cell::SimulationCell) = particlenumber(cell) / cell.density

boxlength(cell::SimulationCell) = cbrt(boxsize(cell))

particlenumber(cell::SimulationCell) = length(cell.particles)

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
