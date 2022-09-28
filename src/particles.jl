using StaticArrays: MVector

export Particle, SimulationCell
export distance

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

function apply_periodic_bc(x, L)
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
        coordinates = map(Base.Fix1(apply_periodic_bc, Δ𝐫), boxlength(cell))
        Particle(coordinates)
    end
end

boxsize(cell::SimulationCell) = length(cell.particles) / cell.density

boxlength(cell::SimulationCell) = cbrt(boxsize(cell))
