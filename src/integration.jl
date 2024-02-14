using OhMyThreads: tmap, tforeach
using ProgressMeter: @showprogress
using ResizableArrays: ResizableArray

import ResizableArrays: ResizableVector

export Step, VelocityVerlet, MetropolisHastings, integrate!, simulation_time

# See https://github.com/emmt/ResizableArrays.jl/pull/3
ResizableVector(A::AbstractArray{T}) where {T} = ResizableArray{T,1}(A)

struct Step{T,S}
    dt::T
    snapshot::Vector{S}
end

abstract type Integrator end
struct VelocityVerlet <: Integrator end
struct MetropolisHastings <: Integrator
    beta::Float64
end

function integrate!(particles, cell::Cell, Δt, ::VelocityVerlet)
    # Parallel computation of initial accelerations
    accelerations = tmap(particles) do particle
        Acceleration(particle)(particles, cell)
    end
    # Parallel update of particle positions and half-step velocities
    tforeach(eachindex(particles)) do i
        particle = particles[i]
        particle.velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt / 2) = 𝐯(t) + 𝐚(t) Δt / 2
        particle.coordinates += particle.velocity * Δt  # 𝐫(t + Δt) = 𝐫(t) + 𝐯(t + Δt / 2) Δt
        movein!(particle, cell)  # Ensure particle is within cell bounds, i.e., back to `0 - L` range
    end
    # Re-compute accelerations after position updates
    accelerations = tmap(particles) do particle
        Acceleration(particle)(particles, cell)
    end
    # Parallel update of final velocities
    tforeach(eachindex(particles)) do i
        particles[i].velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt) = 𝐯(t + Δt / 2) + 𝐚(t + Δt) Δt / 2
    end
    return particles
end
function integrate!(particles, cell::Cell, Δt, n, integrator::Integrator)
    trajectory = @showprogress map(1:n) do _
        integrate!(particles, cell, Δt, integrator)
        Step(Δt, deepcopy(particles))
    end
    return ResizableVector(trajectory)
end
function integrate!(particles, cell::Cell, δv, δr, integrator::MetropolisHastings)
    for (i, particle) in enumerate(particles)
        r = rand(3) .- 0.5  # Random numbers from -0.5 to 0.5
        velocity = particle.velocity .+ δv .* r
        coordinates = particle.coordinates .+ δr .* r
        map!(Base.Fix2(mod, cell.side_length), coordinates, coordinates)  # Move `𝐫` back to `0 - L` range
        new_particle = Particle(coordinates, velocity)
        new_particles = map(enumerate(particles)) do (j, old_particle)
            j == i ? new_particle : old_particle
        end
        ΔU = potential_energy(new_particles) - potential_energy(particles)
        ΔK = kinetic_energy(new_particle) - kinetic_energy(particle)
        ΔE = ΔU + ΔK
        P = exp(-integrator.beta * ΔE)
        if P > rand()
            particles[i] = new_particle
        end
    end
    return particles
end
function integrate!(particles, cell::Cell, δv, δr, Δt, n, integrator::MetropolisHastings)
    trajectory = @showprogress map(1:n) do _
        integrate!(particles, cell, δv, δr, integrator)
        Step(Δt, deepcopy(particles))
    end
    return ResizableVector(trajectory)
end

function relax!(particles::Particles{M,C,V}, cell, Δt, n) where {M,C,V}
    integrate!(particles, cell, Δt, n, VelocityVerlet())
    init_velocities!(particles, Constant(zeros(V, 3)))
    return particles
end

simulation_time(trajectory) = cumsum(step.dt for step in trajectory)
