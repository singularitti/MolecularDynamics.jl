using ProgressMeter: @showprogress

export VelocityVerlet, MetropolisHastings
export take_one_step!, take_n_steps!

abstract type Integrator end
struct VelocityVerlet <: Integrator end
struct MetropolisHastings <: Integrator
    β::Float64
end

function take_one_step!(
    particles::AbstractArray{<:Particle{M,C,V}}, cell::Cell, Δt, ::VelocityVerlet
) where {M,C,V}
    accelerations = Vector{Acceleration{typeof(zero(V) / Δt)}}(undef, length(particles))
    # Parallel computation of initial accelerations
    Threads.@threads for i in eachindex(particles)
        accelerations[i] = Acceleration(particles[i])(particles, cell)
    end
    # Parallel update of particle positions and half-step velocities
    Threads.@threads for i in eachindex(particles)
        particle = particles[i]
        particle.velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt / 2) = 𝐯(t) + 𝐚(t) Δt / 2
        particle.coordinates += particle.velocity * Δt  # 𝐫(t + Δt) = 𝐫(t) + 𝐯(t + Δt / 2) Δt
        movein!(particle, cell)  # Move `𝐫` back to `0 - L` range
    end
    # Re-compute accelerations after position updates
    Threads.@threads for i in eachindex(particles)
        accelerations[i] = Acceleration(particles[i])(particles, cell)  # 𝐚(t + Δt)
    end
    # Parallel update of final velocities
    Threads.@threads for i in eachindex(particles)
        particles[i].velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt) = 𝐯(t + Δt / 2) + 𝐚(t + Δt) Δt / 2
    end
    return particles
end
function take_one_step!(particles, cell::Cell, δv, δr, integrator::MetropolisHastings)
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
        P = exp(-integrator.β * ΔE)
        if P > rand()
            particles[i] = new_particle
        end
    end
    return particles
end

function take_n_steps!(particles, cell::Cell, n, Δt, integrator::Integrator)
    trajectory = @showprogress map(1:n) do _
        take_one_step!(particles, cell, Δt, integrator)
        Step(Δt, deepcopy(particles))
    end
    return trajectory
end
function take_n_steps!(particles, cell::Cell, n, Δt, δv, δr, integrator::MetropolisHastings)
    trajectory = @showprogress map(1:n) do _
        take_one_step!(particles, cell, δv, δr, integrator)
        Step(Δt, deepcopy(particles))
    end
    return trajectory
end
