using ProgressMeter: @showprogress

export VelocityVerlet, MetropolisHastings
export take_one_step!, take_n_steps!

abstract type Integrator end
struct VelocityVerlet <: Integrator end
struct MetropolisHastings <: Integrator
    β::Float64
end

function take_one_step!(particles, cell::Cell, Δt, ::VelocityVerlet)
    accelerations = [Acceleration(particle.velocity ./ Δt) for particle in particles]
    # Parallel computation of initial accelerations
    Threads.@threads for i in eachindex(particles)
        accelerations[i] = Acceleration(particles[i], particles, cell)
    end
    # Parallel update of particle positions and half-step velocities
    Threads.@threads for i in eachindex(particles)
        particle = particles[i]
        particle.velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt / 2) = 𝐯(t) + 𝐚(t) Δt / 2
        particle.coordinates += particle.velocity * Δt  # 𝐫(t + Δt) = 𝐫(t) + 𝐯(t + Δt / 2) Δt
        particle.coordinates = map(Base.Fix2(mod, cell.side_length), particle.coordinates)   # Move `𝐫` back to `0 - L` range
    end
    # for (particle, 𝐚) in zip(particles, accelerations)
    #     particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt / 2) = 𝐯(t) + 𝐚(t) Δt / 2
    #     particle.coordinates += particle.velocity * Δt  # 𝐫(t + Δt) = 𝐫(t) + 𝐯(t + Δt / 2) Δt
    #     particle.coordinates = map(Base.Fix2(mod, cell.side_length), particle.coordinates)  # Move `𝐫` back to `0 - L` range
    # end
    # Re-compute accelerations after position updates
    Threads.@threads for i in eachindex(particles)
        accelerations[i] = Acceleration(particles[i], particles, cell)  # 𝐚(t + Δt)
    end
    # Parallel update of final velocities
    Threads.@threads for i in eachindex(particles)
        particles[i].velocity += accelerations[i] * Δt / 2  # 𝐯(t + Δt) = 𝐯(t + Δt / 2) + 𝐚(t + Δt) Δt / 2
    end
    # for particle in particles
    #     𝐚 = Acceleration(particle, particles, cell)  # 𝐚(t + Δt)
    #     particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt) = 𝐯(t + Δt / 2) + 𝐚(t + Δt) Δt / 2
    # end
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
    @showprogress for _ in 1:n
        take_one_step!(particles, cell, Δt, integrator)
    end
    return particles
end
function take_n_steps!(logger::Logger, particles, cell::Cell, n, Δt, integrator::Integrator)
    if !isempty(logger.trajectory)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    @showprogress for _ in 1:n
        take_one_step!(particles, cell, Δt, integrator)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    return particles
end
function take_n_steps!(
    logger::Logger, particles, cell::Cell, n, Δt, δv, δr, integrator::MetropolisHastings
)
    if !isempty(logger.trajectory)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    @showprogress for _ in 1:n
        take_one_step!(particles, cell, δv, δr, integrator)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    return particles
end
