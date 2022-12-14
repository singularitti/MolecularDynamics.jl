using ProgressMeter: @showprogress

export VelocityVerlet, MetropolisHastings
export take_one_step!, take_n_steps!

abstract type Integrator end
struct VelocityVerlet <: Integrator end
struct MetropolisHastings <: Integrator
    β::Float64
end

function take_one_step!(particles, box::Box, Δt, ::VelocityVerlet)
    for (particle, 𝐟) in zip(particles, force(particles, box))
        particle.velocity += 𝐟 * Δt / 2  # 𝐯(t + Δt / 2)
        particle.coordinates += particle.velocity * Δt  # 𝐫(t + Δt)
        map!(Base.Fix2(mod, box.side_length), particle.coordinates, particle.coordinates)  # Move `𝐫` back to `0 - L` range
    end
    for particle in particles
        𝐟 = force(particle, particles, box)  # 𝐚(t + Δt)
        particle.velocity += 𝐟 * Δt / 2  # 𝐯(t + Δt)
    end
    return particles
end
function take_one_step!(particles, box::Box, δv, δr, integrator::MetropolisHastings)
    for (i, particle) in enumerate(particles)
        r = rand(3) .- 0.5  # Random numbers from -0.5 to 0.5
        velocity = particle.velocity .+ δv .* r
        coordinates = particle.coordinates .+ δr .* r
        map!(Base.Fix2(mod, box.side_length), coordinates, coordinates)  # Move `𝐫` back to `0 - L` range
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

function take_n_steps!(particles, box::Box, n, Δt, integrator::Integrator)
    @showprogress for _ in 1:n
        take_one_step!(particles, box, Δt, integrator)
    end
    return particles
end
function take_n_steps!(logger::Logger, particles, box::Box, n, Δt, integrator::Integrator)
    if !isempty(logger.trajectory)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    @showprogress for _ in 1:n
        take_one_step!(particles, box, Δt, integrator)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    return particles
end
function take_n_steps!(
    logger::Logger, particles, box::Box, n, Δt, δv, δr, integrator::MetropolisHastings
)
    if !isempty(logger.trajectory)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    @showprogress for _ in 1:n
        take_one_step!(particles, box, δv, δr, integrator)
        push!(logger.trajectory, Step(Δt, deepcopy(particles)))
    end
    return particles
end
