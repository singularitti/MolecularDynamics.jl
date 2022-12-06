using ProgressMeter: @showprogress

export VelocityVerlet
export take_one_step!, take_n_steps!

abstract type Integrator end
struct VelocityVerlet <: Integrator end

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

function take_n_steps!(particles, box::Box, n, Δt, ::VelocityVerlet)
    @showprogress for _ in 1:n
        take_one_step!(particles, box, Δt, VelocityVerlet())
    end
    return particles
end
function take_n_steps!(
    logger::Logger{N}, particles, box::Box, n, Δt, ::VelocityVerlet
) where {N}
    if !isempty(logger.trajectory)
        push!(logger.trajectory, Step(Δt, Tuple(deepcopy(particles))))
    end
    @showprogress for _ in 1:n
        take_one_step!(particles, box, Δt, VelocityVerlet())
        push!(logger.trajectory, Step(Δt, Tuple(deepcopy(particles))))
    end
    return particles
end
