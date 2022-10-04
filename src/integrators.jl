using ProgressMeter: @showprogress

export VelocityVerlet
export take_one_step!, take_n_steps!

abstract type Integrator end
struct VelocityVerlet <: Integrator end

function take_one_step!(particles, box::Box, Δt, ::VelocityVerlet)
    new_coordinates = map(particles, acceleration(particles, box)) do particle, 𝐚
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt / 2)
        coordinates = particle.coordinates + particle.velocity * Δt  # 𝐫(t + Δt)
        coordinates = map(Base.Fix2(mod, box.side_length), coordinates)  # Move `𝐫` back to `0 - L` range
    end
    for (particle, coordinates) in zip(particles, new_coordinates)
        𝐚 = acceleration(particle, coordinates, particles, box)  # 𝐚(t + Δt)
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt)
    end
    for (particle, coordinates) in zip(particles, new_coordinates)
        particle.coordinates = coordinates
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
    @showprogress for _ in 1:n
        take_one_step!(particles, box, Δt, VelocityVerlet())
        push!(logger.history, Step(Δt, SVector{N}(deepcopy(particles))))
    end
    return particles
end
