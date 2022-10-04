using ElasticArrays: ElasticMatrix
using ProgressMeter: @showprogress

export VelocityVerlet, TimeStepTracker
export take_one_step!, take_n_steps!, extract_velocities, extract_coordinates

abstract type Integrator end
struct VelocityVerlet <: Integrator end

struct TimeStepTracker
    Δt::Float64
    steps::ElasticMatrix{Particle}
end

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
    data = ElasticMatrix{Particle}(undef, length(particles), n)
    @showprogress for i in 1:n
        # Must use `deepcopy`!
        take_one_step!(particles, box, Δt, VelocityVerlet())
        data[:, i] = deepcopy(particles)
    end
    return TimeStepTracker(Δt, data)
end
function take_n_steps!(tracker::TimeStepTracker, particles, box::Box, n, ::VelocityVerlet)
    @showprogress for _ in 1:n
        take_one_step!(particles, box, tracker.Δt, VelocityVerlet())
        append!(tracker.steps, deepcopy(particles))
    end
    return tracker
end

function extract_velocities(tracker::TimeStepTracker)
    return map(tracker.steps) do particle
        particle.velocity
    end
end

function extract_coordinates(tracker::TimeStepTracker)
    return map(tracker.steps) do particle
        particle.coordinates
    end
end

function Base.show(io::IO, tracker::TimeStepTracker)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(tracker)
        Base.show_default(IOContext(io, :limit => true), tracker)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, summary(tracker))
        println(io, " time step: ", tracker.Δt)
        print(
            io,
            " total ",
            size(tracker.steps, 2),
            " steps for ",
            size(tracker.steps, 1),
            " particles",
        )
    end
end
