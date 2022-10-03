using ElasticArrays: ElasticMatrix
using ProgressMeter: @showprogress

export VelocityVerlet, TimeStepTracker
export take_one_step!, take_n_steps!, velocities, positions

abstract type Integrator end
struct VelocityVerlet <: Integrator end

struct TimeStepTracker
    Δt::Float64
    steps::ElasticMatrix{Particle}
end

function take_one_step!(cell::Cell, Δt, ::VelocityVerlet)
    L = boxlength(cell)
    positions = map(eachparticle(cell), acceleration(cell)) do particle, 𝐚
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt / 2)
        position = particle.position + particle.velocity * Δt  # 𝐫(t + Δt)
        position = map(Base.Fix2(mod, L), position)  # Move `𝐫` back to `0 - L` range
    end
    for (particle, position) in zip(eachparticle(cell), positions)
        𝐚 = acceleration(cell, particle, position)  # 𝐚(t + Δt)
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt)
    end
    for (particle, position) in zip(eachparticle(cell), positions)
        particle.position = position
    end
    return cell
end

function take_n_steps!(cell::Cell, n, Δt, ::VelocityVerlet)
    data = ElasticMatrix{Particle}(undef, particlenumber(cell), n)
    @showprogress for i in 1:n
        # Must use `deepcopy`!
        take_one_step!(cell, Δt, VelocityVerlet())
        data[:, i] = deepcopy(cell.particles)
    end
    return TimeStepTracker(Δt, data)
end
function take_n_steps!(tracker::TimeStepTracker, cell::Cell, n, ::VelocityVerlet)
    @showprogress for _ in 1:n
        take_one_step!(cell, tracker.Δt, VelocityVerlet())
        append!(tracker.steps, deepcopy(cell.particles))
    end
    return tracker
end

function velocities(tracker::TimeStepTracker)
    return map(tracker.steps) do particle
        particle.velocity
    end
end

function positions(tracker::TimeStepTracker)
    return map(tracker.steps) do particle
        particle.position
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
