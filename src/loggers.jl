using ElasticArrays: ElasticMatrix

export TimeStepTracker
export extract_velocities, extract_coordinates

struct TimeStepTracker
    Δt::Float64
    steps::ElasticMatrix{Particle}
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
