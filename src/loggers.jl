using ElasticArrays: ElasticMatrix

export ObservableLogger
export extract_velocities, extract_coordinates

struct ObservableLogger
    Δt::Float64
    history::ElasticMatrix{Particle}
end

function extract_velocities(logger::ObservableLogger)
    return map(logger.history) do particle
        particle.velocity
    end
end

function extract_coordinates(logger::ObservableLogger)
    return map(logger.history) do particle
        particle.coordinates
    end
end

function Base.show(io::IO, logger::ObservableLogger)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(logger)
        Base.show_default(IOContext(io, :limit => true), logger)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, summary(logger))
        println(io, " time step: ", logger.Δt)
        print(
            io,
            " total ",
            size(logger.history, 2),
            " steps for ",
            size(logger.history, 1),
            " particles",
        )
    end
end
