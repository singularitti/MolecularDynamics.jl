using ElasticArrays: ElasticVector
using StaticArrays: SVector

export Step, Logger
export extract, nsteps, simulation_time

struct Step{N}
    Δt::Float64
    snapshot::SVector{N,Particle}
end

struct Logger{N}
    history::ElasticVector{Step{N}}
end

function extract_history(logger::VelocityLogger)
    return map(logger.history) do particle
        particle.velocity
    end
end
function extract_history(logger::CoordinatesLogger)
    return map(logger.history) do particle
        particle.coordinates
    end
end

simulation_time(logger::Logger) = sum(step.Δt for step in logger.history; init=0)

nsteps(logger::Logger) = size(logger.history, 2)

function Base.show(io::IO, logger::Logger)
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
