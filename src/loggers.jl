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
Logger(N::Integer) = Logger{N}(ElasticVector(Step{N}[]))

function extract(::Type{Velocity}, logger::Logger)
    return map(logger.history) do step
        map(step.snapshot) do particle
            extract(Velocity, particle)
        end
    end
end
function extract(::Type{Velocity}, logger::Logger, m::Integer)
    return map(filter(==(m), eachindex(logger.history))) do step
        map(step.snapshot) do particle
            extract(Velocity, particle)
        end
    end
end
function extract(::Type{Velocity}, logger::Logger, m::Integer, n::Integer)
    return map(filter(==(m), eachindex(logger.history))) do step
        particles = step.snapshot
        map(filter(==(n), eachindex(particles))) do i
            extract(Velocity, particles[i])
        end
    end
end
extract(::Type{Velocity}, particle::Particle) = getfield(particle, :velocity)
extract(::Type{Coordinates}, particle::Particle) = getfield(particle, :coordinates)

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
