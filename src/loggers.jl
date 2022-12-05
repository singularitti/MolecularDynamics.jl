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

function extract(::Type{T}, logger::Logger) where {T}
    return map(logger.history) do step
        map(step.snapshot) do particle
            extract(T, particle)
        end
    end
end
function extract(::Type{T}, logger::Logger, m::Integer) where {T}
    particles = logger.history[m].snapshot
    return map(particles) do particle
        extract(T, particle)
    end
end
function extract(::Type{T}, logger::Logger, m::Integer, n::Integer) where {T}
    particle = logger.history[m].snapshot[n]
    return extract(T, particle)
end
extract(::Type{Velocity}, particle::Particle) = getfield(particle, :velocity)
extract(::Type{Coordinates}, particle::Particle) = getfield(particle, :coordinates)
extract(::Type{Particle}, particle::Particle) = particle

simulation_time(logger::Logger) = cumsum(step.Δt for step in logger.history)

nsteps(logger::Logger) = length(logger.trajectory)

function Base.show(io::IO, logger::Logger{N}) where {N}
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(logger)
        Base.show_default(IOContext(io, :limit => true), logger)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, summary(logger))
        print(io, " total ", length(logger.history), " steps for ", N, " particles")
    end
end
