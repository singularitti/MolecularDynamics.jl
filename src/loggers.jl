using ElasticArrays: ElasticVector

export Step, Logger
export extract, nsteps, simulation_time

struct Step
    Δt::Float64
    snapshot::Vector{Particle}
end

struct Logger
    trajectory::ElasticVector{Step}
end
Logger() = Logger(ElasticVector(Step[]))

function extract(::Type{T}, logger::Logger) where {T}
    return map(logger.trajectory) do step
        map(step.snapshot) do particle
            extract(T, particle)
        end
    end
end
function extract(::Type{T}, logger::Logger, m::Integer) where {T}
    particles = logger.trajectory[m].snapshot
    return map(particles) do particle
        extract(T, particle)
    end
end
function extract(::Type{T}, logger::Logger, m::Integer, n::Integer) where {T}
    particle = logger.trajectory[m].snapshot[n]
    return extract(T, particle)
end
extract(::Type{Velocity}, particle::Particle) = getfield(particle, :velocity)
extract(::Type{Coordinates}, particle::Particle) = getfield(particle, :coordinates)
extract(::Type{Particle}, particle::Particle) = particle

simulation_time(logger::Logger) = cumsum(step.Δt for step in logger.trajectory)

nsteps(logger::Logger) = length(logger.trajectory)

function Base.show(io::IO, logger::Logger)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(logger)
        Base.show_default(IOContext(io, :limit => true), logger)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, summary(logger))
        print(io, length(logger.trajectory), " steps")
    end
end
