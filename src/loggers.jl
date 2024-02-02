using ElasticArrays: ElasticVector

export Step, Logger
export extract, simulation_time

struct Step{T,S}
    Δt::T
    snapshot::Vector{S}
end

struct Logger{S,T}
    trajectory::ElasticVector{Step{S,T}}
end
Logger{S,T}() where {S,T} = Logger(ElasticVector(Step{S,T}[]))

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

# Implementing iteration interface for Logger
function Base.iterate(logger::Logger, state=1)
    if state > length(logger.trajectory)
        return nothing
    else
        return (logger.trajectory[state], state + 1)
    end
end

Base.length(logger::Logger) = length(logger.trajectory)

Base.eltype(logger::Logger) = eltype(logger.trajectory)

# Implementing indexing interface for Logger
Base.getindex(logger::Logger, index) = logger.trajectory[index]

Base.setindex!(logger::Logger, step, index) = (logger.trajectory[index] = step)

Base.firstindex(logger::Logger) = 1

Base.lastindex(logger::Logger) = length(logger.trajectory)

function Base.show(io::IO, ::MIME"text/plain", logger::Logger)
    summary(io, logger)
    println(io)
    print(io, length(logger.trajectory), " steps")
    return nothing
end
