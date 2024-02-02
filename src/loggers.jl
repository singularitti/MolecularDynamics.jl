using ElasticArrays: ElasticVector

export Step, Logger
export extract, nsteps, simulation_time

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

nsteps(logger::Logger) = length(logger.trajectory)

function Base.show(io::IO, ::MIME"text/plain", logger::Logger)
    summary(io, logger)
    println(io)
    print(io, length(logger.trajectory), " steps")
    return nothing
end
