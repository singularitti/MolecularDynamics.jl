using ElasticArrays: ElasticVector

export Step, Trajectory
export extract, simulation_time

struct Step{T,S}
    Δt::T
    snapshot::Vector{S}
end

struct Trajectory{S<:Step} <: AbstractVector{S}
    data::ElasticVector{S}
end

function extract(::Type{T}, logger::Trajectory) where {T}
    return map(logger.data) do step
        map(step.snapshot) do particle
            extract(T, particle)
        end
    end
end
function extract(::Type{T}, logger::Trajectory, m::Integer) where {T}
    particles = logger.data[m].snapshot
    return map(particles) do particle
        extract(T, particle)
    end
end
function extract(::Type{T}, logger::Trajectory, m::Integer, n::Integer) where {T}
    particle = logger.data[m].snapshot[n]
    return extract(T, particle)
end
extract(::Type{Velocity}, particle::Particle) = getfield(particle, :velocity)
extract(::Type{Coordinates}, particle::Particle) = getfield(particle, :coordinates)
extract(::Type{Particle}, particle::Particle) = particle

simulation_time(logger::Trajectory) = cumsum(step.Δt for step in logger.data)

# Implementing iteration interface for Logger
function Base.iterate(logger::Trajectory, state=1)
    if state > length(logger.data)
        return nothing
    else
        return (logger.data[state], state + 1)
    end
end

Base.length(logger::Trajectory) = length(logger.data)

Base.eltype(logger::Trajectory) = eltype(logger.data)

# Implementing indexing interface for Logger
Base.getindex(logger::Trajectory, index) = logger.data[index]

Base.setindex!(logger::Trajectory, step, index) = (logger.data[index] = step)

Base.firstindex(logger::Trajectory) = 1

Base.lastindex(logger::Trajectory) = length(logger.data)

function Base.show(io::IO, ::MIME"text/plain", logger::Trajectory)
    summary(io, logger)
    println(io)
    print(io, length(logger.data), " steps")
    return nothing
end
