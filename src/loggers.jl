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

Base.size(trajectory::Trajectory) = size(trajectory.data)

Base.getindex(trajectory::Trajectory, i::Int) = trajectory.data[i]

Base.setindex!(trajectory::Trajectory, v, i::Int) = setindex!(trajectory.data, v, i)

Base.IndexStyle(::Type{<:Trajectory}) = IndexLinear()
