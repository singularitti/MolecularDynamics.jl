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
Trajectory(data::AbstractArray) = Trajectory(ElasticVector(data))

function extract(::Type{T}, trajectory::Trajectory) where {T}
    return map(trajectory.data) do step
        map(step.snapshot) do particle
            extract(T, particle)
        end
    end
end
function extract(::Type{T}, trajectory::Trajectory, m::Integer) where {T}
    particles = trajectory.data[m].snapshot
    return map(particles) do particle
        extract(T, particle)
    end
end
function extract(::Type{T}, trajectory::Trajectory, m::Integer, n::Integer) where {T}
    particle = trajectory.data[m].snapshot[n]
    return extract(T, particle)
end
extract(::Type{Velocity}, particle::Particle) = getfield(particle, :velocity)
extract(::Type{Coordinates}, particle::Particle) = getfield(particle, :coordinates)
extract(::Type{Particle}, particle::Particle) = particle

simulation_time(trajectory::Trajectory) = cumsum(step.Δt for step in trajectory)

Base.size(trajectory::Trajectory) = size(trajectory.data)

Base.getindex(trajectory::Trajectory, i::Int) = trajectory.data[i]

Base.setindex!(trajectory::Trajectory, v, i::Int) = setindex!(trajectory.data, v, i)

Base.IndexStyle(::Type{<:Trajectory}) = IndexLinear()

Base.similar(::Trajectory, ::Type{T}, dims::Dims) where {T} =
    Trajectory{T}(ElasticVector{T}(undef, dims))
