using ElasticArrays: ElasticVector

export Step, Trajectory, simulation_time

struct Step{T,S}
    dt::T
    snapshot::Vector{S}
end

struct Trajectory{S<:Step} <: AbstractVector{S}
    data::ElasticVector{S}
end
Trajectory(data::AbstractArray) = Trajectory(ElasticVector(data))

simulation_time(trajectory::Trajectory) = cumsum(step.dt for step in trajectory)

Base.size(trajectory::Trajectory) = size(trajectory.data)

Base.getindex(trajectory::Trajectory, i::Int) = getindex(trajectory.data, i)

Base.setindex!(trajectory::Trajectory, v, i::Int) = setindex!(trajectory.data, v, i)

Base.IndexStyle(::Type{<:Trajectory}) = IndexLinear()

Base.similar(::Trajectory{T}, dims::Dims) where {T} =
    Trajectory{T}(ElasticVector{T}(undef, dims))
