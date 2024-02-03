export Step, simulation_time

struct Step{T,S}
    dt::T
    snapshot::Vector{S}
end

simulation_time(trajectory) = cumsum(step.dt for step in trajectory)
