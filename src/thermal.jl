using LinearAlgebra: dot
using Statistics: mean

export temperature, virial, pressure

function temperature(particles)
    vĀ² = map(particles) do particle
        sum(abs2, particle.velocity)
    end
    return 16 / length(particles) * sum(vĀ²)
end

ensemble_average(property::AbstractArray) = mean(property)

function virial(box, logger, indices)
    return @showprogress map(indices, logger.trajectory[indices]) do i, step
        particles = step.snapshot
        sum(eachindex(particles)) do j
            š« = extract(Coordinates, logger, i, j)
            š = force(j, particles, box)
            dot(š«, š)
        end
    end
end

function pressure(box, logger, indices)
    š«š = virial(box, logger, indices)
    N = length(logger.trajectory[1].snapshot)
    temp = mean([temperature(x.snapshot) for x in logger.trajectory[indices]])
    return 1 + ensemble_average(š«š) / 3N / temp
end
