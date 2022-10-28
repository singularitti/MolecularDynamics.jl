using LinearAlgebra: dot
using Statistics: mean

export temperature, virial, pressure

function temperature(particles)
    v² = map(particles) do particle
        sum(abs2, particle.velocity)
    end
    return 16 / length(particles) * sum(v²)
end

ensemble_average(property::AbstractArray) = mean(property)

function virial(box, logger, start, stop, step=10)
    range = start:step:stop
    return @showprogress map(range, logger.history[range]) do i, step
        particles = step.snapshot
        sum(eachindex(particles)) do j
            𝐫 = extract(Coordinates, logger, i, j)
            𝐟 = force(j, particles, box)
            dot(𝐫, 𝐟)
        end
    end
end

function pressure(box, logger, start, stop, step=10)
    𝐫𝐅 = virial(box, logger, start, stop, step)
    N = length(logger.history[1].snapshot)
    temp = mean([temperature(x.snapshot) for x in logger.history[start:step:stop]])
    return 1 + ensemble_average(𝐫𝐅) / 3N / temp
end
