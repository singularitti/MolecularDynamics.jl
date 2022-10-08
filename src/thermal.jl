using LinearAlgebra: dot

export temperature, pressure

function temperature(particles)
    v² = map(particles) do particle
        sum(abs2, particle.velocity)
    end
    return 16 / length(particles) * sum(v²)
end

ensemble_average(property::AbstractArray) = sum(property) / length(property)

function pressure(box, logger, start, stop, step=10)
    virial = @showprogress map(enumerate(logger.history[start:step:stop])) do (i, step)
        particles = step.snapshot
        map(eachindex(particles)) do j
            𝐫 = extract(Coordinates, logger, i + start - 1, j)
            𝐟 = force(i, particles, box)
            dot(𝐫, 𝐟)
        end
    end
    N = length(logger.history[1].snapshot)
    return 1 + sum(ensemble_average(virial)) / 3N
end
