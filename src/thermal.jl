using Statistics: mean

export temperature, virial, pressure

function temperature(particles)
    v² = map(particles) do particle
        sum(abs2, particle.velocity)
    end
    return 16 / length(particles) * sum(v²)
end

ensemble_average(property::AbstractArray) = mean(property)

function virial(box, logger, indices)
    return @showprogress map(indices, logger.trajectory[indices]) do i, step
        particles = step.snapshot
        sum(eachindex(particles)) do j
            𝐫 = extract(Coordinates, logger, i, j)
            𝐟 = Force(j, particles, box)
            muladd(𝐫.x, 𝐟.x, muladd(𝐫.y, 𝐟.y, 𝐫.z * 𝐟.z))  # 3 times faster than `dot`
        end
    end
end

function pressure(box, logger, indices)
    𝐫𝐅 = virial(box, logger, indices)
    N = length(logger.trajectory[1].snapshot)
    temp = mean([temperature(x.snapshot) for x in logger.trajectory[indices]])
    return 1 + ensemble_average(𝐫𝐅) / 3N / temp
end
