export temperature, pressure

function temperature(particles)
    v² = map(particles) do particle
        sum(abs2, particle.velocity)
    end
    return 16 / length(particles) * sum(v²)
end

ensemble_average(property::AbstractArray) = sum(property) / length(property)

function pressure(cell::Cell, steps)
    virial = map(steps) do step
        𝐑 = getcoordinates(cell.particles)()
        𝐀 = acceleration(cell)
        dot(𝐑, 𝐀)
    end
    return 1 + ensemble_average(virial) / 3 / length(cell.particles)
end
