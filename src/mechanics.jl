export potential_energy, kinetic_energy, total_energy, accelerationof, accelerations

function potential_energy(a::Particle, b::Particle)
    r = distance(a, b)
    r⁻⁶ = inv(r^6)
    return 4 * (r⁻⁶^2 - r⁻⁶)
end
function potential_energy(particles)
    total = 0
    for (i, particleᵢ) in enumerate(particles[begin:(end - 1)])
        for particleⱼ in particles[(i + 1):end]
            total += potential_energy(particleᵢ, particleⱼ)
        end
    end
    return total
end

kinetic_energy(particle::Particle) = 24 * sum(abs2, particle.velocity)
kinetic_energy(particles) = sum(kinetic_energy, particles)

total_energy(particles) = kinetic_energy(particles) + potential_energy(particles)

function accelerationof(particle::Particle)
    return function (particle′::Particle)
        η = 1 / distance(particle, particle′)
        return (particle.position - particle′.position) * (η^14 - η^8)
    end
end

function accelerations(cell::Cell)
    return map(eachindex(cell.particles)) do i
        neighborsᵢ = list_neighbors(cell, i)
        sum(accelerationof(cell.particles[i]), neighborsᵢ)
    end
end
