function potential_energy(particle::Particle, particle′::Particle)
    r = distance(particle, particle′)
    η = 1 / r^6
    return 4η * (η - 1)
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

function acceleration(particle::Particle)
    return function (particle′::Particle)
        η = 1 / distance(particle, particle′)
        return (particle.position - particle′.position) * (η^14 - η^8)
    end
end

function accelerationof(cell::SimulationCell, i)
    particles = list_interacting_particles(cell, i)
    return sum(acceleration(cell.particles[i]), particles)
end
accelerationof(particles) = map(Base.Fix1(accelerationof, particles), eachindex(particles))
