export potential_energy, kinetic_energy, total_energy, accelerations

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

function Acceleration(a::Particle)
    function by(b::Particle)
        r = distance(a, b)
        𝐚 = (a.position .- b.position) * (inv(r^14) - inv(r^8) / 2)
        return Acceleration(𝐚...)
    end
end

function accelerations(cell::Cell, particle::Particle, position)
    neighbors = list_neighbors(cell, particle)
    return sum(Acceleration(Particle(position)), neighbors)
end
function accelerations(cell::Cell, particle::Particle)
    neighbors = list_neighbors(cell, particle)
    return sum(Acceleration(particle), neighbors)
end
function accelerations(cell::Cell)
    return map(cell.particles) do particle
        accelerations(cell, particle)
    end
end
