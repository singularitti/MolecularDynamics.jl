export potential_energy, kinetic_energy, total_energy, accelerations

function potential_energy(r::Number)
    r⁻⁶ = inv(r^6)
    return 4 * (r⁻⁶^2 - r⁻⁶)
end
function potential_energy(𝐫ᵢⱼ)
    r = norm(𝐫ᵢⱼ)
    return potential_energy(r)
end
potential_energy(𝐫::Position, 𝐫′::Position) = potential_energy(𝐫 .- 𝐫′)
potential_energy(a::Particle, b::Particle) = potential_energy(a.position, b.position)
function potential_energy(particles::AbstractVector{Particle})
    return sum(enumerate(particles[begin:(end - 1)])) do (i, particleᵢ)
        sum(particles[(i + 1):end]) do particleⱼ
            potential_energy(particleᵢ, particleⱼ)
        end
    end
end

kinetic_energy(particle::Particle) = 24 * sum(abs2, particle.velocity)
kinetic_energy(particles) = sum(kinetic_energy, particles)

total_energy(particles) = kinetic_energy(particles) + potential_energy(particles)

function Acceleration(a::Particle)
    return function by(b::Particle)
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
