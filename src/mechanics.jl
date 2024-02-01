export potential_energy, kinetic_energy, total_energy, potential_gradient

function potential_energy(r::Number)
    r⁻⁶ = inv(r^6)
    r⁻¹² = r⁻⁶^2
    return 4 * (r⁻¹² - r⁻⁶)
end
potential_energy(a::Particle, b::Particle) = potential_energy(distance(a, b))
function potential_energy(particles)
    return sum(enumerate(particles[begin:(end - 1)])) do (i, particleᵢ)
        sum(particles[(i + 1):end]) do particleⱼ
            potential_energy(particleᵢ, particleⱼ)
        end
    end
end

function potential_gradient(𝐫)
    r = norm(𝐫)
    return 𝐫 * (inv(r)^8 / 2 - inv(r)^14)
end
potential_gradient(a::Particle, b::Particle) =
    potential_gradient(b.coordinates .- a.coordinates)

kinetic_energy(particle::Particle) = 24 * sum(abs2, particle.velocity)
kinetic_energy(particles) = sum(kinetic_energy, particles)

total_energy(particles) = kinetic_energy(particles) .+ potential_energy(particles)

"""
    Force(a::Particle)(b::Particle)

Calculate the force particle `b` exerts on particle `a` (direction: from `b` to `a`).
"""
Force(a::Particle) = (b::Particle) -> Force(potential_gradient(a, b))
function Force(i::Integer, particles, box::Box)
    neighbors = find_neighbors(i, particles, box)
    return sum(Force(particles[i]), neighbors)
end
function Force(particle::Particle, particles, box::Box)
    neighbors = find_neighbors(particle, particles, box)
    return sum(Force(particle), neighbors)
end
function Force(particles, box)
    return map(particles) do particle
        Force(particle, particles, box)
    end
end

"""
    Acceleration(a::Particle)(b::Particle)

Calculate the acceleration of particle `a` due to particle `b` (direction: from `b` to `a`).
"""
Acceleration(a::Particle) = (b::Particle) -> Force(a)(b) / a.mass
