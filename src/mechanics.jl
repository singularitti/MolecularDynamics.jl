using LinearAlgebra: dot

export potential_energy, kinetic_energy, total_energy, acceleration, potential_gradient

function potential_energy(r::Number)
    r⁻⁶ = inv(r^6)
    r⁻¹² = r⁻⁶^2
    return 4 * (r⁻¹² - r⁻⁶)
end
potential_energy(𝐫ᵢⱼ) = potential_energy(norm(𝐫ᵢⱼ))
potential_energy(𝐫::Position, 𝐫′::Position) = potential_energy(𝐫 .- 𝐫′)
potential_energy(a::Particle, b::Particle) = potential_energy(a.position, b.position)
function potential_energy(particles::AbstractVector{Particle})
    return sum(eachindex(particles)) do i
        sum(filter(!=(i), eachindex(particles))) do j
            potential_energy(particles[i], particles[j])
        end
    end / 2
end

function potential_gradient(𝐫)
    r = norm(𝐫)
    return 48𝐫 * (inv(r^8) / 2 - inv(r^14))
end

kinetic_energy(particle::Particle) = 24 * sum(abs2, particle.velocity)
kinetic_energy(particles) = sum(kinetic_energy, particles)

total_energy(particles) = kinetic_energy(particles) + potential_energy(particles)

"""
    Acceleration(a::Particle)(b::Particle)

Calculate the acceleration particle `b` induces on particle `a` (direction: from `b` to `a`).
"""
function Acceleration(a::Particle)
    return function (b::Particle)
        return Acceleration(potential_gradient(b.position .- a.position))
    end
end

function acceleration(cell::Cell, i::Integer, new_position)
    particle = cell.particles[i]
    neighbors = find_neighbors(cell, particle, new_position)
    return sum(Acceleration(Particle(new_position, particle.velocity)), neighbors)
end
function acceleration(cell::Cell, particle::Particle, new_position)
    neighbors = find_neighbors(cell, particle, new_position)
    return sum(Acceleration(Particle(new_position, particle.velocity)), neighbors)
end
function acceleration(cell::Cell, i::Integer)
    neighbors = find_neighbors(cell, i)
    return sum(Acceleration(cell.particles[i]), neighbors)
end
function acceleration(cell::Cell, particle::Particle)
    neighbors = find_neighbors(cell, particle)
    return sum(Acceleration(particle), neighbors)
end
function acceleration(cell::Cell)
    return map(cell.particles) do particle
        acceleration(cell, particle)
    end
end
