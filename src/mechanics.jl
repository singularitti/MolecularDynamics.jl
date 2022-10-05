using LinearAlgebra: dot

export potential_energy, kinetic_energy, total_energy, acceleration, potential_gradient

function potential_energy(r::Number)
    r⁻⁶ = inv(r^6)
    r⁻¹² = r⁻⁶^2
    return 4 * (r⁻¹² - r⁻⁶)
end
potential_energy(𝐫ᵢⱼ) = potential_energy(norm(𝐫ᵢⱼ))
potential_energy(𝐫::Coordinates, 𝐫′::Coordinates) = potential_energy(𝐫 .- 𝐫′)
potential_energy(a::Particle, b::Particle) = potential_energy(a.coordinates, b.coordinates)
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
    Force(a::Particle)(b::Particle)

Calculate the force particle `b` exerts on particle `a` (direction: from `b` to `a`).
"""
function Force(a::Particle)
    return function (b::Particle)
        return Force(potential_gradient(b.coordinates .- a.coordinates))
    end
end

function acceleration(i::Integer, new_coordinates, particles, box::Box)
    neighbors = find_neighbors(i, particles, box)
    return sum(Force(Particle(new_coordinates, particles[i].velocity)), neighbors)
end
function acceleration(particle::Particle, new_coordinates, particles, box::Box)
    neighbors = find_neighbors(particle, particles, box)
    return sum(Force(Particle(new_coordinates, particle.velocity)), neighbors)
end
function acceleration(i::Integer, particles, box::Box)
    neighbors = find_neighbors(i, particles, box)
    return sum(Force(particles[i]), neighbors)
end
function acceleration(particle::Particle, particles, box::Box)
    neighbors = find_neighbors(particle, particles, box)
    return sum(Force(particle), neighbors)
end
function acceleration(particles, box)
    return map(particles) do particle
        acceleration(particle, particles, box)
    end
end
