using LinearAlgebra: dot

export potential_energy,
    kinetic_energy,
    total_energy,
    accelerations,
    potential_directional_derivative,
    potential_gradient

function potential_energy(r::Number)
    r⁻⁶ = inv(r^6)
    r⁻¹² = r⁻⁶^2
    return 4 * (r⁻¹² - r⁻⁶)
end
potential_energy(𝐫ᵢⱼ) = potential_energy(norm(𝐫ᵢⱼ))
potential_energy(𝐫::Position, 𝐫′::Position) = potential_energy(𝐫 .- 𝐫′)
potential_energy(a::Particle, b::Particle) = potential_energy(a.position, b.position)
function potential_energy(particles::AbstractVector{Particle})
    return sum(enumerate(particles[begin:(end - 1)])) do (i, particleᵢ)
        sum(particles[(i + 1):end]) do particleⱼ
            potential_energy(particleᵢ, particleⱼ)
        end
    end
end

function potential_gradient(𝐫)
    r = norm(𝐫)
    return (inv(r^8) / 2 - inv(r^14)) * 𝐫
end

kinetic_energy(particle::Particle) = 24 * sum(abs2, particle.velocity)
kinetic_energy(particles) = sum(kinetic_energy, particles)

total_energy(particles) = kinetic_energy(particles) + potential_energy(particles)

# function potential_directional_derivative(a::Particle, b::Particle)
#     𝐫 = a.position - b.position
#     u₀ = potential_energy(𝐫)
#     return function in_direction(Δ𝐫)
#         u₁ = potential_energy(𝐫 .+ Δ𝐫)
#         Δu = u₁ - u₀
#         ∇u = Δu ./ Δ𝐫
#         return dot(∇u, Δ𝐫) / norm(Δ𝐫)
#     end
# end
function potential_directional_derivative(a::Particle, b::Particle, δ=0.01)
    Δ𝐫 = (a.position - b.position) * δ
    ∇uₐ = potential_gradient(a.position)
    u₀ = potential_energy(a.position)
    u₁ = potential_energy(a.position .+ Δ𝐫)
    Δu = u₁ - u₀
    return Δu ./ Δ𝐫, ∇uₐ
end
function potential_directional_derivative(cell::Cell)
    U₀ = potential_energy(cell.particles)
    L = boxlength(cell)
    return function (particle::Particle, Δ𝐫)
        new_position = particle.position + Δ𝐫
        map!(Base.Fix2(mod, L), new_position, new_position)
        new_particle = Particle(new_position)
        particles = push!(list_neighbors(cell, particle), new_particle)
        U₁ = potential_energy(particles)
        Δu = U₁ - U₀
        return Δu ./ Δ𝐫
    end
end

"""
    Acceleration(a::Particle)(b::Particle)

Calculate the acceleration particle `b` induces on particle `a` (direction: from `b` to `a`).
"""
function Acceleration(a::Particle)
    return function (b::Particle)
        return Acceleration(potential_gradient(b.position .- a.position))
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
