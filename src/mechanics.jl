export potential_energy, kinetic_energy, total_energy, accelerations, potential_gradient

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

function potential_gradient(a::Particle, b::Particle)
    𝐫 = a.position - b.position
    u₀ = potential_energy(𝐫)
    return function (Δ𝐫)
        u₁ = potential_energy(𝐫 .+ Δ𝐫)
        Δu = u₁ - u₀
        return Δu ./ Δ𝐫
    end
end
function potential_gradient(cell::Cell)
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
