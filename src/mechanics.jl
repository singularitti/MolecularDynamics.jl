export LennardJones, kinetic_energy, potential_gradient

abstract type PairPotential end
(u::PairPotential)(a::Particle, b::Particle) = u(distance(a, b))
function (u::PairPotential)(particles)
    return sum(enumerate(particles[begin:(end - 1)])) do (i, particleᵢ)
        sum(particles[(i + 1):end]) do particleⱼ
            u(particleᵢ, particleⱼ)
        end
    end
end

struct LennardJones{S,T} <: PairPotential
    ε::S
    σ::T
end
function (uₗⱼ::LennardJones)(r::Number)
    σ⁶r⁻⁶ = (uₗⱼ.σ / r)^6
    σ¹²r⁻¹² = σ⁶r⁻⁶^2
    return 4uₗⱼ.ε * (σ¹²r⁻¹² - σ⁶r⁻⁶)
end

abstract type PairPotentialGradient end
(∇u::PairPotentialGradient)(a::Particle, b::Particle) = ∇u(b.coordinates .- a.coordinates)

struct LennardJonesGradient{S,T} <: PairPotentialGradient
    ε::S
    σ::T
end
function (::LennardJonesGradient)(𝐫)
    r = norm(𝐫)
    return 𝐫 * (inv(r)^8 / 2 - inv(r)^14)
end

kinetic_energy(particle::Particle) = sum(abs2, particle.velocity) * particle.mass / 2
kinetic_energy(particles) = sum(kinetic_energy, particles)

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
function Acceleration(i::Integer, particles, box::Box)
    neighbors = find_neighbors(i, particles, box)
    return sum(Acceleration(particles[i]), neighbors)
end
function Acceleration(particle::Particle, particles, box::Box)
    neighbors = find_neighbors(particle, particles, box)
    return sum(Acceleration(particle), neighbors)
end
function Acceleration(particles, box)
    return map(particles) do particle
        Acceleration(particle, particles, box)
    end
end
