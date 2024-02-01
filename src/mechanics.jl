export LennardJones,
    LennardJonesGradient, potential_energy, kinetic_energy, potential_gradient

function potential_energy end

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
    epsilon::S
    sigma::T
end
function (u::LennardJones)(r::Number)
    σ⁶r⁻⁶ = (u.sigma / r)^6
    σ¹²r⁻¹² = σ⁶r⁻⁶^2
    return 4u.epsilon * (σ¹²r⁻¹² - σ⁶r⁻⁶)
end

function potential_gradient end

abstract type PairPotentialGradient end
(∇u::PairPotentialGradient)(a::Particle, b::Particle) = ∇u(b.coordinates .- a.coordinates)

struct LennardJonesGradient{S,T} <: PairPotentialGradient
    epsilon::S
    sigma::T
end
function (∇u::LennardJonesGradient)(𝐫)  # 𝐫 = 𝐫ᵢ - 𝐫ⱼ
    r = norm(𝐫)
    σr⁻¹ = ∇u.sigma / r
    return 48∇u.epsilon / ∇u.sigma^2 * 𝐫 * (σr⁻¹^8 / 2 - σr⁻¹^14)
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
