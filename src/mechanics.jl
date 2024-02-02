using ThreadsX

export LennardJones,
    LennardJonesGradient, potential_energy, kinetic_energy, potential_gradient

function potential_energy end

abstract type PairPotential end
(u::PairPotential)(a::Particle, b::Particle) = u(distance(a, b))
function (u::PairPotential)(particles)
    n = length(particles)
    return ThreadsX.sum(u(particles[i], particles[j]) for i in 1:(n - 1) for j in (i + 1):n)
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
(∇u::PairPotentialGradient)(particleᵢ::Particle, particleⱼ::Particle) =
    ∇u(particleᵢ.coordinates .- particleⱼ.coordinates)

struct LennardJonesGradient{S,T} <: PairPotentialGradient
    epsilon::S
    sigma::T
end
function (∇u::LennardJonesGradient)(𝐫)  # 𝐫 = 𝐫ᵢ - 𝐫ⱼ
    r = sqrt(sum(abs2, 𝐫))  # rᵢⱼ
    σr⁻¹ = ∇u.sigma / r
    return 48∇u.epsilon / ∇u.sigma^2 * 𝐫 * (σr⁻¹^8 / 2 - σr⁻¹^14)
end

kinetic_energy(particle::Particle) = sum(abs2, particle.velocity) * particle.mass / 2
kinetic_energy(particles) = sum(kinetic_energy, particles)

"""
    Force(particleᵢ::Particle)(particleⱼ::Particle)

Calculate the force particle `j` exerts on particle `i` (direction: from `j` to `i`).
"""
Force(particleᵢ::Particle) =
    (particleⱼ::Particle) -> Force(-potential_gradient(particleᵢ, particleⱼ))
function Force(i::Integer, particles, cell::Cell)
    neighbors = find_neighbors(i, particles, cell)
    return sum(Force(particles[i]), neighbors)
end
function Force(particle::Particle, particles, cell::Cell)
    neighbors = find_neighbors(particle, particles, cell)
    return sum(Force(particle), neighbors)
end

"""
    Acceleration(particleᵢ::Particle)(particleⱼ::Particle)

Calculate the acceleration of particle `i` due to particle `j` (direction: from `j` to `i`).
"""
Acceleration(particleᵢ::Particle) =
    (particleⱼ::Particle) -> Force(particleᵢ)(particleⱼ) / particleᵢ.mass
function Acceleration(i::Integer, particles, cell::Cell)
    neighbors = find_neighbors(i, particles, cell)
    return sum(Acceleration(particles[i]), neighbors)
end
function Acceleration(particle::Particle, particles, cell::Cell)
    neighbors = find_neighbors(particle, particles, cell)
    return sum(Acceleration(particle), neighbors)
end
