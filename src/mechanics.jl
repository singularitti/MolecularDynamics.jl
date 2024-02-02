using ThreadsX

export LennardJones,
    LennardJonesGradient,
    Force,
    Acceleration,
    potential_energy,
    kinetic_energy,
    potential_gradient

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
kinetic_energy(particles) = ThreadsX.sum(kinetic_energy, particles)

struct Force{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end
"""
    Force(particleᵢ::Particle)(particleⱼ::Particle)

Calculate the force particle `j` exerts on particle `i` (direction: from `j` to `i`).
"""
function Force(particle::Particle)
    force(particle′::Particle) = Force(-potential_gradient(particle, particle′))
    function force(particles, cell::Cell)
        neighbors = generate_neighbors(particle, particles, cell)
        return sum(Force(particle), neighbors)
    end
    return force
end

struct Acceleration{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end
"""
    Acceleration(particleᵢ::Particle)(particleⱼ::Particle)

Calculate the acceleration of particle `i` due to particle `j` (direction: from `j` to `i`).
"""
function Acceleration(particle::Particle)
    acceleration(particle′::Particle) = Force(particle)(particle′) / particle.mass
    function acceleration(particles, cell::Cell)
        neighbors = generate_neighbors(particle, particles, cell)
        return sum(Acceleration(particle), neighbors)
    end
    return acceleration
end

similar_type(::Type{<:Force}, ::Type{T}, s::Size{(3,)}) where {T} = Force{T}
similar_type(::Type{<:Acceleration}, ::Type{T}, s::Size{(3,)}) where {T} = Acceleration{T}
