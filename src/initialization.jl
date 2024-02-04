export Random, Even, Constant
export init_coordinates!, init_velocities!, relax!

abstract type CoordinatesDistribution end
struct Random <: CoordinatesDistribution end
struct Even <: CoordinatesDistribution end

abstract type VelocityDistribution end
struct Constant <: VelocityDistribution
    velocity::Velocity
end
struct MaxwellBoltzmann <: VelocityDistribution end

function init_coordinates!(particles, ::Cell, ::Even)
    for (particle, r) in zip(particles, vec(collect(Iterators.product(1:10, 1:10, 1:10))))
        particle.coordinates = collect(r) * 1.1
    end
    @assert unique(particles) == particles
    return particles
end
function init_coordinates!(particles, cell::Cell, ::Random)
    for particle in particles
        particle.coordinates = dimensions(cell) .* rand(3)
    end
    @assert unique(particles) == particles
    return particles
end

function init_velocities!(particles, dist::Constant)
    for particle in particles
        particle.velocity = dist.velocity
    end
    return particles
end
