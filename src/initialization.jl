export Random, Even, Uniform, relax!

abstract type CoordinatesDistribution end
struct Random <: CoordinatesDistribution end
struct Even <: CoordinatesDistribution end

abstract type VelocityDistribution end
struct Uniform <: VelocityDistribution
    velocity::Velocity
end
struct MaxwellBoltzmann <: VelocityDistribution end

function init_coordinates!(particles, ::Box, ::Even)
    for (particle, r) in zip(particles, vec(collect(Iterators.product(1:10, 1:10, 1:10))))
        particle.coordinates = collect(r) * 1.1
    end
    @assert unique(particles) == particles
    return particles
end
function init_coordinates!(particles, box::Box, ::Random)
    for particle in particles
        particle.coordinates = boxsize(box) .* rand(3)
    end
    @assert unique(particles) == particles
    return particles
end

function init_velocities!(particles, dist::Uniform=Uniform(zeros(Velocity)))
    for particle in particles
        particle.velocity = dist.velocity
    end
    return particles
end

function init!(
    particles, box::Box, xdist::CoordinatesDistribution, vdist::VelocityDistribution
)
    init_coordinates!(particles, box, xdist)
    init_velocities!(particles, vdist)
    return particles
end

function damp!(particles, box, n, Δt)
    take_n_steps!(particles, box, n, Δt, VelocityVerlet())
    init_velocities!(particles)
    return particles
end

function relax!(particles, box, Δt)
    for (particle, 𝐟) in zip(particles, force(particles, box))
        particle.coordinates .+= 𝐟 .* Δt
        map!(Base.Fix2(mod, box.side_length), particle.coordinates, particle.coordinates)  # Move `𝐫` back to `0 - L` range
    end
    return particles
end
