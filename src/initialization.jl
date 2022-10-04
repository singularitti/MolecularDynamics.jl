function init_coordinates!(particles, box::Box)
    # for particle in particles
    #     particle.coordinates = boxsize(box) .* rand(3)
    # end
    for (particle, r) in zip(particles, vec(collect(Iterators.product(1:10, 1:10, 1:10))))
        particle.coordinates = collect(r) * 1.1
    end
    @assert unique(particles) == particles
    return particles
end

function init_velocities!(particles)
    for particle in particles
        particle.velocity = zeros(Velocity)
    end
    return particles
end

function init!(particles, box::Box)
    init_coordinates!(particles, box)
    init_velocities!(particles)
    return particles
end

function damp!(particles, box, n, Δt)
    take_n_steps!(particles, box, n, Δt, VelocityVerlet())
    init_velocities!(particles)
    return particles
end
