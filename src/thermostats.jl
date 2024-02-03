export VelocityRescaling
export thermostat!

abstract type Thermostat end
struct VelocityRescaling{T} <: Thermostat
    target_temperature::T
end

function thermostat!(particles, thermostat::VelocityRescaling)
    t = temperature(particles)
    for particle in particles
        particle.velocity *= sqrt(thermostat.target_temperature / t)
    end
    return particles
end
