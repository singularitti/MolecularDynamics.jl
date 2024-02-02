export VelocityRescaling
export thermostat!

abstract type Thermostat end
struct VelocityRescaling <: Thermostat
    target_temperature::Float64
end

function thermostat!(particles, thermostat::VelocityRescaling)
    t = temperature(particles)
    for particle in particles
        particle.velocity *= sqrt(thermostat.target_temperature / t)
    end
    return particles
end
