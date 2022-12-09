export VelocityRescaling
export thermostat!

abstract type Thermostat end
struct VelocityRescaling <: Thermostat
    desired_temperature::Float64
end

function thermostat!(particles, thermostat::VelocityRescaling)
    t = temperature(particles)
    for particle in particles
        particle.velocity *= sqrt(thermostat.desired_temperature / t)
    end
    return particles
end
