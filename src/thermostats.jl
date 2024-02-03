export VelocityRescaling
export thermostat!

abstract type Thermostat end
struct VelocityRescaling{T} <: Thermostat
    target_temperature::T
end

function thermostat!(particles, kB, thermostat::VelocityRescaling)
    T = temperature(particles, kB)
    λ = sqrt(thermostat.target_temperature / T)
    # Parallel in-place modification of particle velocities
    ThreadsX.foreach(particles) do particle
        particle.velocity *= λ
    end
    return particles
end
