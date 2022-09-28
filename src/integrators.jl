export VelocityVerlet
export damp!, take_one_step!, take_n_steps!, velocities, positions

abstract type Integrator end
struct VelocityVerlet <: Integrator end

struct StepTracker
    data::Matrix{Particle}
end

function damp!(particles, n, Δt)
    take_n_steps!(particles, n, Δt, VelocityVerlet())
    init_velocities!(particles)
    return particles
end

function take_one_step!(cell::SimulationCell, i, Δt, ::VelocityVerlet)
    particles = cell.particles
    particles[i].velocity += accelerationof(cell, i) * Δt / 2  # 𝐯(t + Δt / 2)
    particles[i].position += particles[i].velocity * Δt  # 𝐫(t + Δt)
    map!(
        Base.Fix2(apply_pbc, boxlength(cell)), particles[i].position, particles[i].position
    )
    𝐚 = accelerationof(cell, i)  # 𝐚(t + Δt)
    particles[i].velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt)
    return particles
end
function take_one_step!(cell::SimulationCell, Δt, ::VelocityVerlet)
    for i in eachindex(cell.particles)
        take_one_step!(cell, i, Δt, VelocityVerlet())
    end
    return cell
end

function take_n_steps!(cell::SimulationCell, n, Δt, ::VelocityVerlet)
    data = Matrix{Particle}(undef, particlenumber(cell), n)
    for i in 1:n
        # Must use `deepcopy`!
        take_one_step!(cell, Δt, VelocityVerlet())
        data[:, i] = deepcopy(cell.particles)
    end
    return StepTracker(data)
end

function velocities(tracker::StepTracker)
    return map(tracker.data) do particle
        particle.velocity
    end
end

function positions(tracker::StepTracker)
    return map(tracker.data) do particle
        particle.position
    end
end
