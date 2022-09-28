export VelocityVerlet
export take_one_step!, take_n_steps!, velocities, positions

abstract type Integrator end
struct VelocityVerlet <: Integrator end

struct StepTracker
    data::Matrix{Particle}
end

function take_one_step!(cell::Cell, i, Δt, ::VelocityVerlet)
    particle = cell.particles[i]
    particle.velocity += accelerationof(cell, i) * Δt / 2  # 𝐯(t + Δt / 2)
    particle.position += particle.velocity * Δt  # 𝐫(t + Δt)
    map!(Base.Fix2(mod, boxlength(cell)), particle.position, particle.position)
    𝐚 = accelerationof(cell, i)  # 𝐚(t + Δt)
    particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt)
    return cell
end
function take_one_step!(cell::Cell, Δt, ::VelocityVerlet)
    for i in eachindex(cell.particles)
        take_one_step!(cell, i, Δt, VelocityVerlet())
    end
    return cell
end

function take_n_steps!(cell::Cell, n, Δt, ::VelocityVerlet)
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
