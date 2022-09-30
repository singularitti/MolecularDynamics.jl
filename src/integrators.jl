export VelocityVerlet
export take_one_step!, take_n_steps!, velocities, positions

abstract type Integrator end
struct VelocityVerlet <: Integrator end

struct StepTracker
    data::Matrix{Particle}
end

function take_half_step!(particle::Particle, 𝐚, Δt, L)
    particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt / 2)
    particle.position += particle.velocity * Δt  # 𝐫(t + Δt)
    mapinto!(Base.Fix2(mod, L), particle.position)  # Move `𝐫` back to `0 - L` range
    return particle
end

function take_one_step!(cell::Cell, Δt, ::VelocityVerlet)
    L = boxlength(cell)
    positions = map(eachparticle(cell), accelerations(cell)) do particle, 𝐚
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt / 2)
        position = particle.position + particle.velocity * Δt  # 𝐫(t + Δt)
        position = map(Base.Fix2(mod, L), position)  # Move `𝐫` back to `0 - L` range
    end
    for (particle, position) in zip(eachparticle(cell), positions)
        𝐚 = accelerations(cell, particle, position)  # 𝐚(t + Δt)
        particle.velocity += 𝐚 * Δt / 2  # 𝐯(t + Δt)
    end
    for (particle, position) in zip(eachparticle(cell), positions)
        particle.position = position
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

# See https://discourse.julialang.org/t/why-there-is-no-in-place-map-on-arrays/87956
mapinto!(f, arg) = map!(f, arg, arg)
