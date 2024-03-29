using MolecularDynamics
using ProgressMeter: @showprogress
using Unitful: k, @u_str
using Plots

import MolecularDynamics: potential_energy, potential_gradient

ε = 119.8u"K" * k  # 119.8 kB
σ = 3.405u"angstrom"  # angstrom
mass = 39.948u"u"  # atomic mass unit
N = 864  # number of particles
target_T = 1385.42u"K"  # Corresponds to 1.069

u = LennardJones(ε, σ)
∇u = LennardJonesGradient(ε, σ)

potential_energy(particleᵢ::Particle, particleⱼ::Particle) = u(particleᵢ, particleⱼ)
potential_energy(particles) = u(particles)

potential_gradient(particleᵢ::Particle, particleⱼ::Particle) = ∇u(particleᵢ, particleⱼ)

particles = [
    Particle(
        mass,
        zeros(Coordinates{typeof(1.0u"angstrom")}),
        zeros(Velocity{typeof(1.0u"angstrom/fs")}),
    ) for _ in 1:N
];
cell = CubicCell(length(particles), 0.75u"angstrom^-3")
init_coordinates!(particles, cell, Even());
init_velocities!(particles, Constant(zeros(Velocity{typeof(1.0u"angstrom/s")})));
Δt = 1e-18u"s"

for _ in 1:50
    relax!(particles, cell, Δt, 1)
    println(potential_energy(particles))
end

integrator = VelocityVerlet()
trajectory = integrate!(particles, cell, Δt, 500, integrator);

@showprogress dt = 5 for round in 1:100
    error = abs(temperature(particles, k) - target_T)
    println("round: ", round, ", error: ", error)
    if error <= 20u"K"
        break
    else
        thermostat!(particles, k, VelocityRescaling(target_T))
        trajectory′ = integrate!(particles, cell, Δt, 800, integrator)
        append!(trajectory, trajectory′)
    end
end

energyplot(trajectory; yunit=u"eV")
savefig("e-t.pdf")

plot()
let indices = 100:100:1000
    for index in indices
        traceplot!(trajectory[12200:50:end], index)
    end
end

steps = (19000, 2000)
function plothist(step)
    velocities = norm.(extract(Velocity, trajectory, step))
    velocityhist(velocities)
    savefig("velocityhist_$step.pdf")
    return current()
end

temperatureplot(trajectory, k)
