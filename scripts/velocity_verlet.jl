using MolecularDynamics
using Unitful: k, @u_str
using Plots

import MolecularDynamics: potential_energy, potential_gradient

ε = u"J"(0.0103u"eV")  # 119.8 kB
σ = 3.405u"angstrom"  # angstrom
mass = 39.948u"u"  # atomic mass unit

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
    ) for _ in 1:864
];
cell = CubicCell(length(particles), 0.75u"angstrom^-3")
init!(particles, cell, Random(), Uniform(zeros(Velocity{typeof(1.0u"angstrom/fs")})));
Δt = 0.05u"fs"

for _ in 1:100
    relax!(particles, cell, 1, Δt)
    println(potential_energy(particles))
end

integrator = VelocityVerlet()
trajectory = run!(particles, cell, 500, Δt, integrator);

while abs(temperature(particles) - 1.069) >= 0.01
    trajectory2 = run!(particles, cell, 2, Δt, integrator)
    thermostat!(particles, VelocityRescaling(1.069))
end

energyplot(trajectory)
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

temperatureplot(trajectory)
