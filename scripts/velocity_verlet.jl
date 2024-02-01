using MolecularDynamics
using Unitful: k, @u_str
using Plots

import MolecularDynamics: potential_energy, potential_gradient

ε = 0.0103u"eV"  # 119.8 kB
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
box = CubicBox(length(particles), 0.75u"angstrom^-3")
init!(particles, box, Random(), Uniform(zeros(Velocity{typeof(1.0u"angstrom/fs")})));
Δt = 1e-16u"fs"

for _ in 1:100
    relax!(particles, box, 1, Δt)
    println(potential_energy(particles))
end

integrator = VelocityVerlet()
logger = Logger{typeof(Δt),eltype(particles)}()
take_n_steps!(logger, particles, box, 500, Δt, integrator)

take_n_steps!(logger, particles, box, 400, Δt, integrator)
while abs(temperature(particles) - 1.069) >= 0.01
    take_n_steps!(logger, particles, box, 2, Δt, integrator)
    thermostat!(particles, VelocityRescaling(1.069))
end

energyplot(logger)
savefig("e-t.pdf")

plot()
let indices = 100:100:1000
    for index in indices
        traceplot!(logger.trajectory[12200:50:end], index)
    end
end

steps = (19000, 2000)
function plothist(step)
    velocities = norm.(extract(Velocity, logger, step))
    velocityhist(velocities)
    savefig("velocityhist_$step.pdf")
    return current()
end

temperatureplot(logger)
