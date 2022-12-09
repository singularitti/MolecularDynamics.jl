using LinearAlgebra: norm
using MolecularDynamics
using Plots

particles = [Particle(rand(3), rand(3)) for _ in 1:864];
box = CubicBox(length(particles), 0.75)
init!(particles, box, Random(), Uniform(zeros(Velocity)));
logger = Logger()
Δt = 0.02

for _ in 1:100
    relax!(particles, box, 1, Δt)
    println(potential_energy(particles))
end

# integrator = VelocityVerlet()
integrator = MetropolisHastings(1 / 1.069)
take_n_steps!(logger, particles, box, 500, Δt, 0.6, 0.05, integrator)

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
