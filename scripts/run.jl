using LinearAlgebra
using MolecularDynamics
using Plots
using Plots.Measures
using ProgressMeter

particles = [Particle(rand(3), rand(3)) for _ in 1:864];
box = CubicBox(length(particles), 0.75)
init!(particles, box, Random(), Uniform(zeros(Velocity)));
logger = Logger()
Δt = 0.0001

for _ in 1:100
    relax!(particles, box, Δt)
    # damp!(particles, box, 1, Δt)
    println(potential_energy(particles))
end

# integrator = VelocityVerlet()
integrator = MetropolisHastings(1 / 1.069)

take_n_steps!(logger, particles, box, 400, Δt, integrator)
while abs(temperature(particles) - 1.069) >= 0.01
    apply_coupling!(particles, VelocityRescaling(1.069))
    take_n_steps!(logger, particles, box, 2, Δt, integrator)
end

U = progress_map(logger.trajectory) do step
    potential_energy(step.snapshot)
end
T = map(logger.trajectory) do step
    kinetic_energy(step.snapshot)
end
E = U .+ T

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
