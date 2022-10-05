using LinearAlgebra
using MolecularDynamics
using Plots

particles = [Particle(rand(3), rand(3)) for _ in 1:1000];
box = CubicBox(length(particles), 0.75)
init!(particles, box, Even(), Uniform(zeros(Velocity)));
logger = Logger(length(particles))
Δt = 0.032

take_n_steps!(logger, particles, box, 5000, Δt, VelocityVerlet())
# sort(norm.(extract(Velocity, lg, length(p))))

U = map(logger.history) do step
    potential_energy(step.snapshot)
end
T = map(logger.history) do step
    kinetic_energy(step.snapshot)
end
E = U .+ T

fig = plot(
    range(0; step=Δt, length=length(logger.history)),
    T;
    label="kinetic energy",
    framestyle=:box,
)
plot!(fig, range(0; step=Δt, length=length(logger.history)), U; label="potential energy")
plot!(fig, range(0; step=Δt, length=length(logger.history)), E; label="total energy")
xlim!((0, simulation_time(logger)))
xlabel!("time")
ylabel!("energy")
savefig("e-t.pdf")
