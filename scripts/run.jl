using LaTeXStrings
using LinearAlgebra
using MolecularDynamics
using Plots
using ProgressMeter

particles = [Particle(rand(3), rand(3)) for _ in 1:1000];
box = CubicBox(length(particles), 0.75)
init!(particles, box, Even(), Uniform(zeros(Velocity)));
logger = Logger(length(particles))
Δt = 0.032

take_n_steps!(logger, particles, box, 4000, Δt, VelocityVerlet())
apply_coupling!(particles, VelocityRescaling(1.069))

U = progress_map(logger.history) do step
    potential_energy(step.snapshot)
end
# append!(
#     U,
#     progress_map((length(U) + 1):length(logger.history)) do i
#         potential_energy(logger.history[i].snapshot)
#     end,
# )
T = map(logger.history) do step
    kinetic_energy(step.snapshot)
end
E = U .+ T

fig = plot(simulation_time(logger), T; label="kinetic energy", framestyle=:box)
plot!(fig, simulation_time(logger), U; label="potential energy")
plot!(fig, simulation_time(logger), E; label="total energy")
xlims!(extrema(simulation_time(logger)))
plot!(; legend=:left)
xlabel!(L"simulation time ($t$)")
ylabel!(L"energy ($\varepsilon$)")
savefig("e-t.pdf")

histogram(norm.(extract(Velocity, logger, length(particles))))

plot(
    simulation_time(logger),
    temperature.(step.snapshot for step in logger.history);
    label="temperature",
    framestyle=:box,
)
xlims!(extrema(simulation_time(logger)))
plot!(; legend=:left)
xlabel!(L"time ($t$)")
ylabel!(L"temperature ($T$)")
