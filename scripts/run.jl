using LinearAlgebra
using MolecularDynamics
using Plots
using Plots.Measures
using ProgressMeter

particles = [Particle(rand(3), rand(3)) for _ in 1:864];
box = CubicBox(length(particles), 0.75)
init!(particles, box, Random(), UniformVelocity(zeros(Velocity)));
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

let indices = 100:100:1000
    for index in indices
        traceplot(logger.trajectory[12200:50:end], index)
    end
end

let steps = (500, 1000, 2000, 3000, 5000, 8000, 10000, 12000, 16000)
    plot_array = []
    for step in steps
        velocities = norm.(extract(Velocity, logger, step))
        push!(
            plot_array,
            histogram(
                velocities;
                title=string(step) * "th step",
                legend=:none,
                framestyle=:box,
                xlims=extrema(velocities),
                ylims=(0, Inf),
                xlabel=raw"velocity ($v$)",
                ylabel="frequency",
                titlefontsize=5,
                labelfontsize=5,
                tickfontsize=4,
                bottom_margin=-1mm,
                top_margin=-1mm,
                right_margin=0mm,
                right_margin=-1mm,
            ),
        )
    end
    plot(plot_array...)
    savefig("maxwell.pdf")
end

temperatureplot(logger)
