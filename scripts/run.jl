using LinearAlgebra
using MolecularDynamics
using Plots
using Plots.Measures
using ProgressMeter

cummean(A) = cumsum(A) ./ (1:length(A))  # See https://discourse.julialang.org/t/cummean-cumall-and-cumany/46219/2

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

energyplot(logger)
savefig("e-t.pdf")

let indices = 100:100:1000
    L = box.side_length
    for index in indices
        particle = map(12200:50:length(logger.history)) do step
            extract(Particle, logger, step, index)
        end
        plot!(getcoordinates(particle)(); label="particle $index", palette=:tab20)
        scatter3d!(getcoordinates(particle)(); markersizes=1, markerstrokewidth=0, label="")
    end
    plot!(; legend=:none, framestyle=:box)
    xlims!(0, box.side_length)
    ylims!(0, box.side_length)
    zlims!(0, box.side_length)
    xlabel!(raw"x ($\sigma$)")
    ylabel!(raw"y ($\sigma$)")
    zlabel!(raw"z ($\sigma$)")
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
