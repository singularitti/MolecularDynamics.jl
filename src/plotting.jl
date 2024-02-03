using RecipesBase: RecipesBase, @recipe, @userplot, @series

export energyplot, temperatureplot, traceplot, velocityhist

@recipe function f(particles::AbstractArray{Particle})
    seriestype --> :scatter3d
    markersize --> 1
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xguide --> raw"$x / \sigma$"
    yguide --> raw"$y / \sigma$"
    zguide --> raw"$z / \sigma$"
    guidefontsize --> 10
    tickfontsize --> 8
    legend --> :none
    frame --> :box
    palette --> :tab10
    X, Y, Z = getcoordinates(particles)()
    return X, Y, Z
end

@userplot TracePlot
@recipe function f(plot::TracePlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    trace, index = plot.args
    particles = map(trace) do step
        step.snapshot[index]
    end
    formatter --> :plain
    palette --> :tab20
    for type in (:scatter, :path)
        @series begin
            seriestype --> type
            particles
        end
    end
end

@userplot EnergyPlot
@recipe function f(plot::EnergyPlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    trajectory = plot.args[end]  # Extract `trace` from the args
    time = simulation_time(trajectory)
    U = map(trajectory) do step
        potential_energy(step.snapshot)
    end
    T = map(trajectory) do step
        kinetic_energy(step.snapshot)
    end
    E = U .+ T
    seriestype --> :path
    xlims --> extrema(time)
    xguide --> raw"simulation time ($t$)"
    yguide --> raw"energy ($\varepsilon$)"
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    legend_foreground_color --> nothing
    legend_position --> :topright
    frame --> :box
    palette --> :tab10
    grid --> nothing
    for (energy, label) in
        zip((T, U, E), ("kinetic energy", "potential energy", "total energy"))
        @series begin
            label --> label
            time, energy
        end
    end
end

@userplot TemperaturePlot
@recipe function f(plot::TemperaturePlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    trajectory, kB = plot.args
    time = simulation_time(trajectory)
    T = [temperature(step.snapshot, kB) for step in trajectory]
    seriestype --> :path
    xlims --> extrema(time)
    xguide --> raw"simulation time ($t$)"
    yguide --> raw"temperature ($T$)"
    guidefontsize --> 10
    tickfontsize --> 8
    legend --> :none
    frame --> :box
    palette --> :tab10
    grid --> nothing
    return time, T
end

@userplot VelocityHist
@recipe function f(plot::VelocityHist)
    velocities = plot.args[end]  # Extract `velocities` from the args
    seriestype --> :histogram
    xlims --> extrema(velocities)
    ylims --> (0, Inf)
    xlabel --> raw"velocity ($v$)"
    ylabel --> "frequency"
    guidefontsize --> 10
    tickfontsize --> 8
    legend --> :none
    framestyle --> :box
    grid --> nothing
    return velocities
end
