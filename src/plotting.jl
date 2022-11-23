using RecipesBase: RecipesBase, @userplot, @recipe

export coordplot

@userplot CoordPlot
@recipe function f(plot::CoordPlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    particles = plot.args[end]  # Extract `magnetization` from the args
    seriestype --> :scatter3d
    markersize --> 1
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xguide --> raw"$x / \sigma$"
    yguide --> raw"$y / \sigma$"
    zguide --> raw"$z / \sigma$"
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    legend_foreground_color --> nothing
    frame --> :box
    X, Y, Z = getcoordinates(particles)()
    return X, Y, Z
end
