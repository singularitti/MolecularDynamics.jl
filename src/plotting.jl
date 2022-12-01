using RecipesBase: RecipesBase, @recipe

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
    legendfontsize --> 8
    legend_foreground_color --> nothing
    frame --> :box
    X, Y, Z = getcoordinates(particles)()
    return X, Y, Z
end
