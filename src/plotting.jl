using RecipesBase: RecipesBase, @recipe

@recipe function f(particles::AbstractVector{Particle})
    # Set a default value for an attribute with `-->`
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    X, Y, Z = getcoordinates(particles)()
    return X, Y, Z
end
