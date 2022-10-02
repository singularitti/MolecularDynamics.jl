using RecipesBase: @recipe

@recipe function plot(particles::AbstractVector{Particle})
    # Set a default value for an attribute with `-->`
    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"
    X, Y, Z = ntuple(_ -> Vector{Float64}(undef, length(particles)), 3)
    for (i, particle) in enumerate(particles)
        X[i], Y[i], Z[i] = particle.position
    end
    return X, Y, Z
end
