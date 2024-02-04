export Random, Even, Constant
export init_coordinates!, init_velocities!, relax!

abstract type CoordinatesDistribution end
struct Random <: CoordinatesDistribution end
struct Even <: CoordinatesDistribution end

abstract type VelocityDistribution end
struct Constant <: VelocityDistribution
    velocity::Velocity
end
struct MaxwellBoltzmann <: VelocityDistribution end

function generate_even_positions(n_particles, box_size)
    # Calculate initial grid dimensions based on cube root of number of particles
    cube_root = cbrt(n_particles)
    base_grid_size = floor(Int, cube_root)
    # Find the nearest grid dimensions that can accommodate all particles
    grid_dims = [base_grid_size, base_grid_size, base_grid_size]
    while prod(grid_dims) < n_particles
        # Increment the smallest dimension
        min_dim_index = argmin(grid_dims)
        grid_dims[min_dim_index] += 1
    end
    # Calculate spacing between particles along each axis
    spacing = box_size ./ grid_dims
    positions = []
    for i in 1:grid_dims[1]
        for j in 1:grid_dims[2]
            for k in 1:grid_dims[3]
                if length(positions) >= n_particles
                    break
                end
                # Calculate position within the box, centering particle in its grid cell
                pos = [
                    (i - 0.5) * spacing[1], (j - 0.5) * spacing[2], (k - 0.5) * spacing[3]
                ]
                push!(positions, pos)
            end
        end
    end
    return positions
end

function init_coordinates!(particles, cell::Cell, ::Even)
    positions = generate_even_positions(length(particles), dimensions(cell))
    for (particle, position) in zip(particles, positions)
        particle.coordinates = position
    end
    @assert unique(particles) == particles
    return particles
end
function init_coordinates!(particles, cell::Cell, ::Random)
    for particle in particles
        particle.coordinates = dimensions(cell) .* rand(3)
    end
    @assert unique(particles) == particles
    return particles
end

function init_velocities!(particles, dist::Constant)
    for particle in particles
        particle.velocity = dist.velocity
    end
    return particles
end
