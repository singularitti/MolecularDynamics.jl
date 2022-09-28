using StaticArrays: SVector, MVector

mutable struct Particle
    r::MVector{3,Float64}
    v::MVector{3,Float64}
    Particle(r, v) = new(r, v)
    Particle() = new()  # Incomplete initialization
end

function distance(particle::Particle, particle′::Particle)
    return sqrt(sum(abs2, particle.r - particle′.r))
end

