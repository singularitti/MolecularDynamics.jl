using JSON: JSON

export save

function Base.Dict(particle::Particle)
    return Dict("position" => particle.position, "velocity" => particle.velocity)
end

function Base.print(io::IO, particle::Particle)
    JSON.print(io, Dict(particle))
    return nothing
end
function Base.print(io::IO, particles::AbstractVector{Particle})
    JSON.print(io, map(Dict, particles))
    return nothing
end

function save(filename, data)
    open(filename, "w") do io
        print(io, data)
    end
    return nothing
end
