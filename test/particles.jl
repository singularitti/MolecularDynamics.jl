@testset "Test the `Particle` constructor" begin
    @test Particle([1, 2, 3], [4, 5, 6]) == Particle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    @test Particle([1//1, 2, 9//3], [4, 10//2, 6]) == Particle([1, 2.0, 3], [4.0, 5, 6])
    @test Particle([1.1, 2.3, 3], [4 / 3, 5, 6]) ==
        Particle([1.1, 2.3, 3.0], [4//3, 5.0, 6.0])
    @test_throws MethodError Particle([1, 2, 3])
end

@testset "Test function `isapprox` on `Particle`s" begin
    a = Particle([1, 2, 3], [4, 5, 6])
    b = Particle([1, 2, 3] .+ 1e-10, [4, 5, 6] .- 1e-10)
    @test isapprox(a, b)
end

@testset "Test `find_nearest_image` and `find_neighbors`" begin
    particles = [Particle(rand(3), rand(3)) for _ in 1:100]
    box = CubicBox(11)
    init!(particles, box)
    for particle in particles
        neighbors = find_neighbors(particle, particles, box)
        @test particle ∉ neighbors
        @test length(neighbors) == length(particles) - 1
        @test all(neighbors) do neighbor
            𝐫 = neighbor.position .- particle.position
            all(0 <= abs(rᵢ) <= Lᵢ for (rᵢ, Lᵢ) in zip(𝐫, boxsize(box)))
        end
    end
end
