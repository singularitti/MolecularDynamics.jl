@testset "Test the potential gradient" begin
    𝐫 = 4:6
    δ = 0.00005
    u₀ = potential_energy(𝐫)
    for i in 1:3
        Δ𝐫 = zeros(3)
        Δ𝐫[i] = δ
        u₁ = potential_energy(𝐫 + Δ𝐫)
        @test (u₁ - u₀) / δ - 48 * potential_gradient(𝐫)[i] < 1e-8
    end
end
