@testset "Test the potential gradient" begin
    š« = 4:6
    Ī“ = 0.00005
    uā = potential_energy(š«)
    for i in 1:3
        Īš« = zeros(3)
        Īš«[i] = Ī“
        uā = potential_energy(š« + Īš«)
        @test (uā - uā) / Ī“ - 48 * potential_gradient(š«)[i] < 1e-8
    end
end
