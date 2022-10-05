using LinearAlgebra
using MolecularDynamics
using Plots

p = [Particle(rand(3), rand(3)) for _ in 1:1000];
b = CubicBox(length(p), 0.75)
init!(p, b, Even(), Uniform(zeros(Velocity)));
lg = Logger(length(p))
Δt = 0.032

take_n_steps!(lg, p, b, 5000, Δt, VelocityVerlet())
# sort(norm.(extract(Velocity, lg, length(p))))

u = map(lg.history) do step
    potential_energy(step.snapshot)
end
k = map(lg.history) do step
    kinetic_energy(step.snapshot)
end
e = u .+ k

fig = plot(range(0; step=Δt, length=length(lg.history)), k; label="kinetic energy")
plot!(fig, range(0; step=Δt, length=length(lg.history)), u; label="potential energy")
plot!(fig, range(0; step=Δt, length=length(lg.history)), e; label="total energy")
xlabel!("time")
ylabel!("energy")
plot!(; framestyle=:box)
savefig("e-t.pdf")
