using JLD2, Test
# Compare identical short, seeded CPU/GPU trajectories before chaotic divergence.
lastframe(root)=last(sort(filter(f->endswith(f,".jld2"),readdir(root;join=true))))
a=load(lastframe(ARGS[1])); b=load(lastframe(ARGS[2]))
@test a["t"] ≈ b["t"]
for key in ("theta", "u", "v", "w", "zeta")
    @test isapprox(a[key], b[key]; rtol=1e-5, atol=1e-9)
end
println("Short CPU/GPU trajectory comparison passed")
