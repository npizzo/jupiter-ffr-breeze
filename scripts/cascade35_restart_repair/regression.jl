using JLD2,TOML,Test
include("../cascade35/restart_comparison.jl")
root="data/cascade35-preflight-v1"
a=load(joinpath(root,"restart_full/frame_00002.jld2"));b=load(joinpath(root,"restart_second/frame_00002.jld2"))
H=TOML.parsefile(joinpath(root,"restart_full/metadata.toml"))["height"]
report=Dict{String,Any}(streamfunction_restart_report(a,b,H))
@testset "Integrated diagnostic restart regression" begin
    @test report["passed"]
    @test streamfunction_restart_report(a,a,H)["passed"]
    # Deliberately inconsistent derived diagnostic must be rejected even if
    # a permissive absolute tolerance on the depth-normalized value would pass.
    corrupted=copy(b);corrupted["mass_streamfunction"]=copy(b["mass_streamfunction"])
    corrupted["mass_streamfunction"][2,2]+=100report["propagated_bound"]
    @test !streamfunction_restart_report(a,corrupted,H)["propagated_error_pass"]
    # A changed native flow must fail the unchanged native-flux criterion.
    badflow=copy(b);badflow["native_rhov"]=copy(b["native_rhov"])
    badflow["native_rhov"][2,2]+=1e-6
    @test !streamfunction_restart_report(a,badflow,H)["native_flux_pass"]
end
out="data/cascade35-restart-audit-v1";report["status"]="PASS"
open(joinpath(out,"regression.toml"),"w") do io;TOML.print(io,report);end
println(report)
