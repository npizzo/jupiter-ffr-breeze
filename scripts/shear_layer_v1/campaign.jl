include("run.jl")
root,case=ARGS
case in ["$(a)_$(b)" for a in ("stable","neutral") for b in ("cyclonic","none","anticyclonic")]||error("Unknown case")
p=TOML.parsefile(joinpath(root,"preflight","PASS.toml"));p["status"]=="PASS"||error("Missing gate")
for (file,h) in p["source_hashes"]
 bytes2hex(sha256(read(file)))==h||error("Source changed: $file")
end
run_eddy(["configs/shear-layer-v1/$case.toml",joinpath(root,case),"GPU"])
