using TOML, SHA
const HEATS = (0, 54, 162)
config_path(heat)=joinpath("configs","convection2026","f$(heat).toml")
function source_hashes()
    Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if endswith(f,".jl") || endswith(f,".py") || endswith(f,".sh"))
end
configuration_hashes()=Dict(string(h)=>bytes2hex(sha256(read(config_path(h)))) for h in HEATS)
function verify_gate(path)
    gate=TOML.parsefile(path)
    gate["status"]=="PASS" || error("Gate did not pass: $path")
    gate["source_hashes"]==source_hashes() || error("Source changed after validation")
    gate["configuration_hashes"]==configuration_hashes() || error("Configuration changed after validation")
    return gate
end
function write_gate(path; details...)
    data=Dict{String,Any}("status"=>"PASS","source_hashes"=>source_hashes(),"configuration_hashes"=>configuration_hashes())
    for (k,v) in details; data[string(k)]=v; end
    open(path,"w") do io; TOML.print(io,data);end
end
