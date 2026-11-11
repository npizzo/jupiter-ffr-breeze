using TOML, SHA
const CASES = ("b0","b162","boundary162","axis0","axis162")
config_path(name)=joinpath("configs","ferrel2026","$(name).toml")
function source_hashes()
    Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if endswith(f,".jl") || endswith(f,".py") || endswith(f,".sh"))
end
configuration_hashes()=Dict(c=>bytes2hex(sha256(read(config_path(c)))) for c in CASES)
function verify_gate(path)
    gate=TOML.parsefile(path)
    gate["status"]=="PASS" || error("Gate failed: $path")
    gate["source_hashes"]==source_hashes() || error("Source changed after validation")
    gate["configuration_hashes"]==configuration_hashes() || error("Configuration changed after validation")
    return gate
end
function write_gate(path;details...)
    data=Dict{String,Any}("status"=>"PASS","source_hashes"=>source_hashes(),"configuration_hashes"=>configuration_hashes())
    for (k,v) in details;data[string(k)]=v;end
    open(path,"w") do io;TOML.print(io,data);end
end
