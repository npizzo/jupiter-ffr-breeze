using TOML,SHA
const CONFIG_ROOT="configs/ferrel35-width-v2"
const ORIGINAL_DAY200="data/cascade35-300-v1/cascade35_heated/segment_20"
source_hashes()=Dict(basename(f)=>bytes2hex(sha256(read(f))) for f in readdir(@__DIR__;join=true) if isfile(f))
config_hashes()=Dict(f=>bytes2hex(sha256(read(joinpath(CONFIG_ROOT,f)))) for f in ("diagnostic.toml","wide.toml"))
function gate(path;details...)
    d=Dict{String,Any}("status"=>"PASS","source_hashes"=>source_hashes(),"configuration_hashes"=>config_hashes(),"environment_hashes"=>Dict(f=>bytes2hex(sha256(read(f))) for f in ["Project.toml","Manifest.toml"]))
    for (k,v) in details;d[string(k)]=v;end
    open(path,"w") do io;TOML.print(io,d);end
end
function verify_gate(path)
    d=TOML.parsefile(path);d["status"]=="PASS" || error("Missing successful gate")
    d["source_hashes"]==source_hashes() || error("Campaign source changed")
    d["configuration_hashes"]==config_hashes() || error("Campaign config changed")
    for (f,h) in d["environment_hashes"];bytes2hex(sha256(read(f)))==h || error("Environment changed");end
    d
end
function writeconfig(root,name,cfg)
    mkpath(root);p=joinpath(root,name*".toml")
    open(p,"w") do io;TOML.print(io,cfg);end
    p
end
