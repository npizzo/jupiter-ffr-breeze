# Environment repair only. Execute on an allocated Unity GPU node.
using Pkg, TOML, SHA
root=normpath(joinpath(@__DIR__,"..",".."))
env=joinpath(root,"environments","gamma2026_gpu")
original=joinpath(root,"environments","gamma2026","Manifest.toml")
manifest=joinpath(env,"Manifest.toml")
preferences=joinpath(env,"LocalPreferences.toml")
isfile(original) || error("Original resolved Manifest is missing")
isfile(preferences) || error("Explicit CUDA runtime/compiler preferences are missing")
function signature(path)
    d=TOML.parsefile(path)
    result=String[]
    for (name,records) in d["deps"], record in records
        push!(result,join((name,get(record,"uuid",""),get(record,"version",""),
            get(record,"git-tree-sha1",""),get(record,"repo-rev",""),
            get(record,"repo-url",""),get(record,"path","")),"|"))
    end
    sort!(result)
end
expected=signature(original)
if isfile(manifest)
    signature(manifest)==expected || error("Existing repaired Manifest differs from original package versions")
else
    cp(original,manifest)
end
Pkg.activate(env)
Pkg.instantiate()
Pkg.precompile()
signature(manifest)==expected || error("Repair changed resolved package versions")
Pkg.status()
# Preferences were present before CUDA import/precompilation.
using CUDA
CUDA.functional(true) || error("CUDA still unavailable after isolated runtime repair")
CUDA.versioninfo()
prefs=TOML.parsefile(preferences)
requested=VersionNumber(prefs["CUDA_Runtime_jll"]["version"])
actual=CUDA.runtime_version()
(actual.major,actual.minor)==(requested.major,requested.minor) || error("Selected CUDA runtime differs from explicit preference")
out=ARGS[1];mkpath(out)
metadata=Dict("status"=>"PASSED","repair_scope"=>"CUDA runtime/compiler preference and isolated cache only",
    "runtime_version"=>string(actual),"driver_version"=>string(CUDA.driver_version()),
    "gpu"=>CUDA.name(CUDA.device()),"original_manifest_sha256"=>bytes2hex(sha256(read(original))),
    "manifest_sha256"=>bytes2hex(sha256(read(manifest))),
    "preferences_sha256"=>bytes2hex(sha256(read(preferences))),
    "package_versions_unchanged"=>true,"compute_location"=>"Unity allocated GPU node")
open(joinpath(out,"runtime.toml"),"w") do io;TOML.print(io,metadata);end
