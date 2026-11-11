# Unity allocated GPU node only; no model configuration or original environment edits.
using Pkg
include("environment.jl")
preflight=abspath(only(ARGS));mkpath(preflight)
isfile(joinpath(preflight,"runtime.toml")) && error("Refusing to overwrite a runtime gate")
isfile(joinpath(REPAIR_ENV,"READY.toml")) && error("Isolated environment already initialized; do not rebuild an active environment")
preferences=joinpath(REPAIR_ENV,"LocalPreferences.toml")
prefs=TOML.parsefile(preferences)
prefs["CUDA_Runtime_jll"]["version"]==prefs["CUDA_Compiler_jll"]["version"]=="12.2" || error("Expected matched runtime/compiler preferences before precompilation")

original_project=joinpath(REPAIR_REPO,"Project.toml")
original_manifest=joinpath(REPAIR_REPO,"Manifest.toml")
original_project_sha=repair_sha(original_project)
original_manifest_sha=repair_sha(original_manifest)
expected_signature=package_signature(original_manifest)
function preserve_copy(src,dst)
    if isfile(dst)
        repair_sha(src)==repair_sha(dst) || error("Existing preserved original differs: $dst")
    else
        cp(src,dst)
    end
end
# Keep byte-exact originals before deriving an effective Project with JLL extras.
preserve_copy(original_project,joinpath(REPAIR_ENV,"OriginalProject.toml"))
preserve_copy(original_manifest,joinpath(REPAIR_ENV,"OriginalManifest.toml"))
project=joinpath(REPAIR_ENV,"Project.toml")
expected=expected_project(TOML.parsefile(original_project))
if isfile(project)
    TOML.parsefile(project)==expected || error("Existing effective Project differs")
else
    cp(original_project,project)
    open(project,"w") do io;TOML.print(io,expected);end
end
manifest=joinpath(REPAIR_ENV,"Manifest.toml")
if isfile(manifest)
    package_signature(manifest)==expected_signature || error("Existing Manifest package versions differ")
else
    cp(original_manifest,manifest)
end

Pkg.activate(REPAIR_ENV)
Pkg.instantiate()
Pkg.precompile()
package_signature(manifest)==expected_signature || error("Runtime setup changed resolved package versions")
TOML.parsefile(project)==expected || error("Runtime setup changed Project beyond preference extras")
repair_sha(original_project)==original_project_sha || error("Original Project changed during setup")
repair_sha(original_manifest)==original_manifest_sha || error("Original Manifest changed during setup")

using CUDA
CUDA.functional(true) || error("Explicit runtime did not initialize CUDA")
CUDA.versioninfo()
runtime=CUDA.runtime_version()
(runtime.major,runtime.minor)==(12,2) || error("Unexpected CUDA runtime: $runtime")
# Exercise allocation, compilation, execution and reduction on the assigned GPU.
probe=CUDA.ones(Float64,64)
probe .= 2 .* probe
CUDA.synchronize()
sum(probe)==128.0 || error("GPU runtime/compiler execution probe failed")
report=Dict{String,Any}(
    "status"=>"PASS","runtime_version"=>string(runtime),
    "compiler_preference"=>prefs["CUDA_Compiler_jll"]["version"],
    "driver_version"=>string(CUDA.driver_version()),"gpu_name"=>CUDA.name(CUDA.device()),
    "original_project_sha256"=>original_project_sha,"original_manifest_sha256"=>original_manifest_sha,
    "project_sha256"=>repair_sha(project),"manifest_sha256"=>repair_sha(manifest),
    "preferences_sha256"=>repair_sha(preferences),"repair_source_hashes"=>repair_source_hashes(),
    "package_versions_unchanged"=>true,"gpu_execution_probe"=>"Float64 allocation, broadcast kernel and reduction passed",
    "science_source"=>"Unchanged scripts/ferrel2026_v3; its own source/config gates remain mandatory",
    "scope"=>"Isolated CUDA runtime/compiler selection and writable depot only")
for path in (joinpath(preflight,"runtime.toml"),joinpath(REPAIR_ENV,"READY.toml"))
    open(path,"w") do io;TOML.print(io,report);end
end
println("Isolated Ferrel CUDA setup PASS; root package versions unchanged")
