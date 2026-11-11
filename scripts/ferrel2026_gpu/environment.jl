using TOML, SHA

const REPAIR_REPO=normpath(joinpath(@__DIR__,"..",".."))
const REPAIR_ENV=joinpath(REPAIR_REPO,"environments","ferrel2026_gpu")
repair_sha(path)=bytes2hex(sha256(read(path)))

function package_signature(path)
    manifest=TOML.parsefile(path)
    records=String[]
    for (name,entries) in manifest["deps"], entry in entries
        push!(records,join((name,get(entry,"uuid",""),get(entry,"version",""),
            get(entry,"git-tree-sha1",""),get(entry,"repo-rev",""),
            get(entry,"repo-url",""),get(entry,"path","")),"|"))
    end
    return sort!(records)
end

function expected_project(original)
    project=deepcopy(original)
    extras=get!(project,"extras",Dict{String,Any}())
    extras["CUDA_Runtime_jll"]="76a88914-d11a-5bdc-97e0-2f5a05c973a2"
    extras["CUDA_Compiler_jll"]="d1e2174e-dfdc-576e-b43e-73b79eb1aca8"
    return project
end

function repair_source_hashes()
    # Fixed runtime-guard files allow independent queue wrappers to be prepared
    # while setup runs. The scientific driver has its own complete v3 hash gate.
    Dict(name=>repair_sha(joinpath(@__DIR__,name)) for name in ("environment.jl","setup.jl","launch.jl"))
end

function verify_repaired_environment(preflight)
    marker=joinpath(preflight,"runtime.toml")
    report=TOML.parsefile(marker)
    report["status"]=="PASS" || error("Runtime setup gate did not pass")
    for (key,path) in (
        ("original_project_sha256",joinpath(REPAIR_REPO,"Project.toml")),
        ("original_manifest_sha256",joinpath(REPAIR_REPO,"Manifest.toml")),
        ("project_sha256",joinpath(REPAIR_ENV,"Project.toml")),
        ("manifest_sha256",joinpath(REPAIR_ENV,"Manifest.toml")),
        ("preferences_sha256",joinpath(REPAIR_ENV,"LocalPreferences.toml")))
        report[key]==repair_sha(path) || error("Runtime environment changed after validation: $key")
    end
    report["repair_source_hashes"]==repair_source_hashes() || error("Runtime wrapper source changed after setup")
    package_signature(joinpath(REPAIR_ENV,"Manifest.toml"))==package_signature(joinpath(REPAIR_REPO,"Manifest.toml")) || error("Resolved package versions changed")
    TOML.parsefile(joinpath(REPAIR_ENV,"Project.toml"))==expected_project(TOML.parsefile(joinpath(REPAIR_REPO,"Project.toml"))) || error("Unexpected project dependency changes")
    return report
end
