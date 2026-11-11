"""Read-only completion and provenance checks; no model execution."""
import hashlib
import json
from pathlib import Path
import subprocess
import tomllib

CASES=("zero512_b0","zero512_b162","modes512_b162","modes1024_b162")
SEED=Path("scripts/ferrel2026_seedtest_v1")
CONFIG=Path("configs/ferrel2026_seedtest_v1")


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()


def read(path):
    return tomllib.loads(Path(path).read_text())


def require(condition,message):
    if not condition:raise ValueError(message)


def source_gate(path,new=False):
    g=read(path);require(g["status"]=="PASS",f"Gate failed: {path}")
    folder=SEED if new else Path("scripts/ferrel2026_v3")
    configs=CONFIG if new else Path("configs/ferrel2026")
    expected_sources={p.name for p in folder.iterdir() if p.suffix in ((".jl",".py",".json",".sh") if new else (".jl",".py",".sh"))}
    if new:expected_sources|={f"slurm/ferrel2026-seedtest-v1-{part}.sbatch" for part in ("unit","512","1024")}
    require(set(g["source_hashes"])==expected_sources,"Incomplete source gate inventory")
    expected_configs=set(CASES) if new else {"b0","b162","boundary162","axis0","axis162"}
    require(set(g["configuration_hashes"])==expected_configs,"Incomplete config gate inventory")
    for name,digest in g["source_hashes"].items():
        file=Path(name) if new and name.startswith("slurm/") else folder/name
        require(sha(file)==digest,f"Source changed after validation: {file}")
    for name,digest in g["configuration_hashes"].items():
        require(sha(configs/f"{name}.toml")==digest,f"Configuration changed after validation: {name}")
    return g


def runtime(prior):
    path=prior/"runtime.toml";r=read(path)
    require(r["status"]=="PASS","Repaired runtime gate failed")
    env=Path("environments/ferrel2026_gpu")
    files={"original_project_sha256":Path("Project.toml"),"original_manifest_sha256":Path("Manifest.toml"),
           "project_sha256":env/"Project.toml","manifest_sha256":env/"Manifest.toml","preferences_sha256":env/"LocalPreferences.toml"}
    for key,file in files.items():require(sha(file)==r[key],f"Runtime fingerprint changed: {file}")
    require(set(r["repair_source_hashes"])=={"environment.jl","setup.jl","launch.jl"},"Incomplete runtime guard inventory")
    for name,digest in r["repair_source_hashes"].items():
        require(sha(Path("scripts/ferrel2026_gpu")/name)==digest,"Repaired runtime verifier changed")
    # Same package identity comparison as the frozen Julia verifier.
    def signature(path):
        return sorted((name,*(str(e.get(k,"")) for k in ("uuid","version","git-tree-sha1","repo-rev","repo-url","path")))
                      for name,entries in read(path)["deps"].items() for e in entries)
    require(signature(env/"Manifest.toml")==signature(Path("Manifest.toml")),"Resolved package identities differ")
    expected=read("Project.toml")
    expected.setdefault("extras",{}).update(CUDA_Runtime_jll="76a88914-d11a-5bdc-97e0-2f5a05c973a2",CUDA_Compiler_jll="d1e2174e-dfdc-576e-b43e-73b79eb1aca8")
    require(read(env/"Project.toml")==expected,"Unexpected repaired Project changes")
    return {"file":str(path),"sha256":sha(path),"status":"PASS"}


def successful_jobs(ids):
    ids=sorted(set(map(str,ids)))
    require(bool(ids) and all(j.isdigit() for j in ids),"Missing/malformed enclosing Slurm job IDs")
    output=subprocess.check_output(["sacct","-X","-j",",".join(ids),"-n","-P","-o","JobIDRaw,State,ExitCode"],text=True)
    states={p[0]:p[1:3] for line in output.splitlines() if len(p:=line.strip().split("|"))>=3}
    for job in ids:require(states.get(job)==["COMPLETED","0:0"],f"Enclosing job is not successful: {job}: {states.get(job)}")
    return {j:states[j] for j in ids}


def validate(root,gate_root,prior):
    initialization_path=gate_root/"INITIALIZATION_PASS.toml"
    initialization=source_gate(initialization_path,new=True)
    source_gate(prior/"PASS512.toml");source_gate(prior/"PASS1024.toml")
    runtime_record=runtime(prior)
    require(initialization["prior_pass512_sha256"]==sha(prior/"PASS512.toml"),"Changed original512 gate")
    require(initialization["runtime_marker_sha256"]==runtime_record["sha256"],"Runtime differs from initialization gate")
    reference=gate_root/"thermal_mode_reference.toml"
    require(initialization["thermal_mode_reference_sha256"]==sha(reference),"Changed frozen mode reference")
    audit=json.loads(subprocess.check_output(["python3",str(SEED/"verify_inheritance.py")],text=True))
    require(audit["status"]=="PASS_SOURCE_INHERITANCE_ONLY","Driver inheritance failed")
    ref=read(reference);sources={};records={};jobs=[initialization["slurm_job_id"]]
    for name in CASES:
        template=read(CONFIG/f"{name}.toml")
        folder=root/template["name"];cp=folder/"CASE_PASS.toml";g=source_gate(cp,new=True)
        require(g["case_name"]==name and g["simulated_hours"]==6,"Wrong case completion marker")
        links={"initialization_gate_sha256":initialization_path,"prior_pass512_sha256":prior/"PASS512.toml",
               "prior_pass1024_sha256":prior/"PASS1024.toml","runtime_marker_sha256":prior/"runtime.toml",
               "thermal_seed_record_sha256":folder/"thermal_seed.toml","metadata_sha256":folder/"metadata.toml"}
        for key,path in links.items():require(g[key]==sha(path),f"Broken case provenance: {name}/{key}")
        jobs.append(g["slurm_job_id"])
        cfg_path=root/"configs"/f"{name}.toml";cfg=read(cfg_path);m=read(folder/"metadata.toml")
        require(set(cfg)==set(template)|{"thermal_seed_reference_file","thermal_seed_reference_sha256"},"Unexpected runtime config keys")
        require(all(cfg[k]==v for k,v in template.items()),"Per-run config changed beyond recorded reference routing")
        require(Path(cfg["thermal_seed_reference_file"]).resolve()==reference.resolve(),"Wrong runtime reference file")
        require(cfg["thermal_seed_reference_sha256"]==sha(reference),"Wrong runtime reference hash")
        require(all(m[k]==v for k,v in cfg.items()),"Metadata disagrees with executed config")
        require(m["configuration_sha256"]==sha(cfg_path),"Executed config fingerprint mismatch")
        require(m["script_sha256"]==sha(SEED/"run.jl"),"Executed driver fingerprint mismatch")
        require(m["thermal_seed_source_sha256"]==sha(SEED/"thermal_seed.jl"),"Executed seed fingerprint mismatch")
        require(m["thermal_seed_record_sha256"]==sha(folder/"thermal_seed.toml"),"Seed record mismatch")
        require(m["status"]=="COMPLETE" and m["architecture"]=="GPU" and m["float_type"]=="Float64","Incomplete or wrong execution architecture")
        require(m["stop_seconds"]==m["simulated_seconds"]==21600 and m["segment_start_seconds"]==0 and m["pickup_checkpoint"]=="","Only complete fresh0–6h cases accepted")
        require(m["output_seconds"]==3600 and m["nz"]==96,"Unexpected frame cadence/vertical grid")
        actual=read(folder/"thermal_seed.toml")
        for key in ("reference_grid","reference_height_m","reference_seed","reference_raw_variance","reference_target_variance_K2","amplitude_parameter_K","envelope","phase_denominator","modes"):
            require(actual[key]==ref[key],f"Case does not share frozen seed reference: {name}/{key}")
        sources[name]=folder
        records[name]={"case_pass_file":str(cp),"case_pass_sha256":sha(cp),"enclosing_slurm_job_id":str(g["slurm_job_id"]),
                       "executed_configuration_sha256":sha(cfg_path),"metadata_sha256":sha(folder/"metadata.toml"),
                       "seed_record_sha256":sha(folder/"thermal_seed.toml")}
    states=successful_jobs(jobs)
    return sources,{"status":"PASS_PROVENANCE_AND_ENCLOSING_JOBS","initialization_gate_sha256":sha(initialization_path),
                    "initialization_slurm_job_id":str(initialization["slurm_job_id"]),"runtime":runtime_record,
                    "source_inheritance_audit":audit,"cases":records,"successful_jobs":states}
