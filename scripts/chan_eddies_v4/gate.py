import os,sys,json,hashlib,subprocess,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; HERE=Path(__file__).resolve().parent
os.chdir(ROOT)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(x,indent=2)+"\n")
def manifest():
 files=list(HERE.glob("*.jl"))+list(HERE.glob("*.py"))
 for p in ["slurm/chan-eddies-v4-screen.sbatch","slurm/chan-eddies-v4-case3.sbatch","slurm/chan-eddies-v4-case4.sbatch","slurm/chan-eddies-v4-cpu.sbatch","configs/chan-eddies-v4/case3.toml","configs/chan-eddies-v4/case4.toml","scripts/cuda_repair/environment.sh","Project.toml","Manifest.toml","environments/cascade-cuda122/Project.toml"]: files.append(ROOT/p)
 for p in [ROOT,ROOT/"environments/cascade-cuda122"]:
  files+=list(p.glob("Manifest*.toml"))+list(p.glob("LocalPreferences.toml"))
 return {"files":{str(p.relative_to(ROOT)):sha(p) for p in files},"julia":sha(os.environ["JULIA_BIN"]),"env":{k:os.environ.get(k,"") for k in ["JULIA_LOAD_PATH","JULIA_DEPOT_PATH","JULIA_NUM_THREADS","OPENBLAS_NUM_THREADS"]}}
def verify(p):
 r=json.loads(Path(p).read_text());assert r["status"]=="PASS" and r["manifest"]==manifest();assert sha(r["result"])==r["result_sha256"]
 print("EXACT_DRIVER_RECEIPT_VERIFIED",flush=True)
mode=sys.argv[1]
if mode=="preflight":
 out=Path(sys.argv[2]).resolve();before=manifest();assert not out.exists()
 subprocess.run([os.environ["JULIA_BIN"],"--startup-file=no","-e",'function check(x); if x isa Expr; x.head in (:error,:incomplete) && error(string(x)); foreach(check,x.args); end; end; for f in ARGS; check(Meta.parseall(read(f,String))); end']+[str(p) for p in HERE.glob("*.jl")],check=True)
 subprocess.run([os.environ["JULIA_BIN"],"--project=.",str(HERE/"driver.jl"),"preflight",str(out)],check=True)
 assert before==manifest();p=out/"DRIVER_PASS.toml";assert p.is_file()
 write(out/"LAUNCH_PASS.json",{"status":"PASS","manifest":before,"result":str(p),"result_sha256":sha(p)})
elif mode=="verify":verify(sys.argv[2])
elif mode=="submit":
 receipt,case=sys.argv[2:4];assert case in ["screen","case3","case4"];verify(receipt);state=Path(receipt).with_name("SUBMITTED_"+case+".json")
 with open(str(state)+".lock","a") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if state.exists():raise RuntimeError("Submission already attempted; inspect scheduler before retrying")
  write(state,{"status":"SUBMITTING"})
  p=subprocess.run(["sbatch","--parsable","slurm/chan-eddies-v4-"+case+".sbatch",receipt]+sys.argv[4:],capture_output=True,text=True,check=True)
  job=p.stdout.strip().split(";")[0];assert job.isdigit();write(state,{"status":"SUBMITTED","job_id":job});print("SUBMITTED",case,job)
else:raise RuntimeError(mode)
