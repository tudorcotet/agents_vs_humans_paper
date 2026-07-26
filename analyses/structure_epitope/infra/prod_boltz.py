import json, subprocess, time, os, glob, gemmi, sys, yaml

data = json.load(open("/root/designs_all.json"))
TARGET_A3M = "/root/target.a3m"
YDIR = "/root/prod_yamls"; OUT = "/root/prod_out"
os.makedirs(YDIR, exist_ok=True)

# build 141 YAMLs: chain A target with precomputed MSA; chain B binder single-seq
for did, v in data["designs"].items():
    name = f"design_{int(did):03d}"
    ydict = {"version":1, "sequences":[
        {"protein":{"id":"A","sequence":data["target"],"msa":TARGET_A3M}},
        {"protein":{"id":"B","sequence":v["seq"],"msa":"empty"}},
    ]}
    open(f"{YDIR}/{name}.yaml","w").write(yaml.dump(ydict, sort_keys=False))
print(f"wrote {len(data['designs'])} YAMLs", flush=True)

# batch fold (no --use_msa_server; MSA precomputed). single model, seed 42.
t0=time.time()
p = subprocess.run(f"boltz predict {YDIR} --out_dir {OUT} --seed 42 --write_full_pae --override",
                   shell=True, capture_output=True, text=True)
dt=time.time()-t0
print(f"boltz batch exit={p.returncode} wall={dt:.0f}s", flush=True)
print((p.stdout or "")[-600:]); print((p.stderr or "")[-600:], flush=True)

# integrity check every output: chain B == input binder
def bchain(cif):
    m=gemmi.read_structure(cif)[0]
    chs=list(m)
    tgt=next((c for c in chs if len([r for r in c])==175), None)
    bnd=next(c for c in chs if c is not tgt)
    return gemmi.one_letter_code([r.name for r in bnd])
res=[]
for did,v in data["designs"].items():
    name=f"design_{int(did):03d}"
    cifs=[c for c in glob.glob(f"{OUT}/**/predictions/**/*.cif", recursive=True) if name in c and "model_0" in c]
    ok=bool(cifs); integ=None
    if ok:
        try: integ = (bchain(cifs[0]).replace("X","")==v["seq"])
        except Exception as e: integ=f"err:{e}"
    res.append({"design_id":int(did),"ok":ok,"integrity":integ,"cif":cifs[0] if cifs else None})
json.dump(res, open("/root/prod_integrity.json","w"), indent=2)
nok=sum(1 for r in res if r["ok"]); nint=sum(1 for r in res if r["integrity"] is True)
print(f"PROD_DONE: {nok}/141 folded, {nint}/141 integrity-pass", flush=True)
bad=[r["design_id"] for r in res if r["integrity"] is not True]
if bad: print("INTEGRITY-FAIL or missing:", bad, flush=True)
