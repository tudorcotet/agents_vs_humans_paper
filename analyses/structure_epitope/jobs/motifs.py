"""
Phase 4: binding-motif analysis, Levels 1-3 (addendum Part B). Protenix + Chai binder chains.
L1 interface architecture (which SSEs contact target; beta-augmentation flag).
L2 contact-pair chemistry (binder-restype x target-region), log-odds binders vs non-binders.
L3 convergent contact motifs across independent teams.
"""
import warnings, re, numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from scipy import stats
import pydssp, torch
warnings.filterwarnings("ignore")
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
TARGET_LEN = 175
def cls(a):
    if a in "KR": return "pos"
    if a in "DE": return "neg"
    if a in "STNQHYCG": return "polar"
    if a in "FWY": return "aromatic"
    return "aliphatic"  # AVLIMP
def region(u):  # target UniProt -> functional region
    if u in {44,69,70,71,74,89}: return "tip"
    if u in {46,47,62,76,77,87,122,123}: return "basic"
    if 133<=u<=174: return "stalk"
    if u<19 or u>132: return "other"
    return "apex_other"  # IgSF non-tip/basic
def one(res):
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    try: return index_to_one(three_to_index(res.resname))
    except Exception: return "X"

def analyze(design_id, pred):
    if pred == "af2m":
        p = ROOT / f"data/structures/{pred}/design_{design_id:03d}.pdb"; parser = PDBParser(QUIET=True)
    else:
        p = ROOT / f"data/structures/{pred}/design_{design_id:03d}.cif"; parser = MMCIFParser(QUIET=True)
    if not p.exists(): return None
    model = parser.get_structure("m", str(p))[0]
    chs = list(model.get_chains())
    tgt = next((c for c in chs if len([r for r in c if r.id[0]==" "])==TARGET_LEN), None)
    bnd = next(c for c in chs if c is not tgt)
    bnd_res = [r for r in bnd if r.id[0]==" "]
    # binder SS via pydssp
    coords=[]; idx=[]
    for i,r in enumerate(bnd_res):
        try: coords.append([r["N"].coord,r["CA"].coord,r["C"].coord,r["O"].coord]); idx.append(i)
        except KeyError: pass
    ss_full = ["-"]*len(bnd_res)
    if len(coords)>=5:
        s = pydssp.assign(torch.tensor(np.array(coords,dtype=np.float32)), out_type="c3")
        for j,i in enumerate(idx): ss_full[i]=s[j]
    # inter-chain contacts (4.5A heavy)
    tgt_heavy=[a for a in tgt.get_atoms() if a.element!="H"]
    ns=NeighborSearch(tgt_heavy)
    iface_bidx=set(); pairs=[]  # (binder_type, target_uniprot, target_type)
    for i,r in enumerate(bnd_res):
        contacted=False
        for a in r:
            if a.element=="H": continue
            for x in ns.search(a.coord,4.5,level="A"):
                tr=x.get_parent(); tu=tr.id[1]
                if 19<=tu<=174:
                    pairs.append((one(r), tu, one(tr))); contacted=True
        if contacted: iface_bidx.add(i)
    # L1 interface SSE composition
    iface_ss=[ss_full[i] for i in iface_bidx]
    n=len(iface_ss) or 1
    hf=iface_ss.count("H")/n; ef=iface_ss.count("E")/n; lf=iface_ss.count("-")/n
    # beta-augmentation: binder strand residue contacting target strand region (approx: target 100-131 F/G strands)
    beta_aug = any(ss_full[i]=="E" for i in iface_bidx) and ef>0.15
    return dict(design_id=design_id, pred=pred, iface_helix=round(hf,3), iface_strand=round(ef,3),
                iface_loop=round(lf,3), n_iface=len(iface_bidx), beta_aug=beta_aug), pairs

meta = pd.read_parquet(ROOT/"analyses/structure_epitope/data/designs_canonical.parquet")
teaminfo = meta.set_index("design_id")[["cohort","is_hit","p_expressed","team","kd_arith_mean_nM_all"]]
L1=[]; allpairs={}
for did in meta.design_id:
    for pred in ["protenix","chai","esmfold2","boltz2","af2m"]:
        out=analyze(did,pred)
        if out is None: continue
        d,pairs=out; L1.append(d)
        if pred in ("protenix","chai"):          # L2/L3 contact enrichment stays Protenix+Chai
            allpairs.setdefault(did,[]).extend(pairs)
L1df=pd.DataFrame(L1)
L1df.to_csv(RES/"motifs_L1_architecture_perpredictor.csv",index=False)   # per-predictor helix/strand/loop
avg_src=L1df[L1df.pred.isin(["protenix","chai"])]                          # §4 architecture avg = Protenix+Chai
l1=avg_src.groupby("design_id").mean(numeric_only=True).reset_index()
l1["beta_aug"]=avg_src.groupby("design_id")["beta_aug"].any().values
l1=l1.merge(teaminfo,on="design_id")
# architecture class
def arch(r):
    if r.iface_loop>0.5: return "loop-dominated"
    if r.iface_strand>0.25: return "strand"
    if r.iface_helix>0.5: return "helix"
    return "mixed"
l1["architecture"]=l1.apply(arch,axis=1)
l1.to_csv(RES/"motifs_L1_architecture.csv",index=False)
print("=== L1 interface architecture (expressed): class x binder ===")
exp=l1[l1.p_expressed==True]
ct=pd.crosstab(exp.architecture,exp.is_hit)
print(ct.to_string())
print("beta-augmentation flag: n designs =",int(l1.beta_aug.sum()),"| among binders:",int(l1[l1.is_hit==True].beta_aug.sum()))
print(f"interface helix fraction: binder med {exp[exp.is_hit==True].iface_helix.median():.2f} vs non {exp[exp.is_hit==False].iface_helix.median():.2f}",
      "MWU p=%.3f"%stats.mannwhitneyu(exp[exp.is_hit==True].iface_helix,exp[exp.is_hit==False].iface_helix)[0:2][1])

# ---- L2 contact-pair chemistry: binder-class x target-region, log-odds binder vs non ----
def pair_counts(ids):
    from collections import Counter
    c=Counter()
    for did in ids:
        seen=set()
        for (bt,tu,tt) in allpairs.get(did,[]):
            key=(cls(bt),region(tu));
            if (did,key) not in seen: c[key]+=1; seen.add((did,key))
    return c
exp_ids=set(exp.design_id); bind_ids=set(exp[exp.is_hit==True].design_id); non_ids=set(exp[exp.is_hit==False].design_id)
cb=pair_counts(bind_ids); cn=pair_counts(non_ids); nb=len(bind_ids); nn=len(non_ids)
rows=[]
for key in sorted(set(cb)|set(cn)):
    kb=cb.get(key,0); kn=cn.get(key,0)
    OR,p=stats.fisher_exact([[kb,nb-kb],[kn,nn-kn]])
    rows.append({"binder_class":key[0],"target_region":key[1],"binders_with":kb,"nonbinders_with":kn,
                 "frac_binder":round(kb/nb,2),"frac_non":round(kn/nn,2),"odds_ratio":round(OR,2),"p":round(p,4)})
l2=pd.DataFrame(rows).sort_values(["p","binder_class","target_region"],kind="stable")
l2.to_csv(RES/"motifs_L2_contactpairs.csv",index=False)
print("\n=== L2 contact-pair enrichment (binder vs non, expressed) top rows ===")
print(l2.head(12).to_string(index=False))

# ---- L3 convergent motif across teams: which (binder_class,target_region) motifs recur across >=3 teams among binders ----
print("\n=== L3 convergent contact motifs across independent teams (binders) ===")
rows=[]
for key in sorted(set(cb)):
    carriers=[did for did in sorted(bind_ids) if any((cls(bt),region(tu))==key for (bt,tu,tt) in allpairs.get(did,[]))]
    teams=set(teaminfo.loc[carriers,"team"]); cohorts=set(teaminfo.loc[carriers,"cohort"])
    kd=teaminfo.loc[carriers,"kd_arith_mean_nM_all"].dropna()
    rows.append({"motif":f"{key[0]}->{key[1]}","n_binders":len(carriers),"n_teams":len(teams),
                 "n_cohorts":len(cohorts),"median_kd":round(kd.median(),0) if len(kd) else np.nan})
l3=pd.DataFrame(rows).sort_values(["n_teams","n_binders","motif"],ascending=[False,False,True],kind="stable")
l3.to_csv(RES/"motifs_L3_convergent.csv",index=False)
print(l3.head(12).to_string(index=False))
print("\nwrote motifs_L1/L2/L3 csvs")
