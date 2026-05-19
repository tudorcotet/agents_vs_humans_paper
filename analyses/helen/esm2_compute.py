"""
Compute ESM-2 representations for all 141 designs.

Model: esm2_t33_650M_UR50D — the 650M-parameter ESM-2, the de facto default
for embeddings/scoring. Outputs (cached, idempotent):

- analyses/helen/esm2_cache.npz
    design_id : (141,) int
    emb       : (141, 1280) float32  — mean-pooled final-layer representation
    pll       : (141,) float32       — length-normalized single-pass
                                       log-likelihood (pseudo-log-likelihood
                                       approximation; true masked-marginal is
                                       infeasible on CPU for 650M)

Re-run is a no-op if the cache exists and covers all design_ids.
"""
from __future__ import annotations

import numpy as np
import torch

from scripts.utils import load_designs, repo_root

CACHE = repo_root() / "analyses" / "helen" / "esm2_cache.npz"
MODEL_NAME = "esm2_t33_650M_UR50D"
REPR_LAYER = 33
BATCH = 4


def _device() -> str:
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def main() -> None:
    df = load_designs()[["design_id", "sequence"]].copy()
    df = df.sort_values("design_id").reset_index(drop=True)

    if CACHE.exists():
        cached = np.load(CACHE)
        if set(cached["design_id"].tolist()) == set(df.design_id.tolist()):
            print(f"[esm2] cache hit ({CACHE.name}), {len(cached['design_id'])} designs — skipping")
            return

    import esm  # imported lazily so the figure script doesn't hard-depend on it

    device = _device()
    print(f"[esm2] loading {MODEL_NAME} on {device} …")
    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    ids = df.design_id.to_numpy()
    seqs = df.sequence.tolist()
    embs: list[np.ndarray] = []
    plls: list[float] = []

    for start in range(0, len(seqs), BATCH):
        chunk = seqs[start:start + BATCH]
        data = [(str(i), s) for i, s in enumerate(chunk)]
        _, _, toks = batch_converter(data)
        toks = toks.to(device)
        with torch.no_grad():
            out = model(toks, repr_layers=[REPR_LAYER], return_contacts=False)
        reps = out["representations"][REPR_LAYER]
        logits = out["logits"]
        log_probs = torch.log_softmax(logits, dim=-1)
        for b, seq in enumerate(chunk):
            L = len(seq)
            # tokens: [BOS] residues... [EOS] — residue span is 1 .. L
            res = reps[b, 1:L + 1].mean(0).float().cpu().numpy()
            embs.append(res)
            tok_ids = toks[b, 1:L + 1]
            lp = log_probs[b, 1:L + 1].gather(-1, tok_ids.unsqueeze(-1)).squeeze(-1)
            plls.append(float(lp.mean().cpu()))
        print(f"[esm2] {min(start + BATCH, len(seqs))}/{len(seqs)} sequences")

    np.savez(
        CACHE,
        design_id=ids.astype(int),
        emb=np.vstack(embs).astype(np.float32),
        pll=np.asarray(plls, dtype=np.float32),
    )
    print(f"[esm2] wrote {CACHE} — emb {np.vstack(embs).shape}, pll n={len(plls)}")


if __name__ == "__main__":
    main()
