"""Paso 13: resumen de los análisis de HyPhy (paso 9).

FEL y MEME se mapean al exón humano de cada columna (misma convención que el paso 8).
En FEL, alpha = dS y beta = dN: beta < alpha con p <= 0.05 es selección purificadora, beta > alpha diversificadora.
"""
import itertools, json, os, glob
import numpy as np, pandas as pd
from Bio import AlignIO, SeqIO
from statsmodels.stats.multitest import multipletests

P = 0.05
S = "results/selection"
# FEL y MEME prueban cada columna por separado, así que se corrige por comparaciones múltiples
# con Benjamini-Hochberg (q <= 0.05) además del p sin corregir.
Q = 0.05
aln = AlignIO.read("results/mapt_2N4R_final.aln.faa", "fasta")
hx = pd.read_csv("data/human_exons.tsv", sep="\t")
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
lens = dict(zip(hx.exon, hx.aa))
bounds = list(itertools.accumulate(lens[e] for e in CORE))
hum = str(next(r.seq for r in aln if r.id == "Homo_sapiens"))
col_exon, col_pos, k = [], [], 0
for c in hum:
    if c == "-":
        col_exon.append(col_exon[-1] if col_exon else "E1"); col_pos.append(None); continue
    k += 1
    col_exon.append(CORE[next(i for i, b in enumerate(bounds) if k <= b)]); col_pos.append(k)

def mle(js):
    j = json.load(open(js))
    h = [x[0] for x in j["MLE"]["headers"]]
    return pd.DataFrame(j["MLE"]["content"]["0"], columns=h)

sites = pd.DataFrame({"column": range(1, len(hum) + 1), "exon": col_exon, "human_pos": col_pos})
fel = mle(f"{S}/all.FEL.json")
sites["fel_dS"], sites["fel_dN"], sites["fel_p"] = fel["alpha"], fel["beta"], fel["p-value"]
sites["fel_q"] = multipletests(sites.fel_p, alpha=Q, method="fdr_bh")[1]
sites["fel"] = np.where(sites.fel_p > P, "neutral", np.where(sites.fel_dN < sites.fel_dS, "purifying", "diversifying"))
if os.path.exists(f"{S}/all.MEME.json") and os.path.getsize(f"{S}/all.MEME.json"):
    meme = mle(f"{S}/all.MEME.json")
    sites["meme_p"] = meme["p-value"]
    sites["meme_q"] = multipletests(sites.meme_p, alpha=Q, method="fdr_bh")[1]
    sites["meme"] = sites.meme_p <= P
# codones idénticos en todas las especies (ignorando huecos): FEL no tiene poder en ellos
cds = [str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")]
sites["codon_invariant"] = [len({x[3*i:3*i+3] for x in cds} - {"---"}) == 1 for i in range(len(sites))]
sites.to_csv(f"{S}/sites_summary.tsv", sep="\t", index=False)

# solo columnas con residuo humano; las 122 columnas de huecos en humano se omiten
by_exon = sites[sites.human_pos.notna()].groupby("exon", sort=False).agg(n_columns=("column", "size"),
            codon_invariant=("codon_invariant", "sum"), median_dS=("fel_dS", "median"),
            purifying_q=("fel_q", lambda x: (x <= Q).sum()),
            purifying=("fel", lambda x: (x == "purifying").sum()),
            diversifying=("fel", lambda x: (x == "diversifying").sum()))
by_exon["frac_purifying"] = (by_exon.purifying / by_exon.n_columns).round(2)
if "meme" in sites:
    by_exon["meme"] = sites[sites.human_pos.notna()].groupby("exon", sort=False).meme.sum()
by_exon.to_csv(f"{S}/fel_by_exon.tsv", sep="\t")
print(by_exon.to_string())

rows = []
for js in sorted(glob.glob(f"{S}/*.SLAC.json")):
    j = json.load(open(js))
    w = j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"][0][0]
    rows.append((os.path.basename(js).split(".")[0], j["input"]["number of sequences"], round(w, 3)))
slac = pd.DataFrame(rows, columns=["order", "n", "omega"]); slac.to_csv(f"{S}/slac_omega_by_order.tsv", sep="\t", index=False)
print(slac.to_string(index=False))

rows = []
for tag in ("primates", "catarrhini", "hominoidea"):
    js = f"{S}/{tag}.RELAX.json"
    if os.path.exists(js) and os.path.getsize(js):
        r = json.load(open(js))["test results"]
        rows.append((tag, round(r["relaxation or intensification parameter"], 3), round(r["LRT"], 2), r["p-value"]))
relax = pd.DataFrame(rows, columns=["test_clade", "k", "LRT", "p"]); relax.to_csv(f"{S}/relax_summary.tsv", sep="\t", index=False)
print(relax.to_string(index=False))

if "meme" in sites:
    h = sites[sites.human_pos.notna()]
    print(f"\nMEME: {(h.meme_p <= P).sum()} sitios con p <= {P} de {len(h)} "
          f"(esperados por azar {P * len(h):.0f}), {(h.meme_q <= Q).sum()} con q <= {Q}")
    sig = h[h.meme_q <= Q][["column", "exon", "human_pos", "meme_p", "meme_q", "fel_dS", "fel_dN"]]
    print(sig.to_string(index=False))
