"""Paso 8: conservación por exón (identidad con humano y entre pares dentro de cada orden)
y presencia de exones alternativos (4aL, 4a, 6, 8) en las isoformas anotadas, por orden."""
import itertools
import numpy as np, pandas as pd
from Bio import AlignIO

aln = AlignIO.read("results/mapt_2N4R_final.aln.faa", "fasta")
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
order_of = dict(zip(sel.species.str.replace(" ", "_"), sel.order))
hx = pd.read_csv("data/human_exons.tsv", sep="\t")
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]

# columna del alineamiento -> exón humano
hum = str(next(r.seq for r in aln if r.id == "Homo_sapiens"))
lens = dict(zip(hx.exon, hx.aa)); col_exon, k = [], 0
bounds = list(itertools.accumulate(lens[e] for e in CORE)); 
for c in hum:
    if c == "-":
        col_exon.append(col_exon[-1] if col_exon else "E1"); continue
    col_exon.append(CORE[next(i for i, b in enumerate(bounds) if k < b)]); k += 1
col_exon = np.array(col_exon)
M = np.array([list(str(r.seq)) for r in aln]); ids = [r.id for r in aln]
H = M[ids.index("Homo_sapiens")]

def ident(a, b, mask):
    """Identidad sobre columnas donde ambas tienen residuo (gaps no cuentan)."""
    m = mask & (a != "-") & (b != "-")
    return (a[m] == b[m]).mean() if m.any() else np.nan

rows = []
for i, sp in enumerate(ids):
    for e in CORE:
        rows.append([sp, order_of[sp], e, ident(M[i], H, col_exon == e)])
    rows.append([sp, order_of[sp], "total", ident(M[i], H, np.ones(len(H), bool))])
idh = pd.DataFrame(rows, columns=["species", "order", "exon", "identity_to_human"])
idh.to_csv("results/exon_identity_to_human.tsv", sep="\t", index=False)

# identidad media entre pares dentro de cada orden (n >= 3)
rows = []
for o, g in pd.Series(order_of).groupby(lambda s: order_of[s]):
    idx = [ids.index(s) for s in g.index]
    if len(idx) < 3:
        continue
    for e in CORE + ["total"]:
        mask = np.ones(len(H), bool) if e == "total" else col_exon == e
        v = [ident(M[a], M[b], mask) for a, b in itertools.combinations(idx, 2)]
        rows.append([o, len(idx), e, np.nanmean(v)])
wo = pd.DataFrame(rows, columns=["order", "n", "exon", "mean_pairwise_identity"])
wo.to_csv("results/exon_identity_within_order.tsv", sep="\t", index=False)

# presencia de exones alternativos en alguna isoforma anotada de cada especie
pres = pd.read_csv("results/exon_presence_all_isoforms.tsv", sep="\t")
alt = pres.groupby("species")[["E4aL", "E4a", "E6", "E8"]].max().ge(0.5)
alt["order"] = alt.index.map(dict(zip(sel.species, sel.order)))
alt = alt.dropna(subset=["order"])
alt.groupby("order").agg(["mean"]).round(2).assign(n=alt.groupby("order").size()) \
   .to_csv("results/alt_exon_presence_by_order.tsv", sep="\t")

pd.set_option("display.width", 200)
print(idh.pivot_table(index="order", columns="exon", values="identity_to_human", aggfunc="mean")[CORE + ["total"]].round(2))
print(wo.pivot(index="order", columns="exon", values="mean_pairwise_identity")[CORE + ["total"]].round(2))
print(pd.read_csv("results/alt_exon_presence_by_order.tsv", sep="\t", header=[0, 1], index_col=0))
