"""Paso 21: divergencia de los exones alternativos, para justificar su exclusión.

Los cuatro exones alternativos (4aL, 4a, 6, 8) no entran en el análisis de selección porque el
conjunto de datos está definido como 2N4R, una isoforma nombrada que todas las especies comparten.
Este paso mide cuánto difieren igual, de modo que la exclusión quede sustentada con números y no
solo afirmada. También deja ver que E6 está más conservado que E1 y E5, que sí entran, lo que
confirma que el criterio de exclusión es la composición de isoforma y no la calidad de secuencia.

Para cada exón y especie se toma la isoforma anotada con mayor cobertura de ese exón y se alinea a
la super-isoforma humana de 880 residuos con los mismos parámetros del paso 4.
"""
import collections
import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

ALT = ["E4aL", "E4a", "E6", "E8"]
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]

hx = pd.read_csv("data/human_exons.tsv", sep="\t")
bounds, k = {}, 0
for e in hx.exon:                         # super-isoforma en orden genómico
    aa = int(hx.loc[hx.exon == e, "aa"].iloc[0]); bounds[e] = (k, k + aa); k += aa
SUP = "".join(hx.peptide)
assert len(SUP) == k

al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                           substitution_matrix=substitution_matrices.load("BLOSUM62"))
al.target_end_gap_score = al.query_end_gap_score = 0.0

sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
order_of = dict(zip(sel.species, sel["order"]))
pres = pd.read_csv("results/exon_presence_all_isoforms.tsv", sep="\t")
pres = pres[pres.species.isin(order_of)]
prot = {r.id.split()[0]: str(r.seq) for r in SeqIO.parse("data/mammal_protein.faa", "fasta")}

def mapping(aa):
    m = {}
    for (a0, a1), (b0, b1) in zip(*al.align(SUP, aa)[0].aligned):
        for off in range(a1 - a0):
            m[a0 + off] = b0 + off
    return m

rows = []
for e in ALT:
    s, t = bounds[e]
    best = pres.loc[pres.groupby("species")[e].idxmax()]
    for _, r in best[best[e] >= 0.5].iterrows():
        aa = prot.get(r.protein)
        if not aa:
            continue
        m = mapping(aa)
        pos = [i for i in range(s, t) if i in m]
        if len(pos) < 0.5 * (t - s):
            continue
        rows.append({"exon": e, "species": r.species, "order": order_of[r.species],
                     "protein": r.protein, "aligned_aa": len(pos),
                     "identity_to_human": round(sum(1 for i in pos if SUP[i] == aa[m[i]]) / len(pos), 4)})

# los exones núcleo, medidos igual, como referencia
for e in CORE:
    s, t = bounds[e]
    for sp, p in zip(sel.species, sel.protein):
        aa = prot.get(p)
        if not aa:
            continue
        m = mapping(aa)
        pos = [i for i in range(s, t) if i in m]
        if len(pos) < 0.9 * (t - s):
            continue
        rows.append({"exon": e, "species": sp, "order": order_of[sp], "protein": p,
                     "aligned_aa": len(pos),
                     "identity_to_human": round(sum(1 for i in pos if SUP[i] == aa[m[i]]) / len(pos), 4)})

df = pd.DataFrame(rows)
from esm import write_esm
write_esm(df, 6, "alt_exon_divergence")

summ = (df.groupby("exon")
          .agg(set_=("exon", lambda x: "alternative" if x.iloc[0] in ALT else "core"),
               species=("species", "nunique"),
               mean_identity=("identity_to_human", "mean"),
               min_identity=("identity_to_human", "min"),
               min_aligned_aa=("aligned_aa", "min"),
               max_aligned_aa=("aligned_aa", "max"))
          .rename(columns={"set_": "set"}))
summ["human_aa"] = [bounds[e][1] - bounds[e][0] for e in summ.index]
summ["species_pct"] = (summ.species / len(sel) * 100).round(0).astype(int)
summ = summ.loc[[e for e in ALT if e in summ.index] + [e for e in CORE if e in summ.index]]
summ = summ.round({"mean_identity": 3, "min_identity": 3})
summ.reset_index().to_csv("results/tables/tableS2_alt_exon_summary.tsv", sep="\t", index=False)

cap = ("Table S2. Divergence of the four alternative exons of MAPT, which the selection analyses do "
       "not cover, against the eleven core exons that they do. For each species the annotated "
       "isoform with the highest coverage of the exon was aligned to the human 880 residue super "
       "isoform. Identity is computed over positions where both sequences carry a residue, so a "
       "partially recovered exon still contributes.")
with open("manuscript/suplementarios/tables/tableS2_alt_exon_summary.md", "w") as f:
    f.write(f"**{cap}**\n\n")
    out = summ.reset_index()[["exon", "set", "human_aa", "species", "species_pct",
                              "mean_identity", "min_identity", "min_aligned_aa", "max_aligned_aa"]]
    out.columns = ["Exon", "Set", "Residues in human", "Species", "% of 259",
                   "Mean identity", "Minimum", "Shortest aligned", "Longest aligned"]
    f.write(out.to_markdown(index=False, floatfmt=".3f") + "\n")
print(summ.to_string())
print("\nfilas en el recurso en línea 6:", len(df))
