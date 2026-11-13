"""Paso 25: ¿los codones invariantes del exón 10 coinciden con potenciadores de splicing predichos?

La discusión reconoce que el trabajo muestra que ciertos codones no cambian y que hay elementos
reguladores publicados en el mismo exón, sin conectar las dos cosas. Esta prueba no reemplaza un
experimento, pero convierte la coincidencia en una prueba estadística.

Conjunto de hexámeros, RESCUE-ESE, 238 secuencias de Fairbrother et al. 2002, Science 297:1007.
Se recuperó del archivo del servidor del laboratorio Burge y está en data/ese/.

La nula es que los codones invariantes se ubican al azar dentro del exón. Se permutan codones
completos, no nucleótidos sueltos, porque la invariancia se define por codón.

Advertencia que va al manuscrito: RESCUE-ESE se derivó de exones constitutivos humanos, así que
predice potenciadores genéricos y no los específicos de MAPT descritos en minigenes.
"""
import numpy as np
import pandas as pd
from Bio import SeqIO

RNG = np.random.default_rng(20261022)
N_PERM = 10000

hexamers = {l.strip() for l in open("data/ese/RESCUE-ESE_hexamers.txt") if l.strip()}
print(f"hexámeros RESCUE-ESE cargados: {len(hexamers)}")

sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()].copy()
sites["human_pos"] = sites.human_pos.astype(int)
aln = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}
hum = aln["Homo_sapiens"]

def exon_profile(exon):
    """Secuencia del exón y qué codones son invariantes, en orden."""
    sub = sites[sites.exon == exon].sort_values("human_pos")
    seq = "".join(hum[3 * (c - 1):3 * c] for c in sub.column)
    return seq, sub.codon_invariant.values.astype(bool)

def covered(seq):
    """Nucleótidos cubiertos por al menos un hexámero del conjunto."""
    cov = np.zeros(len(seq), dtype=bool)
    for i in range(len(seq) - 5):
        if seq[i:i + 6] in hexamers:
            cov[i:i + 6] = True
    return cov

def overlap(inv_mask, cov):
    """Nucleótidos que son a la vez invariantes y están cubiertos."""
    nt = np.repeat(inv_mask, 3)
    return int((nt & cov).sum())

rows = []
for exon in ["E10", "E9", "E11", "E12", "E13"]:
    seq, inv = exon_profile(exon)
    if "-" in seq or inv.sum() == 0:
        print(f"{exon}: sin codones invariantes, se omite")
        continue
    cov = covered(seq)
    obs = overlap(inv, cov)
    k = int(inv.sum())
    draws = np.empty(N_PERM, dtype=int)
    idx = np.arange(len(inv))
    for i in range(N_PERM):
        m = np.zeros(len(inv), dtype=bool)
        m[RNG.choice(idx, size=k, replace=False)] = True
        draws[i] = overlap(m, cov)
    p = (draws >= obs).mean()
    rows.append({"exon": exon, "codons": len(inv), "invariant": k,
                 "nt_covered_by_ESE": int(cov.sum()), "nt_total": len(seq),
                 "coverage": round(cov.mean(), 3),
                 "observed_overlap_nt": obs, "expected_mean": round(draws.mean(), 1),
                 "p_empirical": round(float(p), 4)})
    print(f"{exon}: {cov.sum()}/{len(seq)} nt cubiertos ({cov.mean():.0%}); "
          f"solapamiento observado {obs} nt, esperado {draws.mean():.1f}, p = {p:.4f}")

res = pd.DataFrame(rows)
res.to_csv("results/tables/table5_ese_overlap.tsv", sep="\t", index=False)
cap = ("Table 5. Overlap between nucleotide invariant codons and predicted exonic splicing "
       "enhancers. Enhancers are the 238 RESCUE-ESE hexamers of Fairbrother et al. (2002), matched "
       "by sliding a six nucleotide window along each exon. The null distribution comes from 10,000 "
       "permutations that relocate the invariant codons within the same exon, keeping their number "
       "and the enhancer map fixed. RESCUE-ESE was derived from constitutively spliced human exons, "
       "so it predicts generic enhancers rather than the elements mapped in MAPT minigenes.")
with open("manuscript/tables/table5_ese_overlap.md", "w") as f:
    o = res.copy()
    o.columns = ["Exon", "Codons", "Invariant", "Nucleotides in an ESE", "Nucleotides total",
                 "ESE coverage", "Overlap observed (nt)", "Overlap expected (nt)", "p"]
    f.write(f"**{cap}**\n\n" + o.to_markdown(index=False, floatfmt=".3f") + "\n")
print("\ntabla 9 escrita")
