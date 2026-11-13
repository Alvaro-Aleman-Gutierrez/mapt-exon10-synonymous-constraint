"""Paso 24: ¿la invariancia del exón 10 se explica por azar? Prueba de permutación.

El resultado central del paper hoy descansa en comparar conteos, o sea 14 codones invariantes de 31
en el exón 10 contra 2 de 28 en el exón 11. Esto lo vuelve una prueba.

Dos confundidores hay que controlar, porque sin ellos el valor p sale significativo trivialmente.

El primero es la restricción sobre la proteína. Un codón cuyo aminoácido nunca cambia tiene más
chance de conservar también sus nucleótidos, y el exón 10 tiene dN casi nulo en casi todas sus
posiciones. Se controla emparejando por dN.

El segundo es la degeneración del codón. Un codón de Met o Trp no tiene sitios sinónimos, así que
si su aminoácido está fijo queda invariante por construcción. Un codón cuádruple degenerado tiene
tres formas de variar sinónimamente. Se controla contando los sitios sinónimos de cada codón
humano al estilo Nei-Gojobori y emparejando por ese número.

La hipótesis nula es que, entre codones con la misma restricción proteica y la misma degeneración,
la invariancia se distribuye al azar a lo largo del gen.
"""
import collections
import numpy as np
import pandas as pd
from Bio.Data.CodonTable import standard_dna_table as TBL
from Bio import SeqIO

RNG = np.random.default_rng(20261022)
N_PERM = 10000
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
TEST = ["E10", "E9", "E11", "E12", "E13"]

sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()].copy()
sites["human_pos"] = sites.human_pos.astype(int)

# codón humano de cada posición, desde el alineamiento de codones
aln = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}
hum = aln["Homo_sapiens"]
sites["codon"] = [hum[3 * (c - 1):3 * c] for c in sites.column]

def syn_sites(cod):
    """Sitios sinónimos de un codón, al estilo Nei-Gojobori.

    Cada una de las tres posiciones aporta la fracción de sus tres cambios posibles que no alteran
    el aminoácido. Met y Trp dan 0, un codón cuádruple degenerado da 1 en la tercera posición.
    """
    if cod in TBL.stop_codons or "-" in cod or "N" in cod:
        return np.nan
    aa = TBL.forward_table.get(cod)
    if aa is None:
        return np.nan
    tot = 0.0
    for i in range(3):
        same = 0
        for b in "ACGT":
            if b == cod[i]:
                continue
            alt = cod[:i] + b + cod[i + 1:]
            if alt not in TBL.stop_codons and TBL.forward_table.get(alt) == aa:
                same += 1
        tot += same / 3
    return tot

sites["syn_sites"] = [syn_sites(c) for c in sites.codon]
ok = sites.dropna(subset=["syn_sites"]).copy()
print(f"codones utilizables {len(ok)} de {len(sites)}; invariantes {int(ok.codon_invariant.sum())}")
print("\nsitios sinónimos por exón, media")
print(ok.groupby("exon").syn_sites.mean().round(2).reindex(CORE).to_string())

# estratos: cuartil de dN por rango completo, cruzado con sitios sinónimos redondeados a 0.5
ok["dn_bin"] = pd.cut(ok.fel_dN, bins=[-1, 1e-3, 0.05, 0.3, 10], labels=["0", "muy bajo", "bajo", "alto"])
ok["ss_bin"] = (ok.syn_sites * 2).round() / 2
ok["stratum"] = ok.dn_bin.astype(str) + "|" + ok.ss_bin.astype(str)

rows = []
for target in TEST:
    tgt = ok[ok.exon == target]
    obs = int(tgt.codon_invariant.sum())
    pool = ok[ok.exon != target]
    # un sorteo replica la composición de estratos del exón evaluado
    need = collections.Counter(tgt.stratum)
    by = {s: g.codon_invariant.values for s, g in pool.groupby("stratum")}
    usable = sum(k for s, k in need.items() if s in by and len(by[s]) >= 1)
    draws = np.empty(N_PERM, dtype=int)
    for i in range(N_PERM):
        tot = 0
        for s, k in need.items():
            v = by.get(s)
            if v is None or len(v) == 0:
                continue
            tot += int(RNG.choice(v, size=k, replace=len(v) < k).sum())
        draws[i] = tot
    p = (draws >= obs).mean()
    rows.append({"exon": target, "codons": len(tgt), "observed_invariant": obs,
                 "matched_codons": usable, "expected_mean": round(draws.mean(), 2),
                 "expected_p95": int(np.quantile(draws, 0.95)),
                 "p_empirical": p if p > 0 else f"<{1/N_PERM:g}"})
    print(f"\n{target}: observado {obs} de {len(tgt)} | esperado {draws.mean():.2f} "
          f"(p95 {np.quantile(draws,0.95):.0f}) | p {'<'+str(1/N_PERM) if p==0 else f'{p:.4f}'}")

res = pd.DataFrame(rows)
res.to_csv("results/tables/table3_permutation.tsv", sep="\t", index=False)

cap = ("Table 3. Permutation test of nucleotide invariance in exons 9 to 13, which span the microtubule binding region and the carboxy terminal region. "
       "For each exon, 10,000 replicate sets of codons were drawn from the rest of the gene, matched "
       "on the per site nonsynonymous rate and on the number of synonymous sites in the human codon, "
       "so that neither protein level constraint nor codon degeneracy can produce the result. The "
       "empirical p value is the fraction of replicates reaching the observed count or more.")
with open("manuscript/tables/table3_permutation.md", "w") as f:
    o = res.copy()
    o.columns = ["Exon", "Codons", "Invariant observed", "Matched codons available",
                 "Expected mean", "Expected 95th percentile", "p"]
    f.write(f"**{cap}**\n\n" + o.to_markdown(index=False, floatfmt=".2f") + "\n")

# comprobación exacta por hipergeométrica sobre el subconjunto de dN casi nulo
from scipy.stats import hypergeom
low = ok[ok.fel_dN < 1e-3]
K, N = int(low.codon_invariant.sum()), len(low)
k10, n10 = int(low[low.exon == "E10"].codon_invariant.sum()), len(low[low.exon == "E10"])
ph = hypergeom.sf(k10 - 1, N, K, n10)
print(f"\ncomprobacion hipergeometrica sobre dN<1e-3: {k10} de {n10} en E10, "
      f"{K} de {N} en el gen, p = {ph:.3e}")
