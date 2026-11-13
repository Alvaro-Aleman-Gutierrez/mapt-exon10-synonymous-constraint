"""Paso 35: ¿la restricción sinónima del exón 10 se explica por composición de bases o por CpG?

La objeción. Una tasa sinónima baja no tiene por qué ser selección. Puede reflejar una tasa de
mutación local baja, y la mutación depende de la composición: los dinucleótidos CpG mutan mucho
más rápido que el resto por desaminación de la citosina metilada, y un exón pobre en CpG o rico en
GC o AT en la tercera posición cambiaría menos en sitios sinónimos sin que nadie lo restrinja.

Hasta ahora el trabajo respondía con un argumento: los exones vecinos tienen tasas normales. Este
paso lo convierte en medición.

1. Composición por exón: GC en tercera posición (GC3) y CpG, en el codón humano y como promedio
   de las 259 especies.
2. Modelos que preguntan si el exón 10 sigue siendo diferente una vez descontadas esas variables,
   junto con la restricción proteica (dN) y el número de sitios sinónimos del codón:
   - lineal sobre el log de la tasa sinónima posterior de FUBAR.
3. La prueba de permutación del paso 24, con GC3 y CpG añadidos a los estratos.
"""
import collections

import numpy as np
import pandas as pd
import statsmodels.api as sm
from Bio import SeqIO
from Bio.Data.CodonTable import standard_dna_table as TBL
from scipy import stats

RNG = np.random.default_rng(20261023)
N_PERM = 10000
REP = ["E9", "E10", "E11", "E12", "E13"]
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]

fub = pd.read_csv("results/selection/fubar_sites.tsv", sep="\t")
fub["human_pos"] = fub.human_pos.astype(int)
fub = fub.sort_values("human_pos").reset_index(drop=True)
aln = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}
hum = aln["Homo_sapiens"]
fub["codon"] = [hum[3 * (c - 1):3 * c] for c in fub.column]
cds = "".join(fub.codon)
assert len(cds) == 3 * len(fub), "el CDS humano no tiene un codón por posición"


def syn_sites(cod):
    """Sitios sinónimos de un codón, al estilo Nei y Gojobori (1986), como en el paso 24."""
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


fub["syn_sites"] = [syn_sites(c) for c in fub.codon]
fub["gc3"] = fub.codon.str[2].isin(list("GC")).astype(int)
# CpG que toca al codón: dentro de él o cruzando uno de sus bordes
cg = [i for i in range(len(cds) - 1) if cds[i:i + 2] == "CG"]
touch = np.zeros(len(fub), dtype=int)
for i in cg:
    for pos in (i, i + 1):
        touch[pos // 3] = 1
fub["cpg"] = touch
fub["log_dn"] = np.log10(fub.fel_dN + 1e-3)

# ---------------------------------------------------------------- 1. composición por exón
gc3_mam = {}
for ex, g in fub.groupby("exon"):
    vals = []
    for sp, seq in aln.items():
        cods = [seq[3 * (c - 1):3 * c] for c in g.column]
        cods = [c for c in cods if "-" not in c and len(c) == 3]
        if cods:
            vals.append(np.mean([c[2] in "GC" for c in cods]))
    gc3_mam[ex] = np.mean(vals)
comp = []
for ex in CORE:
    g = fub[fub.exon == ex]
    seq = "".join(g.codon)
    cg_n = seq.count("CG")
    exp_cg = seq.count("C") * seq.count("G") / len(seq)
    comp.append({"exon": ex, "codons": len(g), "GC3_human": g.gc3.mean(),
                 "GC3_mean_259_mammals": gc3_mam[ex], "GC_human": (seq.count("G") + seq.count("C")) / len(seq),
                 "CpG_per_100_codons": 100 * cg_n / len(g), "CpG_obs_over_exp": cg_n / exp_cg if exp_cg else np.nan,
                 "codons_touching_CpG": g.cpg.mean()})
comp = pd.DataFrame(comp)
print(comp.round(3).to_string(index=False))
comp.to_csv("results/composition_by_exon.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 2. modelos
ok = fub.dropna(subset=["syn_sites"]).copy()
ok["E10"] = (ok.exon == "E10").astype(int)
ok["log_alpha"] = np.log10(ok.fubar_alpha)
out = []


def fit_ols(formula_cols, label, data):
    X = sm.add_constant(data[formula_cols]); m = sm.OLS(data.log_alpha, X).fit(cov_type="HC3")
    b = m.params["E10"]; lo, hi = m.conf_int().loc["E10"]
    print(f"  {label}: exón 10 desplaza log10(alfa) en {b:+.2f} ({10 ** b:.2f} veces), "
          f"IC95 {10 ** lo:.2f} a {10 ** hi:.2f} veces, p = {m.pvalues['E10']:.1e}")
    out.append({"model": label, "outcome": "log10 FUBAR alpha", "effect_fold": 10 ** b,
                "ci_low": 10 ** lo, "ci_high": 10 ** hi, "p": m.pvalues["E10"], "n": int(m.nobs)})


print("\nTasa sinónima posterior de FUBAR, todos los codones")
fit_ols(["E10"], "sin ajustar", ok)
fit_ols(["E10", "log_dn", "syn_sites"], "ajustado por dN y sitios sinónimos", ok)
fit_ols(["E10", "log_dn", "syn_sites", "gc3", "cpg"], "ajustado además por GC3 y CpG", ok)
rep = ok[ok.exon.isin(REP)]
fit_ols(["E10", "log_dn", "syn_sites", "gc3", "cpg"], "solo exones 9 a 13, ajustado por todo", rep)

# La regresión logística del codón invariante con todas las covariables no converge, porque el
# conteo de invariantes (21) es pequeño frente a seis términos. La permutación de abajo, que
# empareja por los mismos factores, es la prueba del codón invariante.

# ---------------------------------------------------------------- 3. permutación con GC3 y CpG
ok["dn_bin"] = pd.cut(ok.fel_dN, bins=[-1, 1e-3, 0.05, 0.3, 10], labels=["0", "muy bajo", "bajo", "alto"])
ok["ss_bin"] = (ok.syn_sites * 2).round() / 2
ok["stratum"] = (ok.dn_bin.astype(str) + "|" + ok.ss_bin.astype(str) + "|gc3=" + ok.gc3.astype(str)
                 + "|cpg=" + ok.cpg.astype(str))
tgt = ok[ok.exon == "E10"]; pool = ok[ok.exon != "E10"]
need = collections.Counter(tgt.stratum)
by = {s: g.codon_invariant.values for s, g in pool.groupby("stratum")}
matched = sum(k for s, k in need.items() if s in by)
draws = np.empty(N_PERM, dtype=int)
for i in range(N_PERM):
    tot = 0
    for s, k in need.items():
        v = by.get(s)
        if v is None:
            continue
        tot += int(RNG.choice(v, size=k, replace=len(v) < k).sum())
    draws[i] = tot
obs = int(tgt.codon_invariant.sum())
p = (draws >= obs).mean()
print(f"\nPermutación con GC3 y CpG en los estratos: observados {obs}, codones emparejados {matched} "
      f"de {len(tgt)}, esperados {draws.mean():.2f} (p95 {np.quantile(draws, 0.95):.0f}), "
      f"p {'< ' + str(1 / N_PERM) if p == 0 else f'{p:.4f}'}")
out.append({"model": "permutación con GC3 y CpG en los estratos", "outcome": "codones invariantes",
            "effect_fold": obs / max(draws.mean(), 1e-9), "ci_low": np.nan, "ci_high": np.nan,
            "p": p if p > 0 else 1 / N_PERM, "n": matched})
pd.DataFrame(out).to_csv("results/composition_control.tsv", sep="\t", index=False)
from esm import write_esm
write_esm(comp, 11, "composition_by_exon")
