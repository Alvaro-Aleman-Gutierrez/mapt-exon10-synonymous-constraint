"""Paso 34: contexto genómico del exón 10 de MAPT.

Pregunta. ¿Cuán extremo es el exón 10 entre los exones codificantes humanos? Hasta aquí el trabajo
mostró que su restricción sinónima supera la de sus exones vecinos, pero no si eso es raro o común
en el genoma. Sin esa referencia, "extremo" no se puede cuantificar.

Datos. Los elementos de restricción sinónima (SCE) y la tasa sinónima local (lambda_s) de Lin et al.
(2011), que cubren todo el genoma en coordenadas hg18, y los exones codificantes de RefSeq en hg18.
No se vuelven a estimar tasas: se usan las de Lin et al. tal como se publicaron, de modo que esto
es una comparación de posición dentro de su distribución y no un análisis independiente.

Métricas por exón.
- cobertura: fracción de los nucleótidos del exón que cubre un SCE, a resolución de 9, 15 y 30
  codones. A resolución de 30 codones un elemento exige restricción sostenida sobre unos 90 nt.
- lambda_s medio: media de las ventanas de tasa sinónima local que tocan el exón. Valores bajos
  significan menos cambio sinónimo del esperado para el marco de lectura de su gen.

Comparación. El exón 10 mide 93 nt. Se compara con los exones de largo parecido, porque la
cobertura de un exón corto o muy largo no es comparable con la de uno de 93 nt, y por separado con
los alternativos, definidos como los que faltan en al menos un transcrito NM_ del mismo gen.
"""
import gzip
import os
from collections import defaultdict

import numpy as np
import pandas as pd

D = "data/lin2011_genomewide"
OUT = "results"

# ---------------------------------------------------------------- exones codificantes de RefSeq
exons = defaultdict(set)         # (cromosoma, inicio, fin) -> {transcritos}
gene_of = {}
tx_of_gene = defaultdict(set)
with gzip.open(f"{D}/refGene.hg18.txt.gz", "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        name, chrom = p[1], p[2]
        if not name.startswith("NM_") or "_" in chrom[3:]:
            continue
        cds0, cds1 = int(p[6]), int(p[7])
        starts = [int(x) for x in p[9].strip(",").split(",")]
        ends = [int(x) for x in p[10].strip(",").split(",")]
        gene = p[12]
        tx_of_gene[gene].add(name)
        for s, e in zip(starts, ends):
            a, b = max(s, cds0), min(e, cds1)
            if b > a:
                exons[(chrom, a, b)].add(name)
                gene_of[(chrom, a, b)] = gene

rows = []
for (chrom, a, b), tx in exons.items():
    gene = gene_of[(chrom, a, b)]
    nt = len(tx_of_gene[gene])
    kind = "unclassified" if nt < 2 else ("constitutive" if len(tx) == nt else "alternative")
    rows.append((chrom, a, b, gene, b - a, kind))
ex = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "length", "kind"])
print(f"exones codificantes únicos: {len(ex)}  genes: {ex.gene.nunique()}")
print(ex.kind.value_counts().to_string(), "\n")

# ---------------------------------------------------------------- cobertura por SCE
def merged_blocks(path):
    """Bloques de un BED12 por cromosoma, fusionados, con la suma acumulada de largo."""
    per = defaultdict(list)
    for line in open(path):
        if line.startswith("track"):
            continue
        p = line.split("\t")
        s = int(p[1])
        sizes = [int(x) for x in p[10].strip(",").split(",")]
        offs = [int(x) for x in p[11].strip().strip(",").split(",")]
        for o, z in zip(offs, sizes):
            per[p[0]].append((s + o, s + o + z))
    out = {}
    for c, iv in per.items():
        iv.sort()
        m = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        st = np.array([x[0] for x in m]); en = np.array([x[1] for x in m])
        out[c] = (st, en, np.concatenate([[0], np.cumsum(en - st)]))
    return out


def covered_before(x, st, en, cum):
    """Nucleótidos cubiertos con coordenada menor que x."""
    i = np.searchsorted(st, x, side="right")      # intervalos que empiezan antes de x
    if i == 0:
        return 0
    return cum[i - 1] + min(x, en[i - 1]) - st[i - 1]


for res in (9, 15, 30):
    blocks = merged_blocks(f"{D}/SCE{res}.hg18.bed")
    cov = np.zeros(len(ex))
    for k, (c, a, b) in enumerate(zip(ex.chrom, ex.start, ex.end)):
        if c in blocks:
            st, en, cum = blocks[c]
            cov[k] = (covered_before(b, st, en, cum) - covered_before(a, st, en, cum)) / (b - a)
    ex[f"cov{res}"] = cov
    print(f"SCE{res}: exones con cobertura > 0: {int((cov > 0).sum())}")

# ---------------------------------------------------------------- lambda_s medio por exón
bg = pd.read_csv(f"{D}/lambda_s_ORF.hg18.bedGraph", sep="\t", header=None, skiprows=1,
                 names=["chrom", "s", "e", "v"])
lam = np.full(len(ex), np.nan)
nwin = np.zeros(len(ex), dtype=int)
pref = {}
for c, g in bg.groupby("chrom"):
    g = g.sort_values("s")
    pref[c] = (g.s.values, g.e.values, g.v.values)
for k, (c, a, b) in enumerate(zip(ex.chrom, ex.start, ex.end)):
    if c not in pref:
        continue
    s, e, v = pref[c]
    hi = np.searchsorted(s, b, side="left")
    lo = np.searchsorted(s, a - 60, side="left")     # las ventanas son cortas, 60 nt es margen
    idx = np.arange(lo, hi)
    idx = idx[e[idx] > a]
    if len(idx):
        lam[k] = v[idx].mean()
        nwin[k] = len(idx)
ex["lambda_mean"] = lam
ex["n_windows"] = nwin
ex.to_csv(f"{OUT}/genomewide_exons.tsv.gz", sep="\t", index=False, compression="gzip")

# ---------------------------------------------------------------- MAPT
mapt = ex[ex.gene == "MAPT"].sort_values("start")
print("\nMAPT, exones codificantes de RefSeq")
print(mapt[["start", "end", "length", "kind", "cov9", "cov15", "cov30", "lambda_mean"]]
      .round(3).to_string(index=False))
e10 = mapt[(mapt.start == 41443512) & (mapt.end == 41443605)].iloc[0]
e6 = mapt[mapt.length == 198].iloc[0]

# ---------------------------------------------------------------- posición de E10
def report(label, e, lo, hi):
    print(f"\n=== {label}: {int(e.length)} nt  cov9 {e.cov9:.2f}  cov15 {e.cov15:.2f}  "
          f"cov30 {e.cov30:.2f}  lambda medio {e.lambda_mean:.3f}")
    strata = {"todos los exones": ex,
              f"largo parecido ({lo} a {hi} nt)": ex[(ex.length >= lo) & (ex.length <= hi)],
              "largo parecido y alternativos": ex[(ex.length >= lo) & (ex.length <= hi) &
                                                  (ex.kind == "alternative")]}
    res = []
    for name, s in strata.items():
        s = s[s.n_windows > 0]
        r = {"comparación": name, "exones": len(s)}
        for m in ("cov15", "cov30"):
            r[f"{m} >= exón"] = int((s[m] >= e[m] - 1e-9).sum())
        r["lambda <= exón"] = int((s.lambda_mean <= e.lambda_mean + 1e-9).sum())
        res.append(r)
        print(f"  {name}: {len(s)} exones | con cobertura 30 >= la del exón: {r['cov30 >= exón']} "
              f"({r['cov30 >= exón'] / len(s):.4f}) | con cobertura 15 >= : {r['cov15 >= exón']} "
              f"({r['cov15 >= exón'] / len(s):.4f}) | con lambda medio <= : {r['lambda <= exón']} "
              f"({r['lambda <= exón'] / len(s):.4f})")
    return pd.DataFrame(res)


a = report("MAPT exón 10", e10, 75, 115)
b = report("MAPT exón 6", e6, 150, 250)
a.assign(exon="E10").to_csv(f"{OUT}/genomewide_context.tsv", sep="\t", index=False)
b.assign(exon="E6").to_csv(f"{OUT}/genomewide_context.tsv", sep="\t", index=False, mode="a",
                           header=False)

# ---------------------------------------------------------------- figura suplementaria (ESM_10)
import matplotlib.pyplot as plt
from figstyle import FOCUS, REF, MARK, INK, MUTED, MM, W2, panel, despine, save_esm

s = ex[(ex.length >= 75) & (ex.length <= 115) & (ex.n_windows > 0)]
cls = {"alternative": ("alternatively spliced", REF, "-"),
       "constitutive": ("constitutive", MUTED, "--")}
fig, axs = plt.subplots(1, 2, figsize=(W2, 82 * MM))
fig.subplots_adjust(left=0.085, right=0.985, top=0.93, bottom=0.17, wspace=0.34)

ax = axs[0]
grid = np.linspace(0.02, 1.0, 200)
for k, (lab, col, ls) in cls.items():
    v = s[s.kind == k].cov30.values
    ax.plot(grid, [(v >= g).mean() for g in grid], color=col, linestyle=ls, linewidth=1.3)
alt = s[s.kind == "alternative"].cov30.values
for e, col, name, dy in ((e10, FOCUS, "exon 10", 1.35),):
    ax.scatter([e.cov30], [(alt >= e.cov30 - 1e-9).mean()], s=34, color=col, zorder=5,
               edgecolor="white", linewidth=0.6)
    ax.annotate(name, (e.cov30, (alt >= e.cov30 - 1e-9).mean() * dy), ha="right", va="center",
                color=col, fontsize=8, xytext=(-6, 0), textcoords="offset points")
ax.set_yscale("log")
ax.set_yticks([0.01, 0.02, 0.05, 0.1, 0.2])
ax.set_yticklabels(["1 %", "2 %", "5 %", "10 %", "20 %"])
ax.minorticks_off()
ax.set_ylim(0.008, 0.3)
ax.set_xlim(0, 1.03)
ax.set_xlabel("Fraction of the exon covered by 30 codon elements")
ax.set_ylabel("Fraction of exons covered at least this much")
ax.annotate("alternatively spliced", (0.04, 0.205), color=REF, fontsize=8)
ax.annotate("constitutive", (0.04, 0.052), color=MUTED, fontsize=8)
despine(ax)
panel(ax, "a", dx=-0.22, dy=1.02)

ax = axs[1]
gx = np.linspace(0, 1.6, 300)
for k, (lab, col, ls) in cls.items():
    v = np.sort(s[s.kind == k].lambda_mean.values)
    ax.plot(gx, np.searchsorted(v, gx, side="right") / len(v), color=col, linestyle=ls, linewidth=1.3)
va = np.sort(s[s.kind == "alternative"].lambda_mean.values)
for e, col, name, dx_ in ((e10, FOCUS, "exon 10", 6),):
    y = np.searchsorted(va, e.lambda_mean, side="right") / len(va)
    ax.scatter([e.lambda_mean], [y], s=34, color=col, zorder=5, edgecolor="white", linewidth=0.6)
    ax.annotate(name, (e.lambda_mean, y), ha="left" if dx_ > 0 else "right", va="center", color=col,
                fontsize=8, xytext=(dx_, 6 if dx_ > 0 else 8), textcoords="offset points")
ax.set_xlim(0, 1.6)
ax.set_xlabel("Mean local synonymous rate, relative to the gene")
ax.set_ylabel("Fraction of exons at or below")
despine(ax)
panel(ax, "b", dx=-0.22, dy=1.02)
save_esm(fig, 10, "Position of MAPT exon 10 among 45,085 human coding exons of 75 to 115 nt, "
                  "by coverage with the synonymous constraint elements of Lin et al. (2011) (a) and by "
                  "mean local synonymous rate (b).")
plt.close(fig)

# ---------------------------------------------------------------- ¿están enriquecidos los alternativos?
from scipy import stats

print("\nExones cubiertos al menos 95 % por elementos de 30 codones, alternativos contra constitutivos")
enr = []
for lab, lo, hi in (("todos los largos", 0, 10**6), ("75 a 115 nt", 75, 115)):
    z = ex[(ex.n_windows > 0) & (ex.length >= lo) & (ex.length <= hi) &
           ex.kind.isin(["alternative", "constitutive"])]
    a_ = z[z.kind == "alternative"]; c_ = z[z.kind == "constitutive"]
    ta, tc = int((a_.cov30 >= 0.95).sum()), int((c_.cov30 >= 0.95).sum())
    odds, pv = stats.fisher_exact([[ta, len(a_) - ta], [tc, len(c_) - tc]])
    enr.append({"stratum": lab, "alternative": f"{ta}/{len(a_)}", "constitutive": f"{tc}/{len(c_)}",
                "odds_ratio": odds, "fisher_p": pv})
    print(f"  {lab}: alternativos {ta}/{len(a_)} ({ta / len(a_):.1%}), constitutivos {tc}/{len(c_)} "
          f"({tc / len(c_):.1%}), razón de momios {odds:.2f}, Fisher p = {pv:.1e}")
pd.DataFrame(enr).to_csv(f"{OUT}/genomewide_enrichment.tsv", sep="\t", index=False)
