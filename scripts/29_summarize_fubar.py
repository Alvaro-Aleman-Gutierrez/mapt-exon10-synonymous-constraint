"""Paso 29: resume FUBAR y responde la pregunta que FEL no puede responder.

FEL estima la tasa sinónima de cada sitio y la contrasta contra la tasa no sinónima del mismo sitio.
Nunca pregunta si la tasa sinónima de un sitio está por debajo del fondo del gen, que es la
afirmación central de este trabajo. FUBAR ajusta una rejilla de 20 por 20 pares (alfa, beta) con un
prior de Dirichlet común a todos los sitios y deja, en el archivo de caché, la verosimilitud
condicional de cada sitio en cada punto de la rejilla. Con eso y con los pesos de la rejilla se
obtiene la posterior de alfa para cada sitio, y de ahí la probabilidad de que la tasa sinónima de un
codón esté por debajo del fondo.

El fondo se define como la media de alfa bajo la rejilla ajustada, es decir el valor que el gen
tendría en un sitio del que no se supiera nada.
"""
import json
import re

import numpy as np
import pandas as pd
from scipy import stats

J = "results/selection/all.FUBAR.json"
CACHE = "results/mapt_2N4R_final.codon.fna.FUBAR.cache"
REP = ["E9", "E10", "E11", "E12", "E13"]

d = json.load(open(J))
grid = np.array(d["grid"])                      # columnas: alfa, beta, peso
alpha_g, beta_g, w_g = grid[:, 0], grid[:, 1], grid[:, 2]

txt = open(CACHE).read()
i = txt.index("conditionals", txt.index("conditionals") + 5)
# El caché trae, en este orden, la matriz de condicionales (un renglón por punto de la rejilla) y
# luego un renglón de escaladores del mismo largo, así que se toman los primeros len(grid).
rows = [r for r in re.findall(r"\{([0-9eE .,+-]{20,})\}", txt[i:])
        if len(r.split(",")) == 563][:len(grid)]
L = np.array([[float(x) for x in r.split(",")] for r in rows])   # rejilla por sitios
assert L.shape[0] == len(grid), (L.shape, len(grid))

post = L * w_g[:, None]
post /= post.sum(axis=0, keepdims=True)         # posterior de cada sitio sobre la rejilla

bg = float((w_g * alpha_g).sum() / w_g.sum())   # alfa de fondo del gen
site = pd.DataFrame({
    "column": np.arange(1, L.shape[1] + 1),
    "fubar_alpha": (post * alpha_g[:, None]).sum(axis=0),
    "fubar_beta": (post * beta_g[:, None]).sum(axis=0),
    "p_alpha_below_bg": post[alpha_g < bg].sum(axis=0),
    "p_alpha_below_tenth": post[alpha_g < bg / 10].sum(axis=0),
})
mle = pd.DataFrame(d["MLE"]["content"]["0"]).iloc[:, :6]
mle.columns = ["alpha", "beta", "beta_minus_alpha", "p_a_gt_b", "p_a_lt_b", "bf_a_lt_b"]
site["p_beta_gt_alpha"] = mle.p_a_lt_b.values

s = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
j = s.merge(site, on="column")
j = j[j.human_pos.notna()].copy()
j["human_pos"] = j.human_pos.astype(int)
j.to_csv("results/selection/fubar_sites.tsv", sep="\t", index=False)

print(f"alfa de fondo del gen bajo la rejilla ajustada: {bg:.3f}")
print(f"sitios con posición humana: {len(j)}\n")

print("media posterior de alfa por exón")
t = j.groupby("exon").agg(n=("fubar_alpha", "size"), alpha=("fubar_alpha", "median"),
                          beta=("fubar_beta", "median"))
print(t.round(3).to_string(), "\n")

e10 = j[j.exon == "E10"]
others = j[j.exon.isin([e for e in REP if e != "E10"])]
u, p = stats.mannwhitneyu(e10.fubar_alpha, others.fubar_alpha, alternative="less")
print(f"alfa en E10 contra las otras repeticiones: medianas {e10.fubar_alpha.median():.3f} y "
      f"{others.fubar_alpha.median():.3f}, Mann-Whitney p = {p:.2e}")
print(f"razón de medianas: {others.fubar_alpha.median() / e10.fubar_alpha.median():.1f} veces\n")

print("codones con posterior >= 0.9 de tener alfa por debajo del fondo del gen")
for ex in REP:
    g = j[j.exon == ex]
    print(f"  {ex}: {int((g.p_alpha_below_bg >= 0.9).sum())}/{len(g)} por debajo del fondo; "
          f"{int((g.p_alpha_below_tenth >= 0.9).sum())}/{len(g)} por debajo de la décima parte")
g = j[~j.exon.isin(REP)]
print(f"  resto del gen: {int((g.p_alpha_below_bg >= 0.9).sum())}/{len(g)} y "
      f"{int((g.p_alpha_below_tenth >= 0.9).sum())}/{len(g)}")

tab = pd.crosstab(j.exon == "E10", j.p_alpha_below_tenth >= 0.9)
if tab.shape == (2, 2):
    odds, pf = stats.fisher_exact(tab.values)
    print(f"\nE10 contra el resto del gen, alfa por debajo de la décima parte del fondo: "
          f"Fisher p = {pf:.2e}")

print("\nsitios que FUBAR marca como diversificantes con posterior >= 0.9")
div = j[j.p_beta_gt_alpha >= 0.9]
print(div[["human_pos", "exon", "fubar_alpha", "fubar_beta", "p_beta_gt_alpha",
           "fel_p", "meme_p", "meme_q"]].round(4).to_string(index=False))
print(f"\nfalsos positivos esperados a ese umbral, según el propio FUBAR: "
      f"{(1 - j.p_beta_gt_alpha).where(j.p_beta_gt_alpha >= 0.9).sum():.2f}")
print("results/selection/fubar_sites.tsv escrito")

# ---- Tabla 10. Resumen por exón, para el manuscrito
rows = []
for ex in [f"E{i}" for i in [1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13]]:
    g = j[j.exon == ex]
    rows.append([ex, len(g), f"{g.fubar_alpha.median():.3f}", f"{g.fubar_beta.median():.3f}",
                 f"{g.fubar_alpha.median() / g.fubar_beta.median():.1f}",
                 f"{int((g.p_alpha_below_bg >= 0.9).sum())}/{len(g)}",
                 f"{int((g.p_alpha_below_tenth >= 0.9).sum())}/{len(g)}"])
t10 = pd.DataFrame(rows, columns=["Exon", "Codons", "alpha", "beta", "alpha/beta",
                                  "alpha below background", "alpha below a tenth of background"])
t10.to_csv("results/tables/table4_fubar.tsv", sep="\t", index=False)
cap = ("Table 4. Site specific synonymous and nonsynonymous rates estimated by FUBAR. Values are "
       "posterior medians of the per site rates, alpha for synonymous and beta for nonsynonymous. "
       "The last two columns count codons whose posterior probability of carrying a synonymous rate "
       "below the gene background, and below a tenth of it, is at least 0.9. The background is the "
       "mean of alpha under the fitted grid, which is the value a site carries when nothing is known "
       "about it.")
with open("manuscript/tables/table4_fubar.md", "w") as f:
    f.write(f"**{cap}**\n\n" + t10.to_markdown(index=False) + "\n")
print("tabla 10 escrita")
