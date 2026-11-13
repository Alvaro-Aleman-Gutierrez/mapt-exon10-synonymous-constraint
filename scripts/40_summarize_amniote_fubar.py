"""Paso 40: resume el FUBAR de cada clado (paso 39) y aplica los criterios fijados en el plan.

Mismo procedimiento que el paso 29. De la caché de FUBAR salen las verosimilitudes condicionales de
cada sitio sobre la rejilla, con los pesos de la rejilla se obtiene la posterior de cada sitio, y el
fondo es la media de alfa bajo la rejilla ajustada. Aquí el fondo es el de los exones 9 a 13 del
clado, no el del gen completo, así que los valores no se comparan con la Tabla 4 en términos
absolutos sino exón contra exón dentro de cada clado.

Solo entran los sitios homólogos a un residuo humano. Para la pregunta de los mismos sitios se
compara, posición por posición del exón 10, la tasa sinónima de cada clado con la de mamíferos, y se
cuenta cuántos de los 14 codones invariantes en mamíferos caen en el tercio más lento del exón 10 de
cada clado, contra 10.000 permutaciones de 14 posiciones dentro del exón.

Criterios del plan (PLAN_punto1_amniotas.md), fijados antes de correr: hay apoyo si en aves y en al
menos uno de escamosos o tortugas el exón 10 tiene la menor tasa sinónima mediana de los cinco exones,
con Mann-Whitney p < 0,01, y una tasa no sinónima que no es menor que la del exón 11.
"""
import json
import re

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats

OUT = "results/amniote"
REP = ["E9", "E10", "E11", "E12", "E13"]
CLADES = ["Mammalia", "Aves", "Lepidosauria", "Testudines", "Archelosauria", "Anura", "Caudata"]
INVARIANT = [275, 276, 279, 280, 281, 283, 285, 287, 289, 292, 294, 296, 297, 304]
E10_START = 275                                   # numeración de 2N4R, exón 10 = 275 a 305
rng = np.random.default_rng(1)

# sitios sinónimos del codón humano en cada posición del exón 10 (Nei y Gojobori), para controlar la
# correlación entre clados por degeneración del código: un codón de cuatro vías tiene más sitios
# sinónimos que uno de dos en cualquier linaje
from Bio import SeqIO
from Bio.Data import CodonTable
_tab = CodonTable.unambiguous_dna_by_id[1].forward_table
_hc = {r.id: str(r.seq) for r in SeqIO.parse("data/mapt_2N4R_cds.fna", "fasta")}["Homo_sapiens"]

def _syn(c):
    return sum(1 / 3 for i in range(3) for x in "ACGT" if x != c[i]
               and (c[:i] + x + c[i + 1:]) in _tab and _tab[c[:i] + x + c[i + 1:]] == _tab[c])

SYN = pd.Series({j: _syn(_hc[3 * (E10_START - 1 + j):3 * (E10_START + j)]) for j in range(31)})

def site_posteriors(clade, tag):
    d = json.load(open(f"{OUT}/{clade}.{tag}.FUBAR.json"))
    grid = np.array(d["grid"])
    a, b, w = grid[:, 0], grid[:, 1], grid[:, 2]
    n = len(pd.read_csv(f"{OUT}/{clade}.sites.tsv", sep="\t"))
    txt = open(f"{OUT}/{clade}.{tag}.FUBAR.cache").read()
    i = txt.index("conditionals", txt.index("conditionals") + 5)
    rows = [r for r in re.findall(r"\{([0-9eE .,+-]{20,})\}", txt[i:])
            if len(r.split(",")) == n][:len(grid)]
    L = np.array([[float(x) for x in r.split(",")] for r in rows])
    assert L.shape == (len(grid), n), (clade, L.shape)
    post = L * w[:, None]
    post /= post.sum(axis=0, keepdims=True)
    bg = float((w * a).sum() / w.sum())
    s = pd.read_csv(f"{OUT}/{clade}.sites.tsv", sep="\t")
    s["alpha"] = (post * a[:, None]).sum(axis=0)
    s["beta"] = (post * b[:, None]).sum(axis=0)
    s["p_below_bg"] = post[a < bg].sum(axis=0)
    s = s[s.human_residue.notna() & (s.human_residue != "-")].copy()
    # posición humana dentro de cada exón, para comparar sitio por sitio entre clados
    s["pos_in_exon"] = s.groupby("exon").cumcount()
    return s, bg

def summarize(tag):
    rows, e10 = [], {}
    for c in CLADES:
        if c == "Mammalia" and tag == "free":
            continue
        try:
            s, bg = site_posteriors(c, tag)
        except FileNotFoundError:
            continue
        med = s.groupby("exon").alpha.median()
        medb = s.groupby("exon").beta.median()
        x, y = s[s.exon == "E10"].alpha, s[s.exon.isin(["E9", "E11", "E12", "E13"])].alpha
        p = stats.mannwhitneyu(x, y, alternative="less").pvalue
        pb = stats.mannwhitneyu(s[s.exon == "E10"].beta, s[s.exon == "E11"].beta).pvalue
        below = {e: f"{int((s[s.exon == e].p_below_bg >= 0.9).sum())}/{(s.exon == e).sum()}"
                 for e in REP + ["E10b"] if (s.exon == e).any()}
        tl = Phylo.read(f"{OUT}/{c}.constrained.treefile" if c != "Mammalia" else
                        "results/phylo_protein/mapt_2N4R_constrained.treefile",
                        "newick").total_branch_length()
        rows.append({"clade": c, "species": sum(1 for l in open(f"{OUT}/{c}.codon.fna")
                                                 if l.startswith(">")),
                     "background_alpha": round(bg, 3),
                     **{f"alpha_{e}": round(med[e], 3) for e in REP + ["E10b"] if e in med},
                     "alpha_others": round(y.median(), 3),
                     "fold": round(y.median() / x.median(), 1),
                     "mw_p_E10_lower": p,
                     "E10_lowest_of_five": med[REP].idxmin() == "E10",
                     "beta_E10": round(medb["E10"], 3), "beta_E11": round(medb["E11"], 3),
                     "beta_E10_vs_E11_p": pb,
                     **{f"below_bg_{e}": v for e, v in below.items()},
                     "tree_length": round(tl, 2)})
        e10[c] = s[s.exon == "E10"].set_index("pos_in_exon").alpha
    t = pd.DataFrame(rows)

    # los mismos sitios: correlación con mamíferos y codones invariantes en el tercio más lento
    if "Mammalia" in e10:
        inv = [p - E10_START for p in INVARIANT]
        same = []
        for c, a in e10.items():
            if c == "Mammalia":
                continue
            j = pd.concat([e10["Mammalia"], a], axis=1, keys=["m", "c"]).dropna()
            rho, prho = stats.spearmanr(j.m, j.c)
            z = SYN[j.index]
            rx, ry, rz = (stats.rankdata(v) for v in (j.m, j.c, z))
            ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
            ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
            prt, pprt = stats.pearsonr(ex, ey)
            slow = set(a.nsmallest(int(round(len(a) / 3))).index)
            k = len(slow & set(inv))
            null = [len(slow & set(rng.choice(a.index, len(inv), replace=False)))
                    for _ in range(10000)]
            same.append({"clade": c, "spearman_rho_with_mammals": round(rho, 3),
                         "spearman_p": prho, "partial_rho_given_syn_sites": round(prt, 3),
                         "partial_p": pprt, "invariant_in_slowest_third": f"{k}/{len(inv)}",
                         "expected": round(float(np.mean(null)), 2),
                         "perm_p": (1 + sum(n >= k for n in null)) / 10001})
        t = t.merge(pd.DataFrame(same), on="clade", how="left")
    return t

import esm
both = []
for tag in ["constrained", "free"]:
    t = summarize(tag)
    if t.empty:
        continue
    t.to_csv(f"{OUT}/fubar_by_clade.{tag}.tsv", sep="\t", index=False)
    both.append(t.assign(tree="taxonomy constrained" if tag == "constrained" else "unconstrained"))
    print(f"\n==== árbol {tag}")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(t.to_string(index=False))

t = pd.read_csv(f"{OUT}/fubar_by_clade.constrained.tsv", sep="\t").set_index("clade")
def ok(c):
    r = t.loc[c]
    return bool(r.E10_lowest_of_five and r.mw_p_E10_lower < 0.01 and r.beta_E10 >= r.beta_E11)
verdict = ok("Aves") and (ok("Lepidosauria") or ok("Testudines"))
print("\ncriterio del plan:", "APOYO" if verdict else
      ("APOYO PARCIAL" if ok("Aves") else "SIN APOYO"),
      {c: ok(c) for c in ["Aves", "Lepidosauria", "Testudines", "Archelosauria", "Anura",
                          "Caudata"] if c in t.index})

# Online Resource 14, los dos árboles
esm.write_esm(pd.concat(both)[["tree"] + [c for c in both[0].columns]], 14, "fubar_by_clade")
