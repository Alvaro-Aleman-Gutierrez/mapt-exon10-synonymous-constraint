"""Paso 32: pruebas sobre unidades independientes, para no tratar especies emparentadas como réplicas.

Los Wilcoxon por especie de la tabla 4 comparan el exón 10 con las otras repeticiones dentro de cada
especie, y eso corrige la distancia al humano, pero las especies de un mismo orden comparten
ancestría, así que 255 mamíferos no son 255 observaciones. Aquí cada orden (mamíferos y aves) o
cada familia (los demás clados) cuenta una vez, con la media de la diferencia entre sus especies, y
se pregunta con una prueba de signos si el exón 10 supera a las otras repeticiones en más unidades
de las que el azar daría. Las unidades no son independientes del todo, pero lo son mucho más que
las especies, y llevan el mismo tipo de pregunta a un nivel donde la prueba es honesta.

También se corre una PGLS de movimiento browniano sobre los mamíferos, que es lo que un revisor con
métodos comparados probaría primero. Se reporta aunque no salga significativa, junto con la razón
por la que creemos que el modelo no se ajusta a este rasgo.
"""
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats

REP = ["E9", "E10", "E11", "E12", "E13"]
OTHERS = [e for e in REP if e != "E10"]
GROUPS = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Caudata", "Anura"]
UNIT_OF = {"Mammalia": "order", "Aves": "order"}      # el resto, por familia

t = pd.read_csv("results/tetrapod_exon_identity.tsv", sep="\t")
orth = pd.read_csv("data/tetrapod_orthologs.tsv", sep="\t")
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
mam_order = dict(zip(sel.species.str.replace(" ", "_"), sel["order"]))
order_of = dict(zip(orth.species, orth["order"]))
family_of = dict(zip(orth.species, orth.family))


def unit(group, sp):
    if group == "Mammalia":
        return mam_order.get(sp.replace(" ", "_"))
    return order_of.get(sp) if UNIT_OF.get(group) == "order" else family_of.get(sp)


def diffs(group):
    w = (t[t.group == group].pivot_table(index="species", columns="exon", values="syn_identity",
                                         aggfunc="first").dropna(subset=REP))
    return w.E10 - w[OTHERS].mean(axis=1)


rows = []
for g in GROUPS:
    d = diffs(g)
    u = pd.DataFrame({"d": d, "u": [unit(g, s) for s in d.index]}).dropna().groupby("u").d.mean()
    pos, n = int((u > 0).sum()), len(u)
    p = stats.binomtest(pos, n, 0.5, alternative="greater").pvalue
    rows.append({"clade": g, "species": len(d), "unit": UNIT_OF.get(g, "family"), "units": n,
                 "units_E10_higher": pos, "sign_test_p": p, "mean_difference": u.mean()})
res = pd.DataFrame(rows)
res.to_csv("results/independent_units.tsv", sep="\t", index=False)
print(res.to_string(index=False))

# ---- el mismo análisis, comparando E10 solo con E11 y E12, que son las repeticiones homólogas
# E9 contiene la región rica en prolina además de la primera repetición y E13 es la cola carboxilo
# terminal, así que no son repeticiones. Si el resultado dependiera de ellos, cambiaría aquí.
def diffs_rep(group):
    w = (t[t.group == group].pivot_table(index="species", columns="exon", values="syn_identity",
                                         aggfunc="first").dropna(subset=REP))
    return w.E10 - w[["E11", "E12"]].mean(axis=1)


rows2 = []
for g in GROUPS:
    d = diffs_rep(g)
    u = pd.DataFrame({"d": d, "u": [unit(g, s) for s in d.index]}).dropna().groupby("u").d.mean()
    pos, n = int((u > 0).sum()), len(u)
    rows2.append({"clade": g, "species": len(d), "unit": UNIT_OF.get(g, "family"), "units": n,
                  "units_E10_higher": pos,
                  "sign_test_p": stats.binomtest(pos, n, 0.5, alternative="greater").pvalue,
                  "mean_difference": u.mean()})
res2 = pd.DataFrame(rows2)
res2.to_csv("results/independent_units_E11_E12.tsv", sep="\t", index=False)
print("\ncontroles E11 y E12 solamente")
print(res2.to_string(index=False))

# ---- comparaciones de exones alternativos y de E10 contra E11, por orden
d6 = pd.read_csv("results/tables/ESM_9_alt_exon_divergence.tsv", sep="\t", comment="#")
w6 = d6.pivot_table(index=["species", "order"], columns="exon", values="identity_to_human").reset_index()
print("\nE6 contra exones cubiertos, por orden")
for other in ("E1", "E5"):
    s = w6[["order", "E6", other]].dropna()
    s = s.assign(d=s.E6 - s[other]).groupby("order").d.agg(["mean", "size"])
    pos, n = int((s["mean"] > 0).sum()), len(s)
    p = stats.binomtest(pos, n, 0.5, alternative="greater").pvalue
    print(f"  E6 contra {other}: {pos} de {n} órdenes con E6 mayor, p = {p:.1e}, "
          f"diferencia media {s['mean'].mean():+.3f}")

d7 = pd.read_csv("results/tables/ESM_10_exon_synonymous_identity.tsv", sep="\t", comment="#")
m = d7[d7.order.isin(["Rodentia", "Artiodactyla"])].pivot_table(
    index=["species", "order"], columns="exon", values="syn_identity").reset_index()
print("\nE10 contra E11 en los dos órdenes igualados")
for o, g in m.groupby("order"):
    g = g[["E10", "E11"]].dropna()
    print(f"  {o}: {len(g)} especies, diferencia media {(g.E10 - g.E11).mean():+.3f}, "
          f"E10 mayor en {int((g.E10 > g.E11).sum())}")

# ---- PGLS browniana sobre los mamíferos
tr = Phylo.read("results/phylo_protein/mapt_2N4R.treefile", "newick")
tr.root_with_outgroup([l for l in tr.get_terminals() if mam_order.get(l.name) == "Monotremata"])
y = diffs("Mammalia")
y.index = [s.replace(" ", "_") for s in y.index]
tips = [l.name for l in tr.get_terminals() if l.name in y.index]
idx = {n: i for i, n in enumerate(tips)}
n = len(tips)
C = np.zeros((n, n))
depth = tr.depths()


def fill(c):
    """Cada par de hojas comparte la profundidad de su ancestro común más reciente."""
    if c.is_terminal():
        return [idx[c.name]] if c.name in idx else []
    parts = [fill(k) for k in c.clades]
    for a in range(len(parts)):
        for b in range(a + 1, len(parts)):
            for i in parts[a]:
                for j in parts[b]:
                    C[i, j] = C[j, i] = depth[c]
    return [i for p_ in parts for i in p_]


fill(tr.root)
for l in tr.get_terminals():
    if l.name in idx:
        C[idx[l.name], idx[l.name]] = depth[l]
C += np.eye(n) * 1e-6
v = y.loc[tips].values
Ci = np.linalg.inv(C)
one = np.ones(n)
den = one @ Ci @ one
mu = (one @ Ci @ v) / den
se = np.sqrt(((v - mu) @ Ci @ (v - mu)) / (n - 1) / den)
p_pgls = 2 * stats.t.sf(abs(mu / se), n - 1)

hum = next(l for l in tr.get_terminals() if l.name == "Homo_sapiens")
dist = pd.Series({l.name: tr.distance(hum, l) for l in tr.get_terminals() if l.name in idx})
rho, p_rho = stats.spearmanr(dist.loc[tips], v)
pd.DataFrame([{"n": n, "pgls_mean_difference": mu, "pgls_se": se, "pgls_p": p_pgls,
               "simple_mean_difference": v.mean(), "rho_distance_to_human": rho,
               "rho_p": p_rho}]).to_csv("results/pgls_mammals.tsv", sep="\t", index=False)
print(f"\nPGLS mamíferos (n = {n}): {mu:+.3f} (SE {se:.3f}), p = {p_pgls:.2f}; "
      f"diferencia contra distancia a humano, rho = {rho:.2f}, p = {p_rho:.1e}")
