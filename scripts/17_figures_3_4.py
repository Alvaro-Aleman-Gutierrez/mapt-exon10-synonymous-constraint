"""Paso 17: figuras 3 y 4 del manuscrito, un archivo por figura.

Figura 3. Tasa sinónima e invariancia codón por codón en el exón 10, contra las demás
repeticiones. Figura 4. Identidad sinónima del exón 10 en mamíferos y en tetrápodos no mamíferos,
y la disociación entre conservación de la proteína y del nucleótido.

El estilo y las rutas de salida están en figstyle.py, con las especificaciones de arte del envío.
Los títulos descriptivos van en los pies de figura, no dentro de la ilustración.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon

sys.path.insert(0, "scripts")
from figstyle import FOCUS, REF, MARK, INK, MUTED, GRID, MM, W2, panel, despine, save

rng = np.random.default_rng(7)
REPS = ["E9", "E10", "E11", "E12", "E13"]
CLADES = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Caudata", "Anura"]
AMNIOTA = {"Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria"}
FTDP = [(279, "N279K", "splicing"), (280, "ΔK280", "splicing"),
        (284, "L284L", "splicing"), (296, "N296N", "splicing"),
        (301, "P301L", "protein"), (303, "G303V", "protein"),
        (305, "S305N", "splicing")]

sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()].copy()
sites["human_pos"] = sites.human_pos.astype(int)
e10 = sites[sites.exon == "E10"].sort_values("human_pos")
tet = pd.read_csv("results/tetrapod_exon_identity.tsv", sep="\t")

rows = []
for g in CLADES:
    sub = tet[tet.group == g]
    syn = sub.pivot(index="species", columns="exon", values="syn_identity").dropna()
    aa = sub.pivot(index="species", columns="exon", values="aa_identity").dropna()
    nt = sub.pivot(index="species", columns="exon", values="nt_identity").dropna()
    other = syn[[c for c in REPS if c != "E10"]].mean(axis=1)
    rows.append({"clade": g, "n": len(syn), "syn_e10": syn.E10.mean(),
                 "syn_other": other.mean(), "higher": int((syn.E10 > other).sum()),
                 "p": wilcoxon(syn.E10, other).pvalue,
                 "d_aa": aa.E10.mean() - aa.E11.mean(),
                 "d_nt": nt.E10.mean() - nt.E11.mean()})
ct = pd.DataFrame(rows).set_index("clade")

# ================================================================ figura 3
fig = plt.figure(figsize=(W2, 132 * MM))
gs = fig.add_gridspec(4, 1, height_ratios=[1.3, 0.1, 0.34, 1.0], hspace=0.07,
                      left=0.175, right=0.985, top=0.955, bottom=0.085)

axA = fig.add_subplot(gs[0])
x = np.arange(len(e10))
axA.bar(x, e10.fel_dS, width=0.66, color=REF, linewidth=0, zorder=3)
axA.set_xlim(-0.8, len(e10) - 0.2)
axA.set_ylim(0, 1.95)
axA.set_yticks([0, 0.5, 1.0, 1.5])
axA.set_ylabel("Synonymous rate (dS)")
axA.yaxis.grid(True, color=GRID, zorder=0)
axA.set_axisbelow(True)
axA.set_xticks([])
despine(axA, keep=("left",))
panel(axA, "a", dx=-0.215)

for pos, name, kind in FTDP:
    i = int(np.where(e10.human_pos.values == pos)[0][0])
    col = MARK if kind == "splicing" else INK
    y = e10.fel_dS.values[i]
    axA.plot([i], [y + 0.09], marker="v", ms=4, color=col, markeredgewidth=0, zorder=4)
    axA.annotate(name, (i, y + 0.17), ha="center", va="bottom", color=col, rotation=90)

leg = axA.legend(handles=[
    Line2D([], [], marker="v", ls="", ms=4, color=MARK, label="acts on splicing"),
    Line2D([], [], marker="v", ls="", ms=4, color=INK, label="acts on the protein"),
], loc="upper left", frameon=False, handletextpad=0.15, borderpad=0.1,
   labelspacing=0.35, bbox_to_anchor=(0.0, 1.02), title="FTDP-17 mutation")
leg.get_title().set_fontsize(8)
leg._legend_box.align = "left"

# pista de codones invariantes
axT = fig.add_subplot(gs[1], sharex=axA)
axT.set_ylim(0, 1)
for s in axT.spines.values():
    s.set_visible(False)
axT.set_yticks([])
for i, inv in enumerate(e10.codon_invariant.values):
    axT.add_patch(plt.Rectangle((i - 0.33, 0.05), 0.66, 0.9, zorder=3,
                                facecolor=FOCUS if inv else "white",
                                edgecolor=FOCUS if inv else "#b0b0b0", linewidth=0.6))
axT.set_ylabel("Codon identical in\nall 259 mammals", rotation=0, ha="right", va="center",
               labelpad=16, linespacing=1.4)
axT.set_xticks(x[::3])
axT.set_xticklabels(e10.human_pos.values[::3])
axT.tick_params(axis="x", length=2.5, pad=2)
axT.set_xlabel("Codon position in human 2N4R tau", labelpad=4)

# distribución de dS por exón de repetición
axB = fig.add_subplot(gs[3])
# Techo en 4.5, que deja fuera 3 sitios de 255 en vez de los 42 que recortaba un techo de 2.
# Los que quedan fuera se marcan con un caret para que no se lean como si valieran 4.5.
TOP = 4.5
over = []
for j, e in enumerate(REPS):
    v = sites[sites.exon == e].fel_dS.values
    col = FOCUS if e == "E10" else REF
    jit = rng.uniform(-0.17, 0.17, len(v))
    over.append((e, int((v > TOP).sum())))
    axB.scatter(j + jit, np.minimum(v, TOP), s=11, facecolor=col, edgecolor="white",
                linewidth=0.35, zorder=3, clip_on=True)
    hi = v > TOP
    if hi.any():
        axB.plot(j + jit[hi], np.full(hi.sum(), TOP), marker="^", ms=3.2, ls="",
                 color=col, markeredgewidth=0, clip_on=False, zorder=5)
    med = np.median(v)
    axB.plot([j - 0.32, j + 0.32], [med, med], color=col, linewidth=1.8, zorder=4,
             solid_capstyle="butt")
    axB.annotate(f"{med:.2f}" if med >= 0.01 else "<0.01", (j + 0.36, med), va="center",
                 color=col)
axB.set_xticks(range(len(REPS)))
axB.set_xticklabels([f"exon {e[1:]}\nn = {(sites.exon == e).sum()}" for e in REPS])
axB.tick_params(axis="x", length=0, pad=4)
axB.set_xlim(-0.62, len(REPS) - 0.2)
axB.set_ylim(-0.12, 4.72)
axB.set_yticks([0, 1, 2, 3, 4])
axB.set_ylabel("Synonymous rate (dS)")
axB.yaxis.grid(True, color=GRID, zorder=0)
axB.set_axisbelow(True)
despine(axB)
panel(axB, "b", dx=-0.215, dy=1.015)

save(fig, 3)
plt.close(fig)
print("  sitios fuera de escala en el panel b:", over)

# ================================================================ figura 4
# Las leyendas van fuera del área de datos, encima de cada panel, para que no se confundan con los
# puntos. Eso obliga a bajar el tope y a dar más alto a la figura.
fig = plt.figure(figsize=(W2, 94 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.44,
                      left=0.115, right=0.975, top=0.845, bottom=0.165)

axA = fig.add_subplot(gs[0])
ys = np.arange(len(CLADES))[::-1]
for y, g in zip(ys, CLADES):
    r = ct.loc[g]
    axA.plot([r.syn_other, r.syn_e10], [y, y], color="#b0b0b0", linewidth=1.2, zorder=2,
             solid_capstyle="round")
    # forma además de color, para que la figura se lea también en blanco y negro
    axA.scatter([r.syn_other], [y], s=32, marker="s", facecolor=REF, edgecolor="white",
                linewidth=0.4, zorder=3)
    axA.scatter([r.syn_e10], [y], s=38, marker="o", facecolor=FOCUS, edgecolor="white",
                linewidth=0.4, zorder=3)
    axA.annotate(f"{int(r.higher)}/{int(r.n)}", (max(r.syn_e10, r.syn_other) + 0.025, y),
                 va="center", color=MUTED)
axA.set_yticks(ys)
axA.set_yticklabels(CLADES)
axA.tick_params(axis="y", length=0)
axA.set_ylim(-1.75, len(CLADES) - 0.2)
axA.set_xlim(0.1, 1.14)
axA.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
axA.set_xlabel("Synonymous identity to human")
axA.xaxis.grid(True, color=GRID, zorder=0)
axA.set_axisbelow(True)
despine(axA, keep=("bottom",))
panel(axA, "a", dx=-0.19, dy=1.10)
axA.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=5.5, markerfacecolor=FOCUS, markeredgecolor="white",
           color="none", label="exon 10"),
    Line2D([], [], marker="s", ls="", ms=5, markerfacecolor=REF, markeredgecolor="white",
           color="none", label="mean of exons 9, 11, 12, 13"),
], loc="lower left", frameon=False, handletextpad=0.15, borderpad=0.0,
   columnspacing=1.8, ncol=2, bbox_to_anchor=(0.0, 1.075))

axB = fig.add_subplot(gs[1])
axB.axhline(0, color="#b0b0b0", linewidth=0.7, zorder=2)
axB.axvline(0, color="#b0b0b0", linewidth=0.7, zorder=2)
# cada etiqueta se aleja lo suficiente de su punto y elige el lado que no cruza el eje ni el cero
OFF = {"Mammalia": (-0.013, 0, "right"), "Aves": (0.013, 0, "left"),
       "Crocodylia": (0.013, 0, "left"), "Testudines": (-0.013, 0, "right"),
       "Lepidosauria": (-0.013, 0, "right"), "Caudata": (0.013, 0, "left"),
       "Anura": (0.013, 0, "left")}
for g in CLADES:
    r = ct.loc[g]
    amn = g in AMNIOTA
    axB.scatter([r.d_aa], [r.d_nt], s=44, zorder=4, linewidth=1.0, marker="o" if amn else "D",
                facecolor=FOCUS if amn else "white", edgecolor=FOCUS if amn else MUTED)
    dx, dy, ha = OFF[g]
    axB.annotate(g, (r.d_aa + dx, r.d_nt + dy), ha=ha, color=INK,
                 va="bottom" if ha == "center" else "center")
axB.set_xlabel("Amino acid identity,\nexon 10 minus exon 11", linespacing=1.4)
axB.set_ylabel("Nucleotide identity,\nexon 10 minus exon 11", linespacing=1.4)
axB.set_xlim(-0.225, 0.105)
axB.set_ylim(-0.05, 0.105)
axB.set_xticks([-0.15, -0.10, -0.05, 0, 0.05])
axB.set_yticks([-0.04, 0, 0.04, 0.08])
despine(axB)
panel(axB, "b", dx=-0.24, dy=1.10)
axB.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=5.5, color=FOCUS, label="Amniota"),
    Line2D([], [], marker="D", ls="", ms=5, markerfacecolor="white", markeredgecolor=MUTED,
           color="none", label="Amphibia"),
], loc="lower left", frameon=False, handletextpad=0.15, borderpad=0.0,
   columnspacing=1.8, ncol=2, bbox_to_anchor=(0.0, 1.075))

save(fig, 4)
plt.close(fig)
print(ct.round(3).to_string())
