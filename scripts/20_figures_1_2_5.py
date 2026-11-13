"""Paso 20: figuras 1, 2 y 5 del manuscrito, y el árbol completo como recurso en línea.

Figura 1. Filogenia de máxima verosimilitud de las 259 secuencias, resumida por orden, con el
soporte de rama y los tres grupos que no se recuperan. Figura 2. Selección a lo largo de la
proteína, con el mapa de exones a escala. Figura 5. Locus de MAPT en los dos damanes y el mapa de
inserciones y deleciones del alineamiento afroterio.

El árbol con las 259 etiquetas no cabe en una figura legible, así que va aparte como recurso en
línea, en PDF de página alta.
"""
import io
import sys
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import AlignIO, Phylo, SeqIO

sys.path.insert(0, "scripts")
from figstyle import FOCUS, REF, MARK, INK, MUTED, GRID, MM, W2, panel, despine, save, save_esm

CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
MTBR = ["E9", "E10", "E11", "E12", "E13"]

sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
order_of = {s.replace(" ", "_"): o for s, o in zip(sel.species, sel["order"])}
hx = pd.read_csv("data/human_exons.tsv", sep="\t").set_index("exon")
sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()].copy()
sites["human_pos"] = sites.human_pos.astype(int)

# ================================================================ figura 1
tree = Phylo.read("results/phylo_protein/mapt_2N4R.treefile", "newick")
tree.root_with_outgroup([l for l in tree.get_terminals()
                         if order_of[l.name] == "Monotremata"])

counts = collections.Counter(order_of.values())
# orden filogenético aproximado, de monotremas a primates
SEQ = ["Monotremata", "Didelphimorphia", "Microbiotheria", "Dasyuromorphia", "Diprotodontia",
       "Afrosoricida", "Macroscelidea", "Tubulidentata", "Proboscidea", "Sirenia",
       "Cingulata", "Pilosa", "Lagomorpha", "Rodentia", "Scandentia", "Dermoptera",
       "Primates", "Eulipotyphla", "Chiroptera", "Perissodactyla", "Pholidota",
       "Carnivora", "Artiodactyla"]
assert set(SEQ) == set(counts), set(SEQ) ^ set(counts)

def clade_support(members):
    leaves = [l for l in tree.get_terminals() if l.name in members]
    if len(leaves) < 2:
        return None, True, []
    mrca = tree.common_ancestor(leaves)
    inside = {l.name for l in mrca.get_terminals()}
    # el árbol se enraíza en Monotremata, así que su ancestro común es la raíz y la prueba de
    # monofilia no aplica al grupo externo
    if mrca is tree.root:
        return None, True, []
    intruders = sorted(inside - set(members))
    sup = mrca.confidence
    if sup is None and mrca.name:
        sup = mrca.name
    return sup, not intruders, intruders

rows = []
for o in SEQ:
    mem = [k for k, v in order_of.items() if v == o]
    sup, mono, intr = clade_support(mem)
    ufb = shl = None
    if isinstance(sup, str) and "/" in sup:
        ufb, shl = sup.split("/")[:2]
    elif sup is not None:
        ufb = sup
    rows.append({"order": o, "n": counts[o], "mono": mono, "ufboot": ufb, "shalrt": shl,
                 "intruders": intr})
f1 = pd.DataFrame(rows)
f1.to_csv("results/phylo_protein/order_support.tsv", sep="\t", index=False)

# El árbol se dibuja completo, sin etiquetas de punta, porque 259 nombres no caben a 8 pt. Cada
# orden lleva una llave en el margen derecho sobre el tramo de puntas que ocupa, y los nombres se
# separan verticalmente para que no se pisen. Un orden que el árbol no recupera aparece partido en
# varios tramos, y eso se ve directamente.
from matplotlib.collections import LineCollection

leaves = tree.get_terminals()
n = len(leaves)
parents = {c: p for p in tree.find_clades() for c in p.clades}

dep = {}
def set_depth(node, d=0.0):
    dep[id(node)] = d
    for c in node.clades:
        set_depth(c, d + (c.branch_length or 0.0))
set_depth(tree.root)

ycoord, segs = {}, []
def walk(node, counter=[0]):
    if node.is_terminal():
        y = n - 1 - counter[0]
        counter[0] += 1
    else:
        kids = [walk(c) for c in node.clades]
        y = (min(kids) + max(kids)) / 2
        x = dep[id(node)]
        segs.append([(x, min(kids)), (x, max(kids))])
        for c, cy in zip(node.clades, kids):
            segs.append([(x, cy), (dep[id(c)], cy)])
    ycoord[id(node)] = y
    return y
walk(tree.root)

def c_support(node):
    """Devuelve (UFBoot, SH-aLRT). IQ-TREE las escribe juntas en el nombre del nodo."""
    s = node.confidence if node.confidence is not None else node.name
    try:
        parts = str(s).split("/")
        return float(parts[0]), (float(parts[1]) if len(parts) > 1 else None)
    except (TypeError, ValueError):
        return None, None

fig = plt.figure(figsize=(W2, 215 * MM))
ax = fig.add_subplot(111)
ax.add_collection(LineCollection(segs, colors=MUTED, linewidths=0.35, zorder=2))

xmax = max(dep[id(l)] for l in leaves)
for node in tree.get_nonterminals():
    u, a = c_support(node)
    if u is not None and u >= 95 and (a is None or a >= 80) and node.count_terminals() >= 5:
        ax.plot([dep[id(node)]], [ycoord[id(node)]], marker="o", ms=3.2, color=INK,
                markeredgecolor="white", markeredgewidth=0.35, zorder=4)

# tramos contiguos de puntas del mismo orden
runs, cur = [], None
for l in leaves:
    o = order_of[l.name]
    if cur and cur[0] == o:
        cur[2] = ycoord[id(l)]
    else:
        cur = [o, ycoord[id(l)], ycoord[id(l)]]
        runs.append(cur)
nruns = collections.Counter(r[0] for r in runs)
# Diprotodontia tiene sus puntas seguidas y aun así no es monofilético, porque su ancestro común
# incluye dasiuromorfos y Dromiciops. La contigüidad no basta, así que el color sale de f1.mono.
mono_of = dict(zip(f1["order"], f1.mono))
OUTGROUP = "Monotremata"

xb = xmax * 1.05
for o, y1, y2 in runs:
    lo, hi = min(y1, y2), max(y1, y2)
    ax.plot([xb, xb], [lo - 0.45, hi + 0.45],
            color=REF if (mono_of[o] or o == OUTGROUP) else FOCUS,
            linewidth=4.2, solid_capstyle="butt", zorder=3)

# un nombre por orden, colocado sobre su tramo mayor y luego separado de los vecinos
labels = []
for o in {r[0] for r in runs}:
    cand = [r for r in runs if r[0] == o]
    big = max(cand, key=lambda r: abs(r[2] - r[1]))
    bad = not mono_of[o] and o != OUTGROUP
    note = ""
    if o == OUTGROUP:
        note = ", outgroup"
    elif bad:
        note = ", not monophyletic"
    labels.append({"y0": (min(big[1], big[2]) + max(big[1], big[2])) / 2,
                   "anchors": [(min(r[1], r[2]) + max(r[1], r[2])) / 2 for r in cand],
                   "text": f"{o} ({counts[o]}){note}",
                   "col": FOCUS if bad else INK})
labels.sort(key=lambda d: d["y0"])
GAP = 4.6
for i in range(1, len(labels)):
    labels[i]["y"] = max(labels[i]["y0"], labels[i - 1].get("y", labels[i - 1]["y0"]) + GAP)
labels[0]["y"] = labels[0]["y0"]
for i in range(1, len(labels)):
    labels[i]["y"] = max(labels[i]["y0"], labels[i - 1]["y"] + GAP)
# no dejar que la pila se salga por arriba
over = labels[-1]["y"] - (n - 1)
if over > 0:
    for d in labels:
        d["y"] -= over
for i in range(len(labels) - 2, -1, -1):
    labels[i]["y"] = min(labels[i]["y"], labels[i + 1]["y"] - GAP)

xt = xb + xmax * 0.055
for d in labels:
    for a in d["anchors"]:
        ax.plot([xb + xmax * 0.012, xt - xmax * 0.012], [a, d["y"]], color="#c8c8c8",
                linewidth=0.4, zorder=2)
    ax.annotate(d["text"], (xt, d["y"]), va="center", color=d["col"])

ax.set_xlim(-xmax * 0.01, xmax * 1.70)
ax.set_ylim(-4, n + 3)
ax.set_yticks([])
ax.set_xlabel("substitutions per site")
despine(ax, keep=("bottom",))
ax.legend(handles=[
    Line2D([], [], color=REF, lw=4.2, label="order recovered as a clade"),
    Line2D([], [], color=FOCUS, lw=4.2, label="order not recovered"),
    Line2D([], [], marker="o", ls="", ms=4.2, color=INK,
           label="UFBoot \u2265 95 and SH-aLRT \u2265 80, five or more tips"),
], loc="upper left", frameon=False, handletextpad=0.6, borderpad=0.2, labelspacing=0.45,
   bbox_to_anchor=(0.0, 1.0))
fig.subplots_adjust(left=0.02, right=0.995, top=0.995, bottom=0.042)

save(fig, 1, tiff_dpi=450)
plt.close(fig)

# -- recurso en línea, árbol completo con etiquetas y soporte de nodo
# Se dibuja a mano en vez de con Phylo.draw, porque esa función no deja poner círculos
# proporcionales al soporte y además imprime los nombres de nodo interno, que en un árbol de
# IQ-TREE son los propios valores de soporte y acaban encima de los nombres de especie.
n = len(leaves)
figT = plt.figure(figsize=(170 * MM, (n * 3.1 + 34) * MM))
axT = figT.add_subplot(111)
axT.add_collection(LineCollection(segs, colors=MUTED, linewidths=0.45, zorder=2))

# El tamaño y la opacidad del círculo codifican UFBoot, de 50 a 100. El relleno añade la segunda
# medida: lleno cuando SH-aLRT alcanza 80 y hueco cuando no. Hace falta, porque entre UFBoot 70 y 94
# la mayoría de los nodos no llega a SH-aLRT 80, y un círculo lleno ahí sobredeclararía el soporte.
sup_pts = [(dep[id(c)], ycoord[id(c)], u, a) for c in tree.get_nonterminals()
           for u, a in [c_support(c)] if u is not None and u >= 50]

def circle(x, y, u, filled):
    f = min(max((u - 50) / 50, 0.0), 1.0)
    s = 4 + 26 * f ** 1.6
    if filled:
        axT.scatter([x], [y], s=s, color=(0.1, 0.1, 0.1, 0.25 + 0.75 * f),
                    linewidths=0, zorder=4)
    else:
        axT.scatter([x], [y], s=s, facecolor="white", linewidths=0.5,
                    edgecolor=(0.1, 0.1, 0.1, 0.35 + 0.65 * f), zorder=4)

for x, y, u, a in sup_pts:
    circle(x, y, u, a is None or a >= 80)

for l in leaves:
    axT.annotate(" " + l.name.replace("_", " "), (dep[id(l)], ycoord[id(l)]),
                 va="center", ha="left", fontsize=5.5, color=INK, style="italic")

xmaxT = max(dep[id(l)] for l in leaves)
axT.set_xlim(-xmaxT * 0.01, xmaxT * 1.62)
axT.set_ylim(-2, n + 1)
axT.set_yticks([])
axT.set_xlabel("substitutions per site")
axT.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
despine(axT, keep=("bottom",))
# leyenda de tamaños, tres clases de referencia
for i, v in enumerate((60, 80, 100)):
    yy = n - 2 - i * 4.5
    circle(xmaxT * 1.08, yy, v, True)
    axT.annotate(f"UFBoot {v}", (xmaxT * 1.125, yy), va="center", ha="left",
                 fontsize=6.5, color=MUTED)
yy = n - 2 - 3 * 4.5
circle(xmaxT * 1.08, yy, 100, False)
axT.annotate("SH-aLRT below 80", (xmaxT * 1.125, yy), va="center", ha="left",
             fontsize=6.5, color=MUTED)
figT.subplots_adjust(left=0.015, right=0.995, top=0.997, bottom=0.016)
save_esm(figT, 6, "The phylogeny of Fig. 1 with all 259 species named, and ultrafast bootstrap "
                  "and SH-aLRT support drawn on every internal branch.")
plt.close(figT)
print("Fig1 y el recurso en linea 4 listos")
print(f1[["order", "n", "mono", "ufboot", "shalrt"]].to_string(index=False))

# ================================================================ figura 2
bounds, k = {}, 0
for e in CORE:
    bounds[e] = (k, k + int(hx.loc[e, "aa"])); k += int(hx.loc[e, "aa"])
NAA = k
pos = sites.sort_values("human_pos")
x = pos.human_pos.values

fig = plt.figure(figsize=(W2, 126 * MM))
# el panel a necesita alto suficiente para que el número del exón quepa dentro de la barra
gs = fig.add_gridspec(3, 1, height_ratios=[0.34, 1, 1], hspace=0.16,
                      left=0.085, right=0.985, top=0.955, bottom=0.09)

# -- mapa de exones a escala
axM = fig.add_subplot(gs[0])
for e in CORE:
    s, t_ = bounds[e]
    axM.add_patch(plt.Rectangle((s + 1, 0.42), t_ - s, 0.52, zorder=3,
                                facecolor=FOCUS if e == "E10" else REF,
                                edgecolor="white", linewidth=0.5))
    if t_ - s >= 18:
        axM.annotate(e[1:], ((s + t_) / 2 + 1, 0.68), ha="center", va="center",
                     color="white", fontsize=8, zorder=4)
# Las cuatro repeticiones de unión a microtúbulos ocupan los residuos 244 a 368 (Trushina et al.
# 2019). El exón 9 empieza antes, en la región rica en prolina, y el 13 termina con la cola
# carboxilo terminal, así que la línea no cubre los exones enteros.
MBR0, MBR1 = 244, 368
axM.plot([MBR0, MBR1], [0.28, 0.28], color=MUTED, linewidth=0.8)
axM.annotate("microtubule binding repeats", ((MBR0 + MBR1) / 2, 0.16), ha="center", va="top",
             color=MUTED, fontsize=8)
axM.set_xlim(0, NAA + 1); axM.set_ylim(-0.08, 1.0)
axM.set_xticks([]); axM.set_yticks([])
for sp in axM.spines.values():
    sp.set_visible(False)
panel(axM, "a", dx=-0.075, dy=0.98)

# Techo en 4.5. Deja fuera 7 codones de 441 en vez de los 76 que recortaba un techo de 2, y
# mantiene escala lineal, que es lo que permite leer "dS alta en todas partes salvo el exón 10".
# Los dos paneles comparten escala para que el contraste entre dS y dN se lea directo.
TOPS = 4.5
for idx, (col, lab, letter) in enumerate([("fel_dS", "Synonymous rate (dS)", "b"),
                                          ("fel_dN", "Nonsynonymous rate (dN)", "c")]):
    axr = fig.add_subplot(gs[idx + 1])
    v = pos[col].values
    inE10 = (x >= 275) & (x <= 305)
    # 0.7 pt deja un hueco visible entre codones contiguos, que ocupan 1.12 pt cada uno
    axr.vlines(x[~inE10], 0, np.clip(v[~inE10], 0, TOPS), color=REF, linewidth=0.7, zorder=3)
    axr.vlines(x[inE10], 0, np.clip(v[inE10], 0, TOPS), color=FOCUS, linewidth=1.1, zorder=4)
    off = x[v > TOPS]
    if len(off):
        axr.plot(off, np.full(len(off), TOPS), marker="^", ms=2.2, ls="", color=REF,
                 markeredgewidth=0, clip_on=False, zorder=5)
    axr.set_xlim(0, NAA + 1)
    axr.set_ylim(0, TOPS * 1.05)
    axr.set_yticks([0, 1, 2, 3, 4])
    axr.set_ylabel(lab)
    axr.yaxis.grid(True, color=GRID, zorder=0)
    axr.set_axisbelow(True)
    for e in CORE:
        axr.axvline(bounds[e][0] + 0.5, color="#ebebeb", linewidth=0.4, zorder=1)
    despine(axr)
    panel(axr, letter, dx=-0.075, dy=1.01)
    nover = int((v > TOPS).sum())
    if nover:
        axr.annotate(f"{nover} codons above {TOPS:g}, marked with a caret", (6, TOPS * 1.04),
                     ha="left", va="top", color=MUTED)
    if letter == "c":
        axr.set_xlabel("Codon position in human 2N4R tau")
        axr.set_xticks([1, 50, 100, 150, 200, 250, 300, 350, 400, 441])
    else:
        axr.set_xticks([])

save(fig, 2)
plt.close(fig)

# ================================================================ figura 5
loc = pd.read_csv("results/hyrax/locus_summary.tsv", sep="\t")
hits = pd.concat([pd.read_csv(f"results/hyrax/{s}.exon_hits.tsv", sep="\t")
                  for s in loc.species])
aln = AlignIO.read("results/hyrax/hyrax_afro.aln.faa", "fasta")
NAMES = {"Procavia_capensis_model": r"$\it{P.\,capensis}$ (model)",
         "Heterohyrax_brucei_model": r"$\it{H.\,brucei}$ (model)",
         "Loxodonta_long": r"$\it{L.\,africana}$", "Trichechus_long": r"$\it{T.\,manatus}$",
         "Homo_776": r"$\it{H.\,sapiens}$ (776 aa)"}
ORDER5 = ["Procavia_capensis_model", "Heterohyrax_brucei_model", "Loxodonta_long",
          "Trichechus_long", "Homo_776"]
rec = {r.id: str(r.seq) for r in aln}
L = aln.get_alignment_length()

fig = plt.figure(figsize=(W2, 92 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.52,
                      left=0.175, right=0.985, top=0.93, bottom=0.11)

# -- panel a, exones recuperados en el locus
axA = fig.add_subplot(gs[0])
for i, sp in enumerate(loc.species):
    y = 1 - i
    r = loc[loc.species == sp].iloc[0]
    axA.plot([0, r.end - r.start], [y, y], color=MUTED, linewidth=0.8, zorder=2)
    sub = hits[hits.species == sp]
    for _, h in sub.iterrows():
        a, b = sorted((h.ss, h.se))
        axA.add_patch(plt.Rectangle((min(a, b) - r.start, y - 0.17),
                                    max(abs(b - a), 180), 0.34, zorder=3,
                                    facecolor=FOCUS if h.q == "human_E10" else REF,
                                    edgecolor="none"))
    g, s = sp.split("_")[:2]
    axA.annotate(f"$\\it{{{g[0]}.\\,{s}}}$\n{r.assembly}", (-600, y), ha="right", va="center",
                 color=INK, linespacing=1.4)
axA.set_xlim(-600, 38500)
axA.set_ylim(-0.75, 1.75)
axA.set_yticks([])
axA.set_xticks([0, 10000, 20000, 30000])
axA.set_xticklabels(["0", "10", "20", "30"])
axA.set_xlabel("kb along the scaffold", labelpad=2)
despine(axA, keep=("bottom",))
panel(axA, "a", dx=-0.175, dy=1.02)
axA.legend(handles=[
    Line2D([], [], color=FOCUS, lw=4, label="exon 10"),
    Line2D([], [], color=REF, lw=4, label="other exons found by tblastn"),
], loc="upper right", frameon=False, handletextpad=0.5, borderpad=0.1, labelspacing=0.35,
   bbox_to_anchor=(1.0, 1.12), ncol=2)

# -- panel b, mapa de huecos del alineamiento
axB = fig.add_subplot(gs[1])
for i, key in enumerate(ORDER5):
    y = len(ORDER5) - 1 - i
    s = rec[key]
    runs_, st = [], None
    for j, ch in enumerate(s):
        if ch != "-" and st is None:
            st = j
        elif ch == "-" and st is not None:
            runs_.append((st, j)); st = None
    if st is not None:
        runs_.append((st, len(s)))
    hyr = key.endswith("_model")
    for a, b in runs_:
        axB.add_patch(plt.Rectangle((a, y - 0.22), b - a, 0.44, zorder=3,
                                    facecolor=FOCUS if hyr else REF, edgecolor="none"))
    axB.annotate(NAMES[key], (-8, y), ha="right", va="center", color=INK)
axB.set_xlim(-8, L)
axB.set_ylim(-0.75, len(ORDER5) - 0.25)
axB.set_yticks([])
axB.set_xlabel("alignment column", labelpad=2)
despine(axB, keep=("bottom",))
panel(axB, "b", dx=-0.175, dy=1.02)
fig.savefig("manuscript/suplementarios/hyrax_figure.png", dpi=300)
plt.close(fig)
print("Fig2, ESM_6 y la figura de los damanes listas")
