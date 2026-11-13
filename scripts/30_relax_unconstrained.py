"""Paso 30: RELAX sobre el árbol sin restringir, como análisis de sensibilidad.

Las tres corridas de RELAX del paso 9 usan la topología restringida por taxonomía, porque el árbol
génico no resuelve los nodos recientes de primates. Eso impone relaciones en lugar de inferirlas, y
un resultado nulo obtenido así podría deberse a la topología y no a la ausencia de un cambio de
régimen. La forma de saberlo es repetir las corridas sobre el árbol de máxima verosimilitud sin
restringir y comparar.

Hay una diferencia que conviene tener presente. En el árbol sin restringir, Catarrhini y Hominoidea
no son clados, así que no se puede etiquetar el ancestro común y sus descendientes. Lo que se
etiqueta entonces son las ramas terminales de las especies que componen cada grupo. Es un conjunto
de prueba legítimo para RELAX, que admite conjuntos arbitrarios de ramas, pero no es la misma
hipótesis que en el paso 9: ahí la prueba incluye las ramas internas del clado. La comparación dice
si el resultado nulo sobrevive, no si las dos corridas son equivalentes.
"""
import io
import os
import subprocess

import pandas as pd
from Bio import Phylo

HYPHY = "hyphy"
OUT = "results/selection"
ALN = "results/mapt_2N4R_final.codon.fna"
TREE = "results/phylo_protein/mapt_2N4R.treefile"

sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
sel["id"] = sel.species.str.replace(" ", "_")
order_of = dict(zip(sel.id, sel["order"]))
clades = pd.read_csv("data/primate_clades.tsv", sep="\t")

base = Phylo.read(TREE, "newick")
for c in base.find_clades():
    c.branch_length = None
    if not c.is_terminal():
        c.name = None
        c.confidence = None


def newick(t):
    s = io.StringIO()
    Phylo.write(t, s, "newick", plain=True)
    return s.getvalue().strip()


def run(tag, members):
    t = Phylo.read(io.StringIO(newick(base)), "newick")
    t.root_with_outgroup([l for l in t.get_terminals() if order_of.get(l.name) == "Monotremata"])
    leaves = [l for l in t.get_terminals() if l.name in members]
    mrca = t.common_ancestor(leaves)
    mono = {l.name for l in mrca.get_terminals()} == set(members)
    if mono:
        for c in mrca.find_clades():
            c.name = f"{c.name or ''}{{Test}}"
    else:
        for l in leaves:
            l.name = f"{l.name}{{Test}}"
    print(f"{tag}: {len(members)} especies, monofilético en el árbol sin restringir: {mono}")
    js = f"{OUT}/{tag}_unconstrained.RELAX.json"
    if os.path.exists(js):
        print(f"  {js} ya existe")
        return
    tf = f"{OUT}/{tag}_unconstrained.nwk"
    open(tf, "w").write(newick(t) + "\n")
    with open(f"{OUT}/{tag}_unconstrained.RELAX.log", "w") as log:
        subprocess.run([HYPHY, "CPU=10", "relax", "--alignment", ALN, "--tree", tf,
                        "--test", "Test", "--models", "Minimal", "--output", js],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"  terminado: {js}")


run("primates", [s for s in order_of if order_of[s] == "Primates"])
for clade in ("Catarrhini", "Hominoidea"):
    run(clade.lower(), list(clades[clades[clade]].id))
