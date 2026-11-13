"""Paso 9: análisis de selección con HyPhy (entorno conda 'hyphy') sobre el alineamiento de codones 2N4R.

- FEL y MEME en todo el conjunto.
- SLAC por orden (n >= 5): omega global MG94 por orden.
- RELAX: ramas de Primates, Catarrhini y Hominoidea (Test), cada uno contra el resto (Reference).
La topología es el árbol ML de proteína restringido por la taxonomía de NCBI (paso 7c), sin longitudes
de rama: los nodos de primates del árbol sin restricción tienen soporte muy bajo y sesgarían las pruebas por ramas.
"""
import io, os, subprocess, sys
import pandas as pd
from Bio import Phylo, SeqIO

HYPHY = "hyphy"
OUT = "results/selection"; os.makedirs(OUT, exist_ok=True)
ALN = "results/mapt_2N4R_final.codon.fna"
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]; sel["id"] = sel.species.str.replace(" ", "_")
order_of = dict(zip(sel.id, sel.order))

tree = Phylo.read("results/phylo_protein/mapt_2N4R_constrained.treefile", "newick")
for c in tree.find_clades():
    c.branch_length = None
    if not c.is_terminal():
        c.name = None; c.confidence = None

def newick(t):
    s = io.StringIO(); Phylo.write(t, s, "newick", plain=True); return s.getvalue().strip()

def run(method, aln, nwk, tag, extra=()):
    js = f"{OUT}/{tag}.{method}.json"
    if os.path.exists(js):
        return js
    tf = f"{OUT}/{tag}.nwk"; open(tf, "w").write(nwk + "\n")
    cmd = [HYPHY, "CPU=10", method.lower(), "--alignment", aln, "--tree", tf, "--output", js, *extra]
    with open(f"{OUT}/{tag}.{method}.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    return js

steps = sys.argv[1:] or ["fel", "slac", "relax", "meme"]
if "fel" in steps:
    run("FEL", ALN, newick(tree), "all")
if "slac" in steps:
    seqs = {r.id: r for r in SeqIO.parse(ALN, "fasta")}
    for o, g in sel.groupby("order"):
        if len(g) < 5:
            continue
        t = Phylo.read(io.StringIO(newick(tree)), "newick")
        for leaf in [l for l in t.get_terminals() if order_of[l.name] != o]:
            t.prune(leaf)
        sub = f"{OUT}/{o}.codon.fna"; SeqIO.write([seqs[i] for i in g.id], sub, "fasta")
        run("SLAC", sub, newick(t), o)
def relax(clade, members, tag):
    """RELAX con todas las ramas del clado (incluida la rama troncal) como Test y el resto como Reference."""
    t = Phylo.read(io.StringIO(newick(tree)), "newick")
    # el treefile no tiene raíz (cae junto a Homo); se enraíza en Monotremata para definir el clado
    t.root_with_outgroup([l for l in t.get_terminals() if order_of[l.name] == "Monotremata"])
    leaves = [l for l in t.get_terminals() if l.name in members]
    mrca = t.common_ancestor(leaves)
    assert {l.name for l in mrca.get_terminals()} == set(members), f"{clade} no monofilético"
    for c in mrca.find_clades():
        c.name = f"{c.name or ''}{{Test}}"
    run("RELAX", ALN, newick(t), tag, ["--test", "Test", "--models", "Minimal"])

if "relax" in steps:
    relax("Primates", [s for s in order_of if order_of[s] == "Primates"], "primates")
# subclados de primates (NCBI Taxonomy), para contrastar con los cambios de splicing de Recinos et al. 2024
clades = pd.read_csv("data/primate_clades.tsv", sep="\t")
for clade in ("Catarrhini", "Hominoidea"):
    if f"relax_{clade}" in steps:
        relax(clade, list(clades[clades[clade]].id), clade.lower())
if "meme" in steps:
    run("MEME", ALN, newick(tree), "all")
