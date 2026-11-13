"""Paso 7b: árbol de restricción a partir de los linajes de NCBI Taxonomy (todos los rangos con nombre).

Se usa con IQ-TREE (-g) para obtener un árbol ML consistente con los clados taxonómicos, que es la
topología de los análisis de selección (paso 9). El árbol del paso 7 (sin restricción) se reporta como
filogenia del gen.
"""
import collections
import pandas as pd
from Bio import Entrez

Entrez.email = "your.email@example.org"
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
tax = pd.read_csv("data/mammal_orthologs.tsv", sep="\t", dtype=str).set_index("species")
ids = {tax.loc[s, "taxid"]: s.replace(" ", "_") for s in sel.species}

lineage = {}
for t in Entrez.read(Entrez.efetch(db="taxonomy", id=",".join(ids), retmode="xml")):
    lin = [x["ScientificName"] for x in t["LineageEx"]]
    lineage[ids[t["TaxId"]]] = lin[lin.index("Mammalia") + 1:]

# árbol anidado: cada taxón con más de un descendiente se vuelve un nodo (politomías permitidas).
# Los nodos se identifican por la ruta completa porque hay nombres repetidos entre rangos
# (p. ej. género y subgénero Mustela).
children = collections.defaultdict(set)
for sp, lin in lineage.items():
    path = [tuple(["Mammalia"] + lin[:i]) for i in range(len(lin) + 1)]
    for a, b in zip(path, path[1:]):
        children[a].add(b)
    children[path[-1]].add(sp)

def newick(node):
    kids = [newick(k) if isinstance(k, tuple) else k for k in sorted(children[node], key=str)]
    return kids[0] if len(kids) == 1 else "(" + ",".join(kids) + ")"

open("results/phylo_protein/taxonomy_constraint.nwk", "w").write(newick(("Mammalia",)) + ";\n")
pd.DataFrame([(s, ";".join(l)) for s, l in lineage.items()], columns=["id", "lineage"]) \
  .to_csv("data/ncbi_lineages.tsv", sep="\t", index=False)
print(len(lineage), "especies en la restricción")
