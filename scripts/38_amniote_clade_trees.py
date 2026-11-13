"""Paso 38: árboles de cada clado para el FUBAR del paso 39.

Igual que en mamíferos (pasos 7b y 28), la topología se restringe a los linajes de NCBI Taxonomy,
con politomías permitidas, y IQ-TREE resuelve los nodos que la taxonomía deja abiertos a partir del
alineamiento de nucleótidos del clado (GTR+F+G4). FUBAR recibe solo la topología y estima de nuevo
las longitudes de rama. También se guarda el árbol sin restringir, para la prueba de sensibilidad.

Los linajes se bajan de NCBI Taxonomy con E-utilities, sin datos de contacto, y quedan en
data/tetrapod_lineages.tsv para que el paso no dependa de la red al repetirse.
"""
import collections
import io
import os
import subprocess
import time
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd
from Bio import Phylo, SeqIO

OUT = "results/amniote"
LIN = "data/tetrapod_lineages.tsv"
IQTREE = "iqtree2"
CLADES = ["Aves", "Lepidosauria", "Testudines", "Archelosauria", "Anura", "Caudata"]
tax = pd.read_csv("data/tetrapod_orthologs.tsv", sep="\t", dtype=str)

if not os.path.exists(LIN):
    rows, ids = [], list(tax.taxid)
    for i in range(0, len(ids), 100):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml"
               "&id=" + ",".join(ids[i:i + 100]))
        for t in ET.fromstring(urllib.request.urlopen(url, timeout=120).read()).findall("Taxon"):
            lin = [x.findtext("ScientificName") for x in t.findall("LineageEx/Taxon")]
            rows.append((t.findtext("TaxId"), ";".join(lin)))
        time.sleep(0.5)
    pd.DataFrame(rows, columns=["taxid", "lineage"]).to_csv(LIN, sep="\t", index=False)
lin = pd.read_csv(LIN, sep="\t", dtype=str).set_index("taxid").lineage
lineage = {r.species.replace(" ", "_"): lin[r.taxid].split(";") for r in tax.itertuples()}

def constraint(names, root):
    children = collections.defaultdict(set)
    for sp in names:
        L = lineage[sp]
        L = L[L.index(root) + 1:] if root in L else L
        path = [tuple([root] + L[:i]) for i in range(len(L) + 1)]
        for a, b in zip(path, path[1:]):
            children[a].add(b)
        children[path[-1]].add(sp)

    def nw(node):
        kids = [nw(k) if isinstance(k, tuple) else k for k in sorted(children[node], key=str)]
        return kids[0] if len(kids) == 1 else "(" + ",".join(kids) + ")"
    return nw((root,)) + ";"

ROOT = {"Aves": "Aves", "Lepidosauria": "Lepidosauria", "Testudines": "Testudines",
        "Archelosauria": "Archelosauria", "Anura": "Anura", "Caudata": "Caudata"}
for c in CLADES:
    aln = f"{OUT}/{c}.codon.fna"
    names = [r.id for r in SeqIO.parse(aln, "fasta")]
    cons = f"{OUT}/{c}.constraint.nwk"
    open(cons, "w").write(constraint(names, ROOT[c]) + "\n")
    for tag, extra in (("constrained", ["-g", cons]), ("free", [])):
        pre = f"{OUT}/{c}.{tag}"
        if not os.path.exists(pre + ".treefile"):
            subprocess.run([IQTREE, "-s", aln, "-st", "DNA", "-m", "GTR+F+G4", "-T", "8",
                            "--prefix", pre, "-quiet", "-redo", *extra], check=True)
    tree = Phylo.read(f"{OUT}/{c}.constrained.treefile", "newick")
    for cl in tree.find_clades():
        cl.branch_length = None
        if not cl.is_terminal():
            cl.name = None
            cl.confidence = None
    s = io.StringIO()
    Phylo.write(tree, s, "newick", plain=True)
    open(f"{OUT}/{c}.fubar.nwk", "w").write(s.getvalue().strip() + "\n")
    print(f"{c:14s} {len(names):4d} especies, árbol listo")
