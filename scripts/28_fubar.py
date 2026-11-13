"""Paso 28: FUBAR, para estimar la tasa sinónima por sitio dentro de un modelo de codones.

Por qué existe este paso. El argumento central del paper descansa en invariancia nucleotídica, en
identidades sinónimas y en las tasas puntuales que devuelve FEL. FEL estima alfa por sitio, pero su
prueba es alfa contra beta en ese mismo sitio, así que nunca pregunta si alfa está por debajo del
fondo del gen, que es justamente la afirmación del paper. FUBAR ajusta una rejilla de pares
(alfa, beta) sobre todos los sitios y devuelve la posterior de cada sitio sobre esa rejilla, de modo
que sí permite preguntar qué probabilidad hay de que la tasa sinónima de un codón esté por debajo de
la media del gen.

Misma topología que el paso 9, el árbol de proteína restringido por taxonomía y sin longitudes de
rama, para que los resultados sean comparables con FEL y MEME.
"""
import io
import os
import subprocess

from Bio import Phylo

HYPHY = "hyphy"
OUT = "results/selection"
ALN = "results/mapt_2N4R_final.codon.fna"
JS = f"{OUT}/all.FUBAR.json"

tree = Phylo.read("results/phylo_protein/mapt_2N4R_constrained.treefile", "newick")
for c in tree.find_clades():
    c.branch_length = None
    if not c.is_terminal():
        c.name = None
        c.confidence = None
s = io.StringIO()
Phylo.write(tree, s, "newick", plain=True)
tf = f"{OUT}/fubar.nwk"
open(tf, "w").write(s.getvalue().strip() + "\n")

if os.path.exists(JS):
    print(f"{JS} ya existe, no se vuelve a correr")
else:
    cmd = [HYPHY, "CPU=10", "fubar", "--alignment", ALN, "--tree", tf, "--output", JS]
    print(" ".join(cmd))
    with open(f"{OUT}/all.FUBAR.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    print("FUBAR terminado")
