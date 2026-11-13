"""Paso 39: FUBAR dentro de cada clado, sobre los exones 9 a 13.

Mismos ajustes que el paso 28 (HyPhy 2.5.103, rejilla por defecto, topología sin longitudes de
rama). Mamíferos usan la topología del paso 28 recortada a los mismos exones; los otros clados, la
topología restringida del paso 38. Con el argumento `free` se corre además sobre los árboles sin
restringir, como prueba de sensibilidad.
"""
import os
import subprocess
import sys

HYPHY = "hyphy"
OUT = "results/amniote"
CLADES = ["Mammalia", "Aves", "Lepidosauria", "Testudines", "Archelosauria", "Anura", "Caudata"]
runs = [(c, "constrained") for c in CLADES]
if "free" in sys.argv:
    runs += [(c, "free") for c in CLADES if c != "Mammalia"]

for c, tag in runs:
    if c == "Mammalia":
        tree = "results/selection/fubar.nwk"
    elif tag == "constrained":
        tree = f"{OUT}/{c}.fubar.nwk"
    else:
        tree = f"{OUT}/{c}.free.treefile"
    js = f"{OUT}/{c}.{tag}.FUBAR.json"
    if os.path.exists(js):
        print(js, "ya existe"); continue
    cache = f"{OUT}/{c}.{tag}.FUBAR.cache"
    # HyPhy resuelve las rutas relativas desde su propia carpeta, así que van absolutas
    A = os.path.abspath
    cmd = [HYPHY, "CPU=10", "fubar", "--alignment", A(f"{OUT}/{c}.codon.fna"), "--tree", A(tree),
           "--output", A(js), "--cache", A(cache)]
    with open(f"{OUT}/{c}.{tag}.FUBAR.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(c, tag, "FUBAR terminado", flush=True)
