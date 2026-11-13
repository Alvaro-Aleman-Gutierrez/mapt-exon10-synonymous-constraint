"""Paso 33: cruce de los elementos de Lin et al. (2011) con los exones humanos de MAPT.

Los elementos de restricción sinónima de Lin et al. vienen en coordenadas hg18, y los exones de este
trabajo en GRCh38, así que primero se pasan de un ensamblaje a otro con la cadena oficial de UCSC
(hg18ToHg38), que se lee aquí directamente para no depender de otra instalación. Un bloque se
convierte solo si cae entero dentro de un tramo alineado de la cadena. Luego se cuentan los
nucleótidos de cada exón cubiertos por elementos de cada resolución, 9, 15 y 30 codones.

Para el valor mínimo de lambda_s se lleva al mismo ensamblaje la pista de la tasa sinónima local.
"""
import gzip

import pandas as pd

CHAIN = "data/lin2011_sce/hg18ToHg38.over.chain.gz"


def load_chain(chrom="chr17"):
    """Tramos alineados (inicio origen, inicio destino, largo) de las cadenas de un cromosoma."""
    blocks = []
    with gzip.open(CHAIN, "rt") as f:
        cur = None
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "chain":
                cur = (p[2] == chrom and p[7] == chrom and p[4] == "+" and p[9] == "+")
                s, t = int(p[5]), int(p[10])
                continue
            if not cur:
                continue
            size = int(p[0])
            blocks.append((s, t, size))
            s += size + (int(p[1]) if len(p) > 1 else 0)
            t += size + (int(p[2]) if len(p) > 2 else 0)
    return sorted(blocks)


BL = load_chain()


def lift(a, b):
    """Convierte [a, b) si cabe entero en un tramo alineado, si no devuelve None."""
    for s, t, n in BL:
        if s <= a and b <= s + n:
            return t + (a - s), t + (b - s)
    return None


def bed12(path):
    """Bloques de un archivo BED12, en coordenadas de origen."""
    out = []
    for line in open(path):
        p = line.split("\t")
        start = int(p[1])
        sizes = [int(x) for x in p[10].strip(",").split(",")]
        offs = [int(x) for x in p[11].strip().strip(",").split(",")]
        out += [(start + o, start + o + s) for o, s in zip(offs, sizes)]
    return out


ex = pd.read_csv("data/human_exons.tsv", sep="\t")
ex["lo"], ex["hi"] = ex.start - 1, ex.end          # a medio abierto
rows = []
for res in (9, 15, 30):
    lifted, failed = [], 0
    for a, b in bed12(f"data/lin2011_sce/SCE{res}.MAPT.hg18.bed"):
        r = lift(a, b)
        if r is None:
            failed += 1
        else:
            lifted.append(r)
    for _, e in ex.iterrows():
        cov = sum(max(0, min(b, e.hi) - max(a, e.lo)) for a, b in lifted)
        rows.append({"resolution_codons": res, "exon": e.exon, "exon_nt": int(e.hi - e.lo),
                     "covered_nt": int(cov)})
    print(f"SCE{res}: {len(lifted)} bloques convertidos, {failed} sin convertir")

res = pd.DataFrame(rows)
res.to_csv("results/lin2011_overlap.tsv", sep="\t", index=False)
hit = res[res.covered_nt > 0]
print(hit.to_string(index=False))

# lambda_s mínimo dentro del exón 10 y mediana del gen
bg = pd.read_csv("data/lin2011_sce/lambda_s.MAPT.hg18.bedGraph", sep="\t", header=None,
                 names=["c", "s", "e", "v"])
bg["lifted"] = [lift(a, b) for a, b in zip(bg.s, bg.e)]
bg = bg.dropna(subset=["lifted"])
bg["s38"] = [x[0] for x in bg.lifted]
bg["e38"] = [x[1] for x in bg.lifted]
e10 = ex[ex.exon == "E10"].iloc[0]
inside = bg[(bg.e38 > e10.lo) & (bg.s38 < e10.hi)]
print(f"\nlambda_s: {len(bg)} ventanas convertidas en el gen, mediana {bg.v.median():.3f}; "
      f"{len(inside)} ventanas dentro del exón 10, mínimo {inside.v.min():.4f}")

# lambda_s mínimo por exón, para comparar el exón 10 con el resto
print("\nlambda_s mínimo por exón (ventanas que lo tocan)")
for _, e in ex.iterrows():
    w = bg[(bg.e38 > e.lo) & (bg.s38 < e.hi)]
    if len(w):
        print(f"  {e.exon:5s} {len(w):2d} ventanas, mínimo {w.v.min():.4f}, mediana {w.v.median():.3f}")
