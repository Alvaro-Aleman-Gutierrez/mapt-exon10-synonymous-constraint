"""Paso 37: alineamientos de codones de los exones 9 a 13 dentro de cada clado de tetrápodos.

Por qué existe este paso. Fuera de mamíferos, la señal del exón 10 se midió como identidad por pares
con el humano. Para medirla con un modelo de codones dentro de cada linaje (paso 39), cada clado
necesita su propio alineamiento y su propio árbol, sin el humano.

El alineamiento se arma por bloques de exón. Cada especie usa la isoforma 4R del paso 15 y el mismo
mapa al humano, corregido por `e10copies` cuando el exón 10 está en tándem: el bloque E10 es la copia
5' y, si hay una segunda copia, entra como bloque aparte, E10b. Cada bloque se alinea con MAFFT
L-INS-i junto con el exón humano correspondiente, que solo sirve para marcar qué columnas son
homólogas a un residuo humano y luego se retira, y los bloques se concatenan. Así ninguna columna
mezcla exones distintos. Mamíferos usan las columnas de los exones 9 a 13 del alineamiento del
paso 6.

Filtros, por especie: el CDS traduce exactamente a la proteína, cada exón está cubierto en al
menos 0,9 (como en el paso 15), el tramo no tiene codones de parada internos y la identidad
proteica con el humano no cae más de 3 desviaciones estándar bajo la media de su grupo. Después
del alineamiento se retiran las especies con más de 10 codones sin homólogo humano en un bloque.
"""
import collections
import io
import json
import os
import subprocess
import tempfile

import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

import e10copies

OUT = "results/amniote"
os.makedirs(OUT, exist_ok=True)
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
REP = ["E9", "E10", "E11", "E12", "E13"]
BLOCKS = ["E9", "E10", "E10b", "E11", "E12", "E13"]
hx = pd.read_csv("data/human_exons.tsv", sep="\t").set_index("exon")
bounds, k = {}, 0
for e in CORE:
    bounds[e] = (k, k + int(hx.loc[e, "aa"])); k += int(hx.loc[e, "aa"])

HP = {r.id: str(r.seq) for r in SeqIO.parse("data/mapt_2N4R_protein.faa", "fasta")}["Homo_sapiens"]
human_block = {e: HP[slice(*bounds[e])] for e in REP}
human_block["E10b"] = human_block["E10"]

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                                substitution_matrix=substitution_matrices.load("BLOSUM62"))
aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0

picked = pd.read_csv("results/tetrapod_isoforms.tsv", sep="\t")

# CDS emparejado con su proteína por el reporte de productos de NCBI, como en el paso 15
tx2prot = {t["accession_version"]: t["protein"]["accession_version"]
           for rep in json.load(open("data/tetrapod_product_reports.json"))
           for t in rep["product"].get("transcripts", []) if "protein" in t}
P = {r.id.split(":")[0]: str(r.seq) for r in SeqIO.parse("data/tetrapod_protein.faa", "fasta")}
C = {}
for r in SeqIO.parse("data/tetrapod_cds.fna", "fasta"):
    p = tx2prot.get(r.id.split(":")[0])
    if p:
        C[p] = str(r.seq)

log, segs = [], collections.defaultdict(dict)
for _, row in picked.iterrows():
    sp, grp = row.species, row.group
    aa, nt = P[row.accession], C.get(row.accession, "")
    tr = str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip("*")
    if tr != aa.rstrip("*"):
        log.append((sp, grp, "CDS does not translate to the protein")); continue
    aln = aligner.align(HP, aa)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for off in range(a1 - a0):
            m[a0 + off] = b0 + off
    hits = e10copies.remap(m, aa, HP, bounds)
    cov = {e: sum(1 for i in range(*bounds[e]) if i in m) / (bounds[e][1] - bounds[e][0])
           for e in REP}
    if min(cov.values()) < 0.9:
        log.append((sp, grp, f"covers less than 0.9 of {min(cov, key=cov.get)}")); continue
    spans = {}
    for e in REP:
        pos = [m[i] for i in range(*bounds[e]) if i in m]
        spans[e] = (min(pos), max(pos) + 1)
    if len(hits) >= 2:
        spans["E10b"] = (hits[1], hits[1] + 31)
    order = [spans[e] for e in BLOCKS if e in spans]
    if any(a[1] > b[0] for a, b in zip(order, order[1:])):
        log.append((sp, grp, "exon blocks overlap in the species sequence")); continue
    blk = {e: (aa[a:b], nt[3 * a:3 * b]) for e, (a, b) in spans.items()}
    if any("*" in v[0] for v in blk.values()):
        log.append((sp, grp, "internal stop codon in exons 9 to 13")); continue
    ident = sum(1 for e in REP for i in range(*bounds[e]) if i in m and HP[i] == aa[m[i]]) / \
        (bounds["E13"][1] - bounds["E9"][0])
    segs[grp][sp.replace(" ", "_")] = (blk, ident, len(hits))

def mafft(records):
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
        for n, s in records:
            f.write(f">{n}\n{s}\n")
    res = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8",
                          "--quiet", f.name], capture_output=True, text=True, check=True)
    os.unlink(f.name)
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(io.StringIO(res.stdout), "fasta")}

def thread(aln_aa, nt):
    out, j = [], 0
    for c in aln_aa:
        if c == "-":
            out.append("---")
        else:
            out.append(nt[3 * j:3 * j + 3]); j += 1
    assert j * 3 == len(nt), "el CDS no se consumió completo"
    return out

CLADES = {"Aves": ["Aves"], "Lepidosauria": ["Lepidosauria"], "Testudines": ["Testudines"],
          "Archelosauria": ["Aves", "Crocodylia", "Testudines"], "Anura": ["Anura"],
          "Caudata": ["Caudata"]}
summary = []
for clade, groups in CLADES.items():
    seqs = {}
    for g in groups:
        ids = pd.Series({sp: v[1] for sp, v in segs[g].items()}, dtype=float)
        cut = ids.mean() - 3 * ids.std() if len(ids) >= 5 else -1
        for sp, v in segs[g].items():
            if v[1] < cut:
                if len(groups) == 1:
                    log.append((sp.replace("_", " "), g,
                                f"protein identity {v[1]:.3f} below mean minus 3 SD"))
                continue
            seqs[sp] = v
    cols, codon = [], {sp: [] for sp in seqs}
    for e in BLOCKS:
        have = [sp for sp in seqs if e in seqs[sp][0]]
        if len(have) < 2:
            continue
        A = mafft([("HUMAN_REF", human_block[e])] + [(sp, seqs[sp][0][e][0]) for sp in have])
        cod = {sp: thread(A[sp], seqs[sp][0][e][1]) for sp in have}
        for j in range(len(A["HUMAN_REF"])):
            if all(cod[sp][j] == "---" for sp in have):
                continue
            cols.append((e, A["HUMAN_REF"][j]))
            for sp in seqs:
                codon[sp].append(cod[sp][j] if sp in cod else "---")
    # inserciones de más de 10 codones sin homólogo humano dentro de un bloque: casi siempre un
    # modelo génico que incluye secuencia ajena al exón. Se retira la especie y las columnas vacías.
    nonhom = [j for j, (e, h) in enumerate(cols) if h == "-"]
    for sp in list(seqs):
        n_ins = sum(1 for j in nonhom if codon[sp][j] != "---")
        if n_ins > 10:
            seqs.pop(sp); codon.pop(sp)
            if len(groups) == 1:
                log.append((sp.replace("_", " "), groups[0],
                            f"insertion of {n_ins} codons without a human homologue"))
    keep = [j for j in range(len(cols)) if any(codon[sp][j] != "---" for sp in seqs)]
    cols = [cols[j] for j in keep]
    codon = {sp: [c[j] for j in keep] for sp, c in codon.items()}
    with open(f"{OUT}/{clade}.codon.fna", "w") as f:
        for sp in seqs:
            f.write(f">{sp}\n{''.join(codon[sp])}\n")
    pd.DataFrame({"site": range(1, len(cols) + 1), "exon": [c[0] for c in cols],
                  "human_residue": [c[1] for c in cols]}).to_csv(
        f"{OUT}/{clade}.sites.tsv", sep="\t", index=False)
    n2 = sum(1 for v in seqs.values() if v[2] >= 2)
    summary.append((clade, len(seqs), len(cols), n2))

# Mammalia: mismas columnas de los exones 9 a 13 en el alineamiento del paso 6
sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.exon.isin(REP)].sort_values("column")
col = [int(c) - 1 for c in sites.column]
mam = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}
haln = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.aln.faa", "fasta")}["Homo_sapiens"]
with open(f"{OUT}/Mammalia.codon.fna", "w") as f:
    for sp, s in mam.items():
        f.write(f">{sp}\n{''.join(s[3 * c:3 * c + 3] for c in col)}\n")
pd.DataFrame({"site": range(1, len(col) + 1), "exon": list(sites.exon),
              "human_residue": [haln[c] for c in col]}).to_csv(
    f"{OUT}/Mammalia.sites.tsv", sep="\t", index=False)
summary.append(("Mammalia", len(mam), len(col), 0))

pd.DataFrame(log, columns=["species", "group", "reason"]).to_csv(
    f"{OUT}/excluded.tsv", sep="\t", index=False)
for c, n, L, n2 in summary:
    print(f"{c:14s} {n:4d} especies  {L:4d} codones  {n2:4d} con exón 10 en tándem")
print(len(log), "especies excluidas, motivos en", f"{OUT}/excluded.tsv")
