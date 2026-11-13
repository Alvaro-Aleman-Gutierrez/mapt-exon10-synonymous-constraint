"""Nivel de identidad sinónima que cabría esperar si los sitios sinónimos fueran libres.

Lo usan el paso 16, para la columna de la tabla 4, y el paso 27, para los valores que cita el texto.
Vive aparte para que las dos salidas no puedan discrepar.

La idea es simple. En cada codón humano del exón se mira con qué frecuencia el clado usa ese mismo
codón entre los sinónimos del aminoácido que codifica. Esa frecuencia es la probabilidad de coincidir
por azar en ese sitio. Promediada sobre los codones del exón da la identidad esperada bajo
saturación. Un clado cuya identidad observada cae por debajo de esa expectativa no está conservando
nada: está midiendo sitios ya saturados.
"""
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id

TBL = unambiguous_dna_by_id[1].forward_table
_cache = {}


def _codon_usage(species, fasta):
    """Uso de codones del clado, tomado de sus propias transcripciones de MAPT."""
    use = {}
    for r in SeqIO.parse(fasta, "fasta"):
        if "[organism=" not in r.description:
            continue
        if r.description.split("[organism=")[1].split("]")[0] not in species:
            continue
        s = str(r.seq)
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3].upper()
            if c in TBL:
                use.setdefault(TBL[c], {})[c] = use.setdefault(TBL[c], {}).get(c, 0) + 1
    return use


def expected(species, fasta, exons, human_codons):
    """Identidad sinónima esperada por exón. `human_codons` mapea exón a lista de codones humanos."""
    key = (tuple(sorted(species)), fasta)
    if key not in _cache:
        _cache[key] = _codon_usage(species, fasta)
    use = _cache[key]
    out = {}
    for ex in exons:
        v = []
        for c in human_codons[ex]:
            c = c.upper()
            if c in TBL and TBL[c] in use:
                f = use[TBL[c]]
                v.append(f.get(c, 0) / sum(f.values()))
        out[ex] = float(np.mean(v)) if v else float("nan")
    return out


def human_codons_by_exon(sites_tsv="results/selection/sites_summary.tsv",
                         codon_fna="results/mapt_2N4R_final.codon.fna"):
    """Codones humanos de cada exón, en el orden del alineamiento."""
    s = pd.read_csv(sites_tsv, sep="\t")
    s = s[s.human_pos.notna()]
    hum = {r.id: str(r.seq) for r in SeqIO.parse(codon_fna, "fasta")}["Homo_sapiens"]
    return {ex: [hum[3 * (c - 1):3 * c] for c in g.column]
            for ex, g in s.groupby("exon")}
