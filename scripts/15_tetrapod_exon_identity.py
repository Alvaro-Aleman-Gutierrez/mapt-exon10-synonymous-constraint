"""Paso 15: identidad del exón 10 y de las repeticiones vecinas en tetrápodos no mamíferos.

Para cada especie se elige la isoforma con mayor cobertura del exón 10 humano, es decir la forma
4R, y se extrae el codón alineado a cada codón humano. El exón 10 se compara con E9, E11, E12 y
E13, que codifican el resto de la región de unión a microtúbulos y sirven de control interno.
"""
import collections, re
import e10copies
import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
REPEATS = ["E9", "E10", "E11", "E12", "E13"]
hx = pd.read_csv("data/human_exons.tsv", sep="\t").set_index("exon")
bounds, k = {}, 0
for e in CORE:
    bounds[e] = (k, k + int(hx.loc[e, "aa"])); k += int(hx.loc[e, "aa"])

prot = {r.id: str(r.seq) for r in SeqIO.parse("data/mapt_2N4R_protein.faa", "fasta")}
cds = {r.id: str(r.seq) for r in SeqIO.parse("data/mapt_2N4R_cds.fna", "fasta")}
HP, HC = prot["Homo_sapiens"], cds["Homo_sapiens"]

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                                substitution_matrix=substitution_matrices.load("BLOSUM62"))
aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0

tax = pd.read_csv("data/tetrapod_orthologs.tsv", sep="\t", dtype=str)
group_of = dict(zip(tax.species, tax.group))
order_of, family_of, gene_of = (dict(zip(tax.species, tax[c])) for c in ("order", "family", "gene_id"))
accessions = []

# isoformas por especie, con el accession de la proteína como clave. El CDS se empareja con su
# proteína por el reporte de productos de NCBI (transcrito XM_ -> proteína XP_), igual que en
# mamíferos (paso 5). La etiqueta [transcript=X1] del FASTA no sirve para esto, porque en algunas
# especies la variante de transcrito X1 codifica otra isoforma que la proteína X1.
import json
tx2prot = {t["accession_version"]: t["protein"]["accession_version"]
           for rep in json.load(open("data/tetrapod_product_reports.json"))
           for t in rep["product"].get("transcripts", []) if "protein" in t}

def load(path, cds=False):
    out = collections.defaultdict(dict)
    for r in SeqIO.parse(path, "fasta"):
        d = dict(re.findall(r"\[(\w+)=([^\]]+)\]", r.description))
        sp = d.get("organism")
        acc = r.id.split(":")[0]
        if sp in group_of:
            key = tx2prot.get(acc) if cds else acc
            if key:
                out[sp][key] = (acc, str(r.seq))
    return out

P = load("data/tetrapod_protein.faa")
C = load("data/tetrapod_cds.fna", cds=True)

def codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]

def measure(sp, grp, exon, pos, hcod, scod, haa, saa):
    """Identidad de aminoácido, de nucleótido y sinónima frente al humano en un exón.

    La identidad sinónima se calcula solo en los codones cuyo aminoácido no cambió, así que mide
    conservación del nucleótido una vez descontada la restricción sobre la proteína.
    """
    aa_id = sum(1 for i in pos if haa(i) == saa(i)) / len(pos)
    nt_id = sum(sum(1 for x, y in zip(hcod(i), scod(i)) if x == y) for i in pos) / (3 * len(pos))
    syn = [i for i in pos if haa(i) == saa(i)]
    syn_id = sum(1 for i in syn if hcod(i) == scod(i)) / len(syn) if syn else float("nan")
    same = [i for i in pos if hcod(i) == scod(i)]
    return (sp, grp, exon, len(pos), round(aa_id, 4), round(nt_id, 4), len(syn),
            round(syn_id, 4), len(same), ";".join(str(i + 1) for i in same))

rows, picked = [], []
for sp in sorted(P):
    best = None
    for lab, (acc, aa) in P[sp].items():
        if lab not in C[sp] or len(aa) < 200:
            continue
        aln = aligner.align(HP, aa)[0]
        # mapa posicion humana -> posicion de la especie
        m = {}
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for off in range(a1 - a0):
                m[a0 + off] = b0 + off
        # con copias en tándem del exón 10 (aves, algunos escamosos), el exón 10 humano se
        # compara con la copia 5' completa y el exón 11 con el que sigue a la última copia, no con la
        # quimera que deja el alineamiento global
        hits = e10copies.remap(m, aa, HP, bounds)
        ncopies = max(len(hits), 1)
        s, e = bounds["E10"]
        cov = sum(1 for i in range(s, e) if i in m) / (e - s)
        ident = sum(1 for i in m if HP[i] == aa[m[i]]) / len(HP)
        if best is None or (cov, ident) > (best[0], best[1]):
            best = (cov, ident, acc, aa, m, lab, ncopies)
    # registro de cada especie, también de las descartadas, para listar sus accesiones en el
    # suplementario: el manuscrito promete las accesiones de todo lo que se analizó
    fam = family_of[sp]
    if best is None:
        accessions.append([sp, group_of[sp], order_of[sp], fam, gene_of[sp], "", "", "", "", False,
                           "no isoform of at least 200 residues paired to a coding sequence"])
        continue
    accessions.append([sp, group_of[sp], order_of[sp], fam, gene_of[sp], best[2],
                       C[sp][best[5]][0], len(best[3]), round(best[0], 2), best[0] >= 0.9,
                       "" if best[0] >= 0.9 else
                       "best isoform covers less than 0.9 of exon 10"])
    if best[0] < 0.9:
        continue
    cov, ident, acc, aa, m, lab, ncopies = best
    sc = codons(C[sp][lab][1]); hc = codons(HC)
    picked.append((sp, group_of[sp], acc, len(aa), round(cov, 2), round(ident, 3), ncopies))
    for e in REPEATS:
        s, t = bounds[e]
        pos = [i for i in range(s, t) if i in m and m[i] < len(sc)]
        if len(pos) < 0.9 * (t - s):
            continue
        rows.append(measure(sp, group_of[sp], e, pos, lambda i: hc[i],
                            lambda i: sc[m[i]], lambda i: HP[i], lambda i: aa[m[i]]))

pd.DataFrame(picked, columns=["species", "group", "accession", "length", "E10_coverage",
                              "identity_to_human", "E10_copies"]).to_csv(
    "results/tetrapod_isoforms.tsv", sep="\t", index=False)
pd.DataFrame(accessions, columns=["species", "clade", "order", "family", "gene_id", "protein",
                                  "transcript", "length", "E10_coverage", "included",
                                  "exclusion_reason"]).to_csv(
    "results/tetrapod_accessions.tsv", sep="\t", index=False)
COLS = ["species", "group", "exon", "n_codons", "aa_identity", "nt_identity", "n_synonymous",
        "syn_identity", "identical_codons", "identical_positions"]
print(len(picked), "especies con isoforma 4R utilizable")

# Mammalia como referencia, medida igual sobre el alineamiento de codones del paso 6
sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()]
acod = {r.id: codons(str(r.seq)) for r in SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}
aprot = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_final.aln.faa", "fasta")}
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
sel = sel[sel.included]
for sp in acod:
    if sp == "Homo_sapiens":
        continue
    hc, sc, sa, ha = acod["Homo_sapiens"], acod[sp], aprot[sp], aprot["Homo_sapiens"]
    for e in REPEATS:
        col = [int(c) - 1 for c in sites[sites.exon == e].column]
        pos = [c for c in col if sc[c] != "---" and sa[c] != "-"]
        if len(pos) < 0.9 * len(col):
            continue
        rows.append(measure(sp.replace("_", " "), "Mammalia", e, pos,
                            lambda i: hc[i], lambda i: sc[i], lambda i: ha[i],
                            lambda i: sa[i]))

df = pd.DataFrame(rows, columns=COLS)
df.to_csv("results/tetrapod_exon_identity.tsv", sep="\t", index=False)

ORDER = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Caudata", "Anura"]
piv = df.pivot_table(index="group", columns="exon",
                     values=["aa_identity", "nt_identity", "syn_identity"]).round(3)
print("\n", piv.reindex(ORDER).to_string(), sep="")

# prueba apareada: identidad sinónima del exón 10 frente a la media de las otras repeticiones
from scipy.stats import wilcoxon
print("\nidentidad sinónima, E10 vs media(E9,E11,E12,E13), por especie:")
out = []
for grp in ORDER:
    sub = df[df.group == grp].pivot(index="species", columns="exon", values="syn_identity")
    sub = sub.dropna()
    if len(sub) < 4:
        continue
    other = sub[[c for c in REPEATS if c != "E10"]].mean(axis=1)
    d = sub.E10 - other
    w = wilcoxon(sub.E10, other)
    out.append((grp, len(sub), round(sub.E10.mean(), 3), round(other.mean(), 3),
                round(d.mean(), 3), int((d > 0).sum()), f"{w.pvalue:.2e}"))
print(pd.DataFrame(out, columns=["group", "n", "E10", "otras", "delta", "n_E10_mayor",
                                 "wilcoxon_p"]).to_string(index=False))
