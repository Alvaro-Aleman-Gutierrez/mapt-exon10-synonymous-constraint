"""Paso 12: traduce los modelos de exonerate (molde: elefante) y los alinea con afroterios y humano,
anotando cada columna con el exón humano, para evaluar inserciones entre E2-E3 y E3-E4."""
import re, subprocess, itertools
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

OUT = "results/hyrax"
recs = []
for sp in ["Procavia_capensis", "Heterohyrax_brucei"]:
    txt = open(f"{OUT}/{sp}.Loxodonta_africana_long.exonerate.txt").read()
    cds = "".join(txt.split("\n>", 1)[1].split("\n")[1:]).split("-- completed")[0].replace("\n", "")
    cds = re.sub("[^ACGTNacgtn]", "", cds)
    aa = str(Seq(cds[:len(cds) // 3 * 3]).translate())
    exons = [l.split("\t") for l in txt.splitlines() if "\texon\t" in l]
    lens = [int(e[4]) - int(e[3]) + 1 for e in exons]
    print(sp, len(cds), "nt;", len(aa), "aa; stops internos:", aa.rstrip("*").count("*"), "; exones (nt):", lens)
    recs.append(f">{sp}_model\n{aa.rstrip('*')}\n")
prot = {r.id: str(r.seq) for r in SeqIO.parse("data/mammal_protein.faa", "fasta")}
for name, pid in [("Loxodonta_long", "XP_064127056.1"), ("Trichechus_long", "XP_004374308.1"),
                  ("Homo_776", "NP_001116538.2")]:
    recs.append(f">{name}\n{prot[pid]}\n")
open(f"{OUT}/hyrax_afro.faa", "w").write("".join(recs))
aln = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", f"{OUT}/hyrax_afro.faa"],
                     capture_output=True, text=True, check=True).stdout
open(f"{OUT}/hyrax_afro.aln.faa", "w").write(aln)
seqs = {r.id: str(r.seq) for r in SeqIO.parse(f"{OUT}/hyrax_afro.aln.faa", "fasta")}

# exón humano por columna (isoforma 776 = todos los exones salvo E8)
hx = pd.read_csv("data/human_exons.tsv", sep="\t").set_index("exon")
order = ["E1", "E2", "E3", "E4", "E4aL", "E4a", "E5", "E6", "E7", "E9", "E10", "E11", "E12", "E13"]
lab = list(itertools.chain(*[[e] * hx.loc[e, "aa"] for e in order]))
col, k = [], 0
for c in seqs["Homo_776"]:
    col.append(lab[k] if c != "-" and k < len(lab) else (lab[k - 1] if k else "E1") + ("" if c != "-" else "+ins"))
    k += c != "-"
# imprime el bloque E1-E4aL con la etiqueta de exón
idx = [i for i, e in enumerate(col) if e.split("+")[0] in ("E1", "E2", "E3", "E4")]
lo, hi = idx[0], idx[-1] + 8
for n, s in seqs.items():
    print(f"{n:22s} {s[lo:hi]}")
print(f"{'exón (humano)':22s} " + "".join((e[1] if "+" not in e else "^") if e.split('+')[0] != "E4aL" else "L" for e in col[lo:hi]))
# inserciones de hyracoideos respecto a elefante y manatí
for n in ["Procavia_capensis_model", "Heterohyrax_brucei_model"]:
    ins = [(i, col[i]) for i in range(len(col)) if seqs[n][i] != "-" and seqs["Loxodonta_long"][i] == "-"
           and seqs["Trichechus_long"][i] == "-" and seqs["Homo_776"][i] == "-"]
    runs = [(g[0][1], len(g)) for _, g in itertools.groupby(ins, key=lambda x: x[1])]
    print(n, "residuos exclusivos de hyracoideos (exón humano precedente, n):", runs)
