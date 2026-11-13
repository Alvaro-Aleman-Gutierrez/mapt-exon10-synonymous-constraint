"""Paso 11: MAPT en genomas de Hyracoidea (Procavia capensis, Heterohyrax brucei).

1) tblastn de exones humanos y de isoformas largas de afroterios contra cada genoma: localiza el locus.
2) exonerate protein2genome sobre el locus (±30 kb): modelo génico con sitios de empalme.
3) Asigna exones del modelo a exones humanos y busca residuos extra entre E2-E3 y E3-E4
   (la "inserción A/B" del manuscrito original).
"""
import collections, itertools, os, subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

GENOMES = {"Procavia_capensis": "GCA_004026925.3", "Heterohyrax_brucei": "GCA_028571685.1"}
EXONERATE = "exonerate"
OUT = "results/hyrax"; os.makedirs(OUT, exist_ok=True)
DB = "data/genomes"

prot = {r.id: r for r in SeqIO.parse("data/mammal_protein.faa", "fasta")}
def longest(sp):
    return max((r for r in prot.values() if f"[organism={sp}]" in r.description), key=len)
queries = {"Homo_sapiens_776": prot["NP_001116538.2"], "Homo_sapiens_2N4R": prot["NP_005901.2"],
           "Loxodonta_africana_long": longest("Loxodonta africana"),
           "Trichechus_manatus_long": longest("Trichechus manatus latirostris")}
hx = pd.read_csv("data/human_exons.tsv", sep="\t")
with open(f"{OUT}/queries.faa", "w") as f:
    for k, r in queries.items():
        f.write(f">{k} {r.id}\n{r.seq}\n")
    for e, pep in zip(hx.exon, hx.peptide):
        f.write(f">human_{e}\n{pep}\n")

summary = []
for sp, acc in GENOMES.items():
    bl = f"{OUT}/{sp}.tblastn.tsv"
    if not os.path.exists(bl):
        subprocess.run(["tblastn", "-query", f"{OUT}/queries.faa", "-db", f"{DB}/{acc}", "-evalue", "1e-5",
                        "-seg", "no", "-num_threads", "10", "-max_target_seqs", "20",
                        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
                        "-out", bl], check=True)
    h = pd.read_csv(bl, sep="\t", header=None,
                    names="q s pid len qs qe ss se ev bits".split())
    # locus: secuencia con mayor bitscore acumulado para las isoformas completas
    best = h[h.q.str.contains("long|776")].groupby("s").bits.sum().sort_values(ascending=False)
    print(sp, "\n", best.head(5))
    scaf = best.index[0]
    loc = h[h.s == scaf]
    lo, hi = min(loc[["ss", "se"]].min()), max(loc[["ss", "se"]].max())
    region = f"{OUT}/{sp}.locus.fna"
    subprocess.run(["blastdbcmd", "-db", f"{DB}/{acc}", "-entry", scaf,
                    "-range", f"{max(1, lo - 30000)}-{hi + 30000}", "-out", region], check=True)
    # exones humanos: dónde caen y con qué identidad (mejor hit en el locus)
    ex = h[h.q.str.startswith("human_E") & (h.s == scaf)].sort_values("bits", ascending=False).groupby("q").head(1)
    ex.assign(species=sp).to_csv(f"{OUT}/{sp}.exon_hits.tsv", sep="\t", index=False)
    other = h[h.q.str.startswith("human_E") & (h.s != scaf)].q.unique()
    summary.append([sp, acc, scaf, lo, hi, ",".join(sorted(ex.q.str[6:])), ",".join(sorted(other))])

    for q in ["Loxodonta_africana_long", "Trichechus_manatus_long", "Homo_sapiens_776"]:
        qf = f"{OUT}/{q}.faa"; SeqIO.write([queries[q]], qf, "fasta")
        res = subprocess.run([EXONERATE, "--model", "protein2genome", "--query", qf, "--target", region,
                              "--showtargetgff", "yes", "--showalignment", "yes", "--bestn", "1",
                              "--ryo", ">%ti_%qi\n%tcs\n", "--refine", "region"],
                             capture_output=True, text=True, check=True).stdout
        open(f"{OUT}/{sp}.{q}.exonerate.txt", "w").write(res)

pd.DataFrame(summary, columns=["species", "assembly", "locus_seq", "start", "end",
                               "human_exons_found_in_locus", "human_exons_hit_elsewhere"]) \
  .to_csv(f"{OUT}/locus_summary.tsv", sep="\t", index=False)
