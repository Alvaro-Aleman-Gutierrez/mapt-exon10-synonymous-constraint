"""Paso 22: ¿la limitación importa? Identidad sinónima de los exones alternativos.

La limitación declarada es que la mitad de la secuencia codificante queda sin probar. Para saber si
eso esconde una señal parecida a la del exón 10, se aplica a los quince exones la misma medida libre
de modelo del panel de amniotas, o sea la fracción de codones idénticos contados solo donde el
aminoácido no cambió. No necesita filogenia ni modelo de codones, así que tolera la divergencia de
los exones alternativos, que es donde FEL no serviría.

La super-isoforma humana es una concatenación de péptidos de exón, no un transcrito, así que no
tiene CDS propio. Por eso tanto la especie como el humano se mapean contra ella, cada uno con la
isoforma anotada que mejor cubre el exón, y la comparación corre sobre las posiciones que ambos
alcanzan. Eso hace la medida simétrica entre los dos lados.

Si los exones alternativos muestran identidad sinónima ordinaria, la limitación es formal y no
esconde ningún resultado. Si alguno mostrara exceso, sería un hallazgo que habría que reportar.
"""
import re
import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

ALT = ["E4aL", "E4a", "E6", "E8"]
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
REF = "Homo sapiens"

hx = pd.read_csv("data/human_exons.tsv", sep="\t")
bounds, k = {}, 0
for e in hx.exon:
    aa = int(hx.loc[hx.exon == e, "aa"].iloc[0]); bounds[e] = (k, k + aa); k += aa
SUP = "".join(hx.peptide)

al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                           substitution_matrix=substitution_matrices.load("BLOSUM62"))
al.target_end_gap_score = al.query_end_gap_score = 0.0

sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t"); sel = sel[sel.included]
order_of = dict(zip(sel.species, sel["order"]))
pres = pd.read_csv("results/exon_presence_all_isoforms.tsv", sep="\t")

def tag(desc):
    d = dict(re.findall(r"\[(\w+)=([^\]]+)\]", desc))
    return d.get("organism"), d.get("isoform") or d.get("transcript")

prot, cds = {}, {}
for r in SeqIO.parse("data/mammal_protein.faa", "fasta"):
    sp, iso = tag(r.description)
    prot[r.id.split()[0]] = (sp, iso, str(r.seq))
for r in SeqIO.parse("data/mammal_cds.fna", "fasta"):
    sp, iso = tag(r.description)
    cds[(sp, iso)] = str(r.seq)

def codons(s):
    return [s[i:i + 3] for i in range(0, len(s) - len(s) % 3, 3)]

def threaded(acc):
    """Devuelve (posicion en la super-isoforma -> (aminoacido, codon)) para una isoforma anotada."""
    rec = prot.get(acc)
    if not rec:
        return None
    sp, iso, aa = rec
    nt = cds.get((sp, iso))
    if not nt:
        return None
    sc = codons(nt)
    if len(sc) < len(aa):
        return None
    out = {}
    for (a0, a1), (b0, b1) in zip(*al.align(SUP, aa)[0].aligned):
        for o in range(a1 - a0):
            j = b0 + o
            if j < len(sc):
                out[a0 + o] = (aa[j], sc[j])
    return out

# isoforma humana de referencia por exón, con la misma regla de mayor cobertura
hum_pres = pres[pres.species == REF]
ref_map = {}
for e in ALT + CORE:
    sub = hum_pres[hum_pres[e] >= 0.5]
    if sub.empty:
        print(f"  aviso, sin isoforma humana anotada para {e}")
        continue
    acc = sub.loc[sub[e].idxmax(), "protein"]
    th = threaded(acc)
    if th:
        ref_map[e] = (acc, th)

rows = []
for e in ALT + CORE:
    if e not in ref_map:
        continue
    _, ref = ref_map[e]
    s, t = bounds[e]
    cand = pres[pres.species.isin(order_of) & (pres.species != REF)]
    best = cand.loc[cand.groupby("species")[e].idxmax()]
    best = best[best[e] >= (0.5 if e in ALT else 0.9)]
    for _, r in best.iterrows():
        th = threaded(r.protein)
        if not th:
            continue
        pos = [i for i in range(s, t) if i in th and i in ref]
        if len(pos) < 0.5 * (t - s):
            continue
        same = [i for i in pos if th[i][0] == ref[i][0]]
        if len(same) < 8:
            continue
        rows.append({
            "exon": e, "set": "alternative" if e in ALT else "core",
            "species": r.species, "order": order_of[r.species],
            "n_compared": len(pos), "n_same_aa": len(same),
            "aa_identity": round(len(same) / len(pos), 4),
            "nt_identity": round(sum(sum(1 for x, y in zip(th[i][1], ref[i][1]) if x == y)
                                     for i in pos) / (3 * len(pos)), 4),
            "syn_identity": round(sum(1 for i in same if th[i][1] == ref[i][1]) / len(same), 4),
        })

df = pd.DataFrame(rows)
from esm import write_esm
write_esm(df, 7, "exon_synonymous_identity")

summ = (df.groupby(["set", "exon"])
          .agg(species=("species", "nunique"), aa=("aa_identity", "mean"),
               nt=("nt_identity", "mean"), syn=("syn_identity", "mean"))
          .round(3))
order = [("alternative", e) for e in ALT] + [("core", e) for e in CORE]
summ = summ.reindex([i for i in order if i in summ.index])
summ["syn_minus_aa"] = (summ.syn - summ.aa).round(3)
print(summ.to_string())
print("\nfilas:", len(df), "| recurso en línea 7 escrito")

core = df[(df.set == "core") & (~df.exon.isin(["E10"]))].syn_identity
print(f"\nreferencia, exones núcleo sin E10: identidad sinónima media {core.mean():.3f}")
print(f"exón 10: {df[df.exon=='E10'].syn_identity.mean():.3f}")
for e in ALT:
    v = df[df.exon == e].syn_identity
    if len(v):
        print(f"{e}: {v.mean():.3f}  (n={len(v)})")

# ---- comparación pareada por distancia, solo Rodentia y Artiodactyla
# El filtro por identidad de aminoácido sesga hacia parientes cercanos justo en los exones
# divergentes, así que se compara dentro de los dos órdenes mejor muestreados y alejados del humano.
sub = df[df.order.isin(["Rodentia", "Artiodactyla"])]
mt = (sub.groupby("exon")
         .agg(species=("species", "nunique"), aa=("aa_identity", "mean"),
              syn=("syn_identity", "mean"))
         .round(3))
mt = mt[mt.species >= 12].sort_values("syn", ascending=False)
mt["set"] = ["alternative" if e in ALT else "core" for e in mt.index]
mt["spliced"] = ["alternative" if e in {"E2", "E3", "E4aL", "E4a", "E6", "E8", "E10"} else "constitutive"
                 for e in mt.index]
mt.reset_index()[["exon", "set", "spliced", "species", "aa", "syn"]] \
  .to_csv("results/tables/tableS3_matched_synonymous.tsv", sep="\t", index=False)
cap = ("Table S3. Synonymous identity of each exon within Rodentia and Artiodactyla, the two best "
       "sampled orders distant from human. Restricting the comparison to these orders removes the "
       "bias that a filter on amino acid identity would introduce, since such a filter retains only "
       "close relatives for the more divergent exons. Values come from the annotated isoform with "
       "the highest coverage of each exon, which for the alternative exons is not the 2N4R sequence "
       "the selection analyses used.")
with open("manuscript/suplementarios/tables/tableS3_matched_synonymous.md", "w") as f:
    f.write(f"**{cap}**\n\n")
    o = mt.reset_index()[["exon", "set", "spliced", "species", "aa", "syn"]]
    o.columns = ["Exon", "Set", "Splicing", "Species", "Amino acid identity", "Synonymous identity"]
    f.write(o.to_markdown(index=False, floatfmt=".3f") + "\n")
print("\n--- Tabla 7, Rodentia + Artiodactyla")
print(mt.to_string())
