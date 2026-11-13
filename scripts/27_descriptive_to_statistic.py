"""Paso 27: pone un número donde el manuscrito afirmaba un patrón de forma solo descriptiva.

Seis cantidades, cada una atada a una frase concreta del texto. Nada de esto cambia un resultado
anterior; lo que hace es medir lo que antes era un adjetivo. Tres sirven a la tesis central, que la
restricción del exón 10 está en el nucleótido y no en la proteína, y tres defienden decisiones de
alcance que un revisor va a cuestionar.

Los valores p de Wilcoxon por especie que este script imprime en los ejes 3 y en el perímetro
quedan como referencia descriptiva y no se citan en el manuscrito. La inferencia sobre unidades
independientes, órdenes y familias, está en el paso 32.

Los conteos por clado usan el mismo filtro que la tabla 4, que exige las cinco repeticiones
medidas en la misma especie, para que las dos no se contradigan.
"""
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from scipy import stats
import statsmodels.api as sm  # noqa: F401  (la dependencia ya está en el entorno de los pasos previos)

REP = ["E9", "E10", "E11", "E12", "E13"]
GROUPS = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Caudata", "Anura"]
TBL = unambiguous_dna_by_id[1].forward_table
out = []


def say(k, v):
    print(f"  {k:<58s} {v}")
    out.append({"quantity": k, "value": v})


print("=" * 80)
print("EJE 1. La inversión entre proteína y nucleótido")
print("=" * 80)
# La prueba va al nivel de clado, no de especie. Las especies de un mismo clado comparten historia,
# así que un binomial sobre 116 aves trataría como independientes observaciones que no lo son. Los
# clados sí son réplicas independientes. El conteo por especie se conserva, pero solo para mostrar
# que la inversión no aparece al promediar, y sobre él no se corre ninguna prueba.
tet = pd.read_csv("results/tetrapod_exon_identity.tsv", sep="\t")
inv_total = inv_n = 0
clade_hits = clade_n = 0
for g in GROUPS:
    d = tet[tet.group == g]
    aa = d.pivot(index="species", columns="exon", values="aa_identity").dropna()
    nt = d.pivot(index="species", columns="exon", values="nt_identity").dropna()
    sp = aa.index.intersection(nt.index)
    aa, nt = aa.loc[sp], nt.loc[sp]
    q = int(((aa.E10 - aa.E11 < 0) & (nt.E10 - nt.E11 > 0)).sum())
    daa, dnt = (aa.E10 - aa.E11).mean(), (nt.E10 - nt.E11).mean()
    inq = daa < 0 and dnt > 0
    say(f"{g}: media del clado en el cuadrante",
        f"{'sí' if inq else 'no'} (Δaa {daa:+.3f}, Δnt {dnt:+.3f}); "
        f"especies {q}/{len(sp)}, descriptivo")
    if g not in ("Caudata", "Anura"):
        clade_n += 1
        clade_hits += int(inq)
    if g not in ("Mammalia", "Caudata", "Anura"):
        inv_total += q
        inv_n += len(sp)
say("amniotas no mamíferos, conteo descriptivo por especie", f"{inv_total}/{inv_n}")
pc = stats.binomtest(clade_hits, clade_n, 0.25, alternative="greater").pvalue
say("clados de amniotas con la media en el cuadrante",
    f"{clade_hits} de {clade_n}, binomial contra 0,25 p = {pc:.4f}")

print()
print("=" * 80)
print("EJE 2. El nivel esperado de identidad sinónima una vez saturados los sitios")
print("=" * 80)
# Para cada codón humano del exón se toma la frecuencia de ese mismo codón entre los sinónimos
# que el clado usa de hecho. Es la identidad que cabría esperar si el sitio sinónimo fuera libre.
sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()]
hum = {r.id: str(r.seq) for r in
       SeqIO.parse("results/mapt_2N4R_final.codon.fna", "fasta")}["Homo_sapiens"]
# Los mamíferos vienen de su propio archivo de CDS; el resto del panel, del de tetrápodos.
cds = list(SeqIO.parse("data/tetrapod_cds.fna", "fasta"))
cds_mam = list(SeqIO.parse("data/mammal_cds.fna", "fasta"))
for g in GROUPS:
    sp = set(tet[tet.group == g].species)
    cds_g = cds_mam if g == "Mammalia" else cds
    use = {}
    for r in cds_g:
        if "[organism=" not in r.description:
            continue
        if r.description.split("[organism=")[1].split("]")[0] not in sp:
            continue
        s = str(r.seq)
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3].upper()
            if c in TBL:
                use.setdefault(TBL[c], {})[c] = use.setdefault(TBL[c], {}).get(c, 0) + 1
    exp = {}
    for ex in REP:
        v = []
        for col in sites[sites.exon == ex].column:
            c = hum[3 * (col - 1):3 * col].upper()
            if c in TBL and TBL[c] in use:
                f = use[TBL[c]]
                v.append(f.get(c, 0) / sum(f.values()))
        exp[ex] = float(np.mean(v))
    # Mismo filtro que la tabla 4: solo especies con las cinco repeticiones medidas, para que
    # los valores que cita el texto coincidan con los de la tabla.
    obs = (tet[tet.group == g]
           .pivot(index="species", columns="exon", values="syn_identity")
           .dropna()[REP].mean())
    say(f"{g}: E10 observada contra esperada bajo saturación",
        f"{obs['E10']:.3f} contra {exp['E10']:.3f}")
    say(f"{g}: otras repeticiones, observada contra esperada",
        f"{obs.reindex([e for e in REP if e != 'E10']).mean():.3f} contra "
        f"{np.mean([exp[e] for e in REP if e != 'E10']):.3f}")
    say(f"{g}: exones con la observada por debajo de la esperada",
        f"{sum(obs[e] < exp[e] for e in REP)}/5")

print()
print("=" * 80)
print("EJE 3. Exón 10 contra exón 11 a distancia filogenética igualada")
print("=" * 80)
d7 = pd.read_csv("results/tables/ESM_7_exon_synonymous_identity.tsv", sep="\t", comment="#")
mt = d7[d7.order.isin(["Rodentia", "Artiodactyla"])]
w = mt.pivot_table(index="species", columns="exon", values="syn_identity")
s = w[["E10", "E11"]].dropna()
st, p = stats.wilcoxon(s.E10, s.E11)
say("E10 contra E11, pareado en Rodentia y Artiodactyla",
    f"n = {len(s)}, medias {s.E10.mean():.3f} y {s.E11.mean():.3f}, Wilcoxon p = {p:.1e}")
s6 = w[["E6", "E10"]].dropna()
st, p6 = stats.wilcoxon(s6.E6, s6.E10)
say("E6 contra E10, pareado",
    f"n = {len(s6)}, diferencia {s6.E6.mean() - s6.E10.mean():+.3f}, Wilcoxon p = {p6:.1e}")

print()
print("=" * 80)
print("PERÍMETRO 1. El exón 6 contra exones que el análisis sí cubre")
print("=" * 80)
d6 = pd.read_csv("results/tables/ESM_6_alt_exon_divergence.tsv", sep="\t", comment="#")
wd = d6.pivot_table(index="species", columns="exon", values="identity_to_human")
for other in ["E1", "E5"]:
    s = wd[["E6", other]].dropna()
    st, p = stats.wilcoxon(s.E6, s[other])
    say(f"E6 contra {other}, pareado por especie",
        f"n = {len(s)}, medias {s.E6.mean():.3f} y {s[other].mean():.3f}, Wilcoxon p = {p:.1e}")

print()
print("=" * 80)
print("PERÍMETRO 2 y 3. Monofilia por orden y señal disponible dentro de Primates")
print("=" * 80)
sup = pd.read_csv("results/phylo_protein/order_support.tsv", sep="\t")
t = sup[sup.order != "Monotremata"]
say("órdenes monofiléticos, excluido el grupo externo",
    f"{int(t.mono.sum())} de {len(t)}")
aln = AlignIO.read("results/mapt_2N4R_final.aln.faa", "fasta")
lin = pd.read_csv("data/ncbi_lineages.tsv", sep="\t", header=None, names=["sp", "lineage"])
prim = set(lin[lin.lineage.str.contains("Primates", na=False)].sp)
idx = [i for i, r in enumerate(aln) if r.id in prim]


def pi_sites(rows):
    """Sitios informativos para parsimonia: dos estados o más con dos copias o más cada uno."""
    n = 0
    for c in range(aln.get_alignment_length()):
        col = [aln[i, c] for i in rows]
        cnt = pd.Series([x for x in col if x not in "-X"]).value_counts()
        if (cnt >= 2).sum() >= 2:
            n += 1
    return n


say("sitios informativos en el alineamiento completo", pi_sites(range(len(aln))))
say(f"sitios informativos solo entre los {len(idx)} primates", pi_sites(idx))

pd.DataFrame(out).to_csv("results/descriptive_to_statistic.tsv", sep="\t", index=False)
print("\nresults/descriptive_to_statistic.tsv escrito")

print()
print("=" * 80)
print("Tasa sinónima de FEL en el exón 10 contra las otras repeticiones, por codón")
print("=" * 80)
_s = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
_s = _s[_s.human_pos.notna()]
_a = _s[_s.exon == "E10"].fel_dS
_b = _s[_s.exon.isin(["E9", "E11", "E12", "E13"])].fel_dS
_p = stats.mannwhitneyu(_a, _b, alternative="two-sided").pvalue
for _e in ["E9", "E10", "E11", "E12", "E13"]:
    say(f"mediana de dS de FEL en {_e}", f"{_s[_s.exon == _e].fel_dS.median():.3g}")
say("razón entre la mediana de las otras cuatro y la de E10",
    f"{_b.median() / _a.median():,.0f}, {np.log10(_b.median() / _a.median()):.2f} órdenes")
say("Mann-Whitney de dos colas, E10 contra E9, E11, E12 y E13", f"p = {_p:.1e}")
pd.DataFrame(out).to_csv("results/descriptive_to_statistic.tsv", sep="\t", index=False)

print()
print("=" * 80)
print("Posiciones de mutaciones de FTDP-17 en el exón 10 y codón invariante")
print("=" * 80)
# Cinco posiciones donde la mutación descrita actúa sobre el splicing y dos donde actúa sobre la
# proteína. Con siete posiciones la prueba no puede ser concluyente y se reporta como anécdota.
_sp = [279, 280, 284, 296, 305]
_pr = [301, 303]
_s2 = sites[sites.human_pos.notna()].copy()
_s2["pos"] = _s2.human_pos.astype(int)
_z = lambda ps: int(_s2[_s2.pos.isin(ps)].codon_invariant.sum())
_tab = [[_z(_sp), len(_sp) - _z(_sp)], [_z(_pr), len(_pr) - _z(_pr)]]
_p = stats.fisher_exact(_tab)[1]
say("codones invariantes en posiciones de splicing / de proteína",
    f"{_z(_sp)} de {len(_sp)} / {_z(_pr)} de {len(_pr)}, Fisher p = {_p:.2f}")
pd.DataFrame(out).to_csv("results/descriptive_to_statistic.tsv", sep="\t", index=False)
