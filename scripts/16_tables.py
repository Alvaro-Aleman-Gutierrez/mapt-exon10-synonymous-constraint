"""Paso 16: tablas 1 a 5 del manuscrito y los recursos en línea.

Las tablas se escriben en markdown (manuscript/tables/) y en TSV (results/tables/) desde las
salidas de los pasos 5, 6, 8, 13 y 15. Los recursos en línea quedan solo en TSV.
"""
import os
import pandas as pd

MD, TSV = "manuscript/tables", "results/tables"
for d in (MD, TSV):
    os.makedirs(d, exist_ok=True)
CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]

def fmt(x):
    """Tasas por sitio. Los valores por debajo de 0.001 se reportan como <0.001."""
    return "<0.001" if x < 1e-3 else f"{x:.3f}"

def write(name, df, caption, floatfmt=".3f", raw=False, md=MD):
    """raw evita que tabulate reinterprete como número las columnas ya formateadas (p. ej. 2.2e-20)."""
    df.to_csv(f"{TSV}/{name}.tsv", sep="\t", index=False)
    with open(f"{md}/{name}.md", "w") as f:
        f.write(f"**{caption}**\n\n")
        f.write(df.to_markdown(index=False, floatfmt=floatfmt,
                               disable_numparse=raw) + "\n")
    print(f"{name}  {len(df)} filas")

# ---- Tabla 1. Conservación por exón
hx = pd.read_csv("data/human_exons.tsv", sep="\t").set_index("exon")
idh = pd.read_csv("results/exon_identity_to_human.tsv", sep="\t")
idh = idh[idh.species != "Homo_sapiens"]
wo = pd.read_csv("results/exon_identity_within_order.tsv", sep="\t")
wo = wo[wo.n >= 5].groupby("exon").mean_pairwise_identity.mean()
t1 = pd.DataFrame({
    "Exon": CORE,
    "Residues": [int(hx.loc[e, "aa"]) for e in CORE],
    "Mean identity to human": [idh[idh.exon == e].identity_to_human.mean() for e in CORE],
    "Minimum": [idh[idh.exon == e].identity_to_human.min() for e in CORE],
    "Species identical to human": [int((idh[idh.exon == e].identity_to_human == 1).sum()) for e in CORE],
    "Mean identity within orders": [wo.get(e) for e in CORE],
})
write("table1_exon_conservation", t1,
      "Table 1. Amino acid conservation of the core exons of MAPT across 258 non-human mammals. "
      "Identity to human is computed over alignment columns where both sequences carry a residue. "
      "The last column averages mean pairwise identity within the eight orders sampled with at "
      "least five species.", floatfmt=".4f")

# ---- Tabla 2. Selección por sitio y por exón
fel = pd.read_csv("results/selection/fel_by_exon.tsv", sep="\t").set_index("exon")
sites = pd.read_csv("results/selection/sites_summary.tsv", sep="\t")
sites = sites[sites.human_pos.notna()]
med = sites.groupby("exon")[["fel_dS", "fel_dN"]].median()
t2 = pd.DataFrame({
    "Exon": CORE,
    "Codons": [int(fel.loc[e, "n_columns"]) for e in CORE],
    "Purifying (p<=0.05)": [int(fel.loc[e, "purifying"]) for e in CORE],
    "Purifying (q<=0.05)": [int(fel.loc[e, "purifying_q"]) for e in CORE],
    "Fraction": [fel.loc[e, "frac_purifying"] for e in CORE],
    "Diversifying": [int(fel.loc[e, "diversifying"]) for e in CORE],
    "MEME (p<=0.05)": [int(fel.loc[e, "meme"]) for e in CORE],
    "Median dS": [fmt(med.loc[e, "fel_dS"]) for e in CORE],
    "Median dN": [fmt(med.loc[e, "fel_dN"]) for e in CORE],
    "Invariant codons": [int(fel.loc[e, "codon_invariant"]) for e in CORE],
})
write("table2_site_selection", t2,
      "Table 2. Site level selection on MAPT by exon, from FEL and MEME on the codon alignment of "
      "259 mammals and the taxonomy constrained topology. Rates are medians of the per site "
      "maximum likelihood estimates. A codon counts as invariant when its three nucleotides are "
      "identical in every species that carries it.")

# ---- Tabla 3. Selección por linaje
slac = pd.read_csv("results/selection/slac_omega_by_order.tsv", sep="\t")
slac.columns = ["Order", "Species", "omega"]
write("tableS4_slac_by_order", slac,
      "Table S4. Global ratio of nonsynonymous to synonymous rates (omega) estimated by SLAC under "
      "MG94xREV, separately for each mammalian order with at least five species.",
      md="manuscript/suplementarios/tables")

relax = pd.read_csv("results/selection/relax_summary.tsv", sep="\t")
cl = pd.read_csv("data/primate_clades.tsv", sep="\t")
n = {"primates": 33, "catarrhini": int(cl.Catarrhini.sum()), "hominoidea": int(cl.Hominoidea.sum())}
# Sensibilidad a la topología: las mismas pruebas sobre el árbol ML sin restringir (paso 30).
import json
unc = {}
for tag in ("primates", "catarrhini", "hominoidea"):
    f = f"results/selection/{tag}_unconstrained.RELAX.json"
    if os.path.exists(f):
        tr = json.load(open(f))["test results"]
        unc[tag] = (tr["relaxation or intensification parameter"], tr["p-value"])
relax = pd.DataFrame({
    "Test clade": ["Primates", "Catarrhini", "Hominoidea"],
    "Species": [n[t] for t in relax.test_clade],
    "k": relax.k.values, "LRT": relax.LRT.values, "p": relax.p.values,
    "k unconstrained": [f"{unc[t][0]:.3f}" if t in unc else "—" for t in relax.test_clade],
    "p unconstrained": [f"{unc[t][1]:.3f}" if t in unc else "—" for t in relax.test_clade],
})
write("table6_relax", relax,
      "Table 6. RELAX tests of a shift in selection intensity. All branches of the test clade, "
      "the stem included, were compared against the rest of the tree. Values of k above 1 indicate "
      "intensified selection and values below 1 relaxed selection. The last two columns repeat each "
      "test on the unconstrained maximum likelihood tree instead of the taxonomy constrained one. "
      "Hominoidea is not a clade in that tree, so there its test set is the terminal branches of "
      "its nine species rather than a clade and its stem.")

# ---- Tabla 4. Panel de tetrápodos
d = pd.read_csv("results/tetrapod_exon_identity.tsv", sep="\t")
REP = ["E9", "E10", "E11", "E12", "E13"]
ORDER = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Caudata", "Anura"]
from scipy.stats import wilcoxon
import satexp
HC = satexp.human_codons_by_exon()
iu = pd.read_csv("results/independent_units.tsv", sep="\t").set_index("clade")   # paso 32
rows = []
for g in ORDER:
    sub = d[d.group == g].pivot(index="species", columns="exon", values="syn_identity").dropna()
    aa = d[d.group == g].pivot(index="species", columns="exon", values="aa_identity").dropna()
    nt = d[d.group == g].pivot(index="species", columns="exon", values="nt_identity").dropna()
    other = sub[[c for c in REP if c != "E10"]].mean(axis=1)
    fa = "data/mammal_cds.fna" if g == "Mammalia" else "data/tetrapod_cds.fna"
    exp = satexp.expected(set(d[d.group == g].species), fa, REP, HC)
    rows.append([g, len(sub), f"{aa.E10.mean():.3f}", f"{aa.E11.mean():.3f}",
                 f"{nt.E10.mean():.3f}", f"{nt.E11.mean():.3f}",
                 f"{sub.E10.mean():.3f}", f"{other.mean():.3f}",
                 f"{sub.E10.mean() - other.mean():+.3f}",
                 f"{exp['E10']:.3f}",
                 f"{int((sub.E10 > other).sum())}/{len(sub)}",
                 iu.loc[g, "unit"], f"{int(iu.loc[g, 'units_E10_higher'])}/{int(iu.loc[g, 'units'])}",
                 (f"{iu.loc[g, 'sign_test_p']:.3f}" if iu.loc[g, "sign_test_p"] >= 0.01
                  else f"{iu.loc[g, 'sign_test_p']:.1e}".replace("e-0", "e-"))])
t4 = pd.DataFrame(rows, columns=["Clade", "Species", "aa E10", "aa E11", "nt E10", "nt E11",
                                 "syn E10", "syn other exons", "Difference",
                                 "syn E10 expected if free", "Species with E10 higher",
                                 "Unit", "Units with E10 higher", "Sign test p"])
write("table7_outgroup_panel", t4,
      "Table 7. Conservation of exon 10 relative to exons 9, 11, 12 and 13, "
      "measured against human in mammals and in non-mammalian tetrapods. Synonymous identity "
      "counts only codons that encode the same amino acid as human, so it reflects conservation of "
      "the nucleotide sequence once the constraint on the protein is removed. Exon 10 is compared "
      "against the mean of exons 9, 11, 12 and 13 within each species. The expected column gives the "
      "synonymous identity expected in exon 10 if its synonymous sites were free, computed from the "
      "codon usage of each clade, so a clade whose observed value falls below it is measuring "
      "saturated sites rather than conservation. Species of a clade share ancestry, so inference "
      "uses independent units, each contributing the mean difference of its species, which are "
      "orders for mammals and birds and families for the other clades. The sign test asks whether "
      "exon 10 is higher in more units than expected by chance, and the count by species is "
      "descriptive.",
      raw=True)

# ---- Tabla 5. Hyracoidea
loc = pd.read_csv("results/hyrax/locus_summary.tsv", sep="\t")
t5 = pd.DataFrame({
    "Species": [s.replace("_", " ") for s in loc.species],
    "Assembly": loc.assembly,
    "Sequence": [s.split("|")[1] for s in loc.locus_seq],
    "Span (bp)": (loc.end - loc.start),
    "Model length (aa)": [733, 727],
    "Identity to elephant": [0.898, 0.900],
    "Identity to manatee": [0.894, 0.897],
    "Identity to human": [0.816, 0.818],
})
write("tableS1_hyracoidea", t5,
      "Table S1. MAPT loci recovered from two hyrax genome assemblies and the proteins modelled from "
      "them by exonerate. Identities come from the MAFFT alignment of the two models with the long "
      "isoforms of Loxodonta africana, Trichechus manatus latirostris and the human 776 residue "
      "isoform.", md="manuscript/suplementarios/tables")

# ---- Recursos en línea
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
# ESM_4 reúne las accesiones de todo lo que se analizó, mamíferos y no mamíferos. El manuscrito
# promete la accesión de cada proteína y transcrito usados, así que ninguna especie puede faltar.
mam = sel[["species", "superorder", "order", "gene_id", "protein", "transcript", "source",
           "length", "core_coverage", "missing_core_exons", "included", "exclusion_reason"]].copy()
mam.insert(1, "clade", "Mammalia")
tet = pd.read_csv("results/tetrapod_accessions.tsv", sep="\t")
or1 = pd.concat([mam, tet], ignore_index=True)
or1["source"] = or1["source"].fillna("NCBI Datasets orthologs of human MAPT (Gene ID 4137)")
or1 = or1[["species", "clade", "superorder", "order", "family", "gene_id", "protein", "transcript",
           "source", "length", "core_coverage", "E10_coverage", "missing_core_exons", "included",
           "exclusion_reason"]]
from esm import write_esm


def esm(df, n, name):
    write_esm(df, n, name, TSV)


esm(or1, 4, "accessions")


alt = pd.read_csv("results/alt_exon_presence_by_order.tsv", sep="\t", header=[0, 1], index_col=0)
alt.columns = [a for a, _ in alt.columns]
esm(alt.reset_index(), 8, "alt_exons")

esm(d, 13, "tetrapod_per_species")

# ESM_7, soporte por orden del árbol de máxima verosimilitud (lo produce el paso 20)
esm(pd.read_csv("results/phylo_protein/order_support.tsv", sep="\t"), 7, "order_support")

# ESM_2, los transcritos de referencia humanos de los que sale el mapa de exones (paso 3)
esm(pd.read_csv("results/human_reference_transcripts.tsv", sep="\t"), 2, "human_reference_transcripts")
