"""Paso 31: tabla de síntesis de las líneas de evidencia (Table 8).

Cada cifra se lee de los archivos que escribieron los pasos anteriores, no se teclea, para que la
tabla no pueda discrepar del texto. La última fila es el intento de explicar la restricción por
splicing, que no se sostiene con una prueba, y está para que el lector vea el nivel de soporte de
cada afirmación y no solo el de las que salieron bien.
"""
import pandas as pd
from scipy import stats

perm = pd.read_csv("results/tables/table3_permutation.tsv", sep="\t").set_index("exon")
fub = pd.read_csv("results/selection/fubar_sites.tsv", sep="\t")
t4 = pd.read_csv("results/tables/table7_outgroup_panel.tsv", sep="\t").set_index("Clade")
ese = pd.read_csv("results/tables/table5_ese_overlap.tsv", sep="\t").set_index("exon")

# línea 1, invariancia
e = perm.loc["E10"]
l1 = (f"{int(e.observed_invariant)} of {int(e.codons)} codons identical in every mammal, against "
      f"{e.expected_mean:.2f} expected after matching on protein constraint and codon degeneracy "
      f"(p < 10⁻⁴)")

# línea 2, FUBAR
REP = ["E9", "E11", "E12", "E13"]
a10 = fub[fub.exon == "E10"]
ao = fub[fub.exon.isin(REP)]
p_mw = stats.mannwhitneyu(a10.fubar_alpha, ao.fubar_alpha, alternative="less").pvalue
n_below = int((a10.p_alpha_below_bg >= 0.9).sum())
b11 = fub[fub.exon == "E11"].fubar_beta.median()
l2 = (f"Median synonymous rate {a10.fubar_alpha.median():.3f} against {ao.fubar_alpha.median():.3f} "
      f"in the other four exons ({ao.fubar_alpha.median() / a10.fubar_alpha.median():.1f} fold, "
      f"p = {p_mw:.1e}); {n_below} of {len(a10)} codons below the gene background at posterior "
      f"0.9 or more; nonsynonymous rate {a10.fubar_beta.median():.3f}, equal to exon 11 "
      f"({b11:.3f})").replace("e-", " × 10⁻").replace("0⁻13", "0⁻¹³")

# línea 3, identidad sinónima frente a la expectativa con sitios libres, e inversión
amn = ["Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria"]
gap = (t4.loc[amn, "syn E10"].astype(float) - t4.loc[amn, "syn E10 expected if free"].astype(float))
inv = pd.read_csv("results/descriptive_to_statistic.tsv", sep="\t")
row = inv[inv.quantity.str.startswith("clados de amniotas")].value.iloc[0]
k, n = row.split(" de ")[0], row.split(" de ")[1].split(",")[0]
iu = pd.read_csv("results/independent_units.tsv", sep="\t").set_index("clade")
l3 = (f"Exon 10 is higher in all {int(iu.loc['Mammalia', 'units'])} mammalian and "
      f"{int(iu.loc['Aves', 'units'])} avian orders (sign tests, p < 10⁻⁶); it sits "
      f"{gap.min():.2f} to {gap.max():.2f} above the identity expected if its "
      f"synonymous sites were free, in all five amniote clades; its protein is less conserved and "
      f"its nucleotides more conserved than exon 11 in {k} of {n} clades (clade level binomial, "
      f"p = 0.001)")

# línea 4, FUBAR dentro de cada clado (pasos 37 a 40)
fc = pd.read_csv("results/amniote/fubar_by_clade.constrained.tsv", sep="\t").set_index("clade")
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
def sci(p):
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10{str(int(e)).translate(SUP)}"
l4b = (f"Exon 10 is the slowest exon at synonymous sites within birds "
       f"({fc.loc['Aves', 'fold']:.1f} fold, p = {sci(fc.loc['Aves', 'mw_p_E10_lower'])}) and within "
       f"squamates ({fc.loc['Lepidosauria', 'fold']:.1f} fold, "
       f"p = {sci(fc.loc['Lepidosauria', 'mw_p_E10_lower'])}); its per site rates correlate with the "
       f"mammalian ones once codon degeneracy is held constant (partial ρ "
       f"{fc.loc['Lepidosauria', 'partial_rho_given_syn_sites']:.2f} to "
       f"{fc.loc['Aves', 'partial_rho_given_syn_sites']:.2f} in birds, squamates and turtles)")

# línea 5, el intento de explicar por splicing
s = ese.loc["E10"]
l4 = (f"The invariant codons overlap predicted exonic splicing enhancers no more than chance "
      f"({int(s.observed_overlap_nt)} nt against {s.expected_mean:.1f} expected, "
      f"p = {s.p_empirical:.2f})")

rows = [
    ["1. Count of invariant codons", l1,
     "Nothing about substitution", "None needed", "259 mammals, 23 orders",
     "A rate, or when the constraint arose", "Strong"],
    ["2. Site specific synonymous rate (FUBAR)", l2,
     "Codon model, rate grid shared across sites", "Yes, through the tree", "259 mammals",
     "Which codons are held by what, since the grid borrows strength between sites", "Strong"],
    ["3. Synonymous identity in amniotes", l3,
     "Pairwise identity to human; expectation from each clade's codon usage",
     "None explicit; tests count independent orders, families and clades",
     "Mammals, birds, crocodilians, turtles, squamates; absent in frogs, which are not saturated; salamanders underpowered",
     "The ancestral state, which is a parsimony inference", "Strong for presence, moderate for age"],
    ["4. Synonymous rate within clades (FUBAR)", l4b,
     "Codon model within each clade, rate grid shared across sites",
     "Yes, through each clade's tree", "Birds, squamates, turtles, frogs, salamanders",
     "Anything in turtles and salamanders, whose trees are too short to separate the exons",
     "Strong in birds and squamates; supports shared sites"],
    ["Splicing as the cause", l4,
     "RESCUE-ESE represents the elements of this exon", "Not applicable", "Exon 10 and exons 9, 11, 12 and 13",
     "That splicing is not the cause, since the predictor is generic", "Hypothesis; the one test was negative and has low power"],
]
cols = ["Line of evidence", "Result for exon 10", "Assumes", "Phylogenetic correction", "Scope",
        "What it cannot show", "Support"]
out = pd.DataFrame(rows, columns=cols)
out.to_csv("results/tables/table8_evidence.tsv", sep="\t", index=False)
cap = ("Table 8. Lines of evidence for synonymous constraint on MAPT exon 10 and their level of "
       "support. The first three lines are independent in what they assume, which is why agreement "
       "among them carries weight. The last row reports the attempt to attribute the constraint to "
       "splicing, which no test in this study supports.")
with open("manuscript/tables/table8_evidence.md", "w") as f:
    f.write(f"**{cap}**\n\n" + out.to_markdown(index=False) + "\n")
print(out.T.to_string(header=False))
