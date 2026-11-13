#!/usr/bin/env bash
# Paso 7: filogenia ML de proteína 2N4R. ModelFinder (BIC) + UFBoot 1000 (-bnni) + SH-aLRT 1000.
set -euo pipefail
cd "$(dirname "$0")/.."
iqtree2 -s results/mapt_2N4R_final.aln.faa -st AA -m MFP -B 1000 --bnni --alrt 1000 \
        -T AUTO --threads-max 10 --seed 12345 --prefix results/phylo_protein/mapt_2N4R -redo

# Paso 7c: árbol ML restringido por la taxonomía de NCBI (paso 7b), topología para HyPhy (paso 9)
iqtree2 -s results/mapt_2N4R_final.aln.faa -st AA -m HIVb+F+R4 -g results/phylo_protein/taxonomy_constraint.nwk \
        -T AUTO --threads-max 10 --seed 12345 --prefix results/phylo_protein/mapt_2N4R_constrained -redo
