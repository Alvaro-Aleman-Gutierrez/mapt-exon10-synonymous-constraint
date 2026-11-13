# Scripts and derived data for the reanalysis of MAPT evolution in mammals

Alvaro Alemán-Gutiérrez

This record holds the code and the derived data behind a study of the evolution of *MAPT* (tau)
in mammals and other tetrapods. Every analysis starts from public accessions, so the scripts
regenerate each table and figure from scratch.

## What the study did

One homologous 2N4R-equivalent sequence was selected per species for 259 mammals from 23 orders,
built where necessary by in silico excision of the exons that the human 2N4R isoform lacks, so that
exon-level comparisons do not depend on which isoform a database happened to annotate. Codon models
in HyPhy then separated synonymous from nonsynonymous divergence. A panel of 236 non-mammalian
tetrapods places the mammalian result in a longer time frame, and *MAPT* was reconstructed from two
hyrax genomes because no annotated ortholog exists for Hyracoidea.

## Layout

| Folder | Contents |
| --- | --- |
| `scripts/` | The numbered Python and shell scripts of the analyses, 01 to 40, in the order they run, with the shared modules they import. Steps 18, 19, 23, 26 and 36 only assemble the manuscript and the deposit and are not included |
| `data/` | Inputs retrieved from NCBI, the human exon map, taxonomy tables, and the subset of Lin et al. (2011) synonymous constraint elements that overlaps *MAPT* |
| `alignments/` | The selected 2N4R sequences, the MAFFT protein alignment and the derived codon alignment |
| `trees/` | IQ-TREE output for the unconstrained gene tree and for the taxonomy-constrained tree, plus the constraint tree itself |
| `selection/` | HyPhy JSON and logs for FEL, MEME, SLAC per order and the three RELAX runs, with the per-order subalignments and pruned trees |
| `amniote/` | Exon 9 to 13 codon alignments, trees and FUBAR output for each non-mammalian clade (steps 37 to 40), with the NCBI Taxonomy lineages used for the constraint trees |
| `hyrax/` | tblastn hits, the extracted loci, exonerate models and the afrotherian alignment |
| `tables/` | The tables of the study and its supplementary tables (Online Resources), as tab-separated files |
| `figures/` | Figures in EPS, TIFF and PNG, and the two supplementary figures in PDF |
| `supplementary/` | The three supplementary documents (Online Resources 1, 3 and 5) |

## Running it

The scripts expect to run from the directory above `scripts/`, in the order their numbers give.
Steps 01, 02, 10, 14 and 38 download from NCBI and need network access. Step 10 downloads two genome
assemblies of about 3.6 GB each and builds BLAST databases from them, which this record does not
include, because they are public and the manuscript cites them by accession. Those assemblies are
GCA_004026925.3 for *Procavia capensis* and GCA_028571685.1 for *Heterohyrax brucei*.

Software versions are the ones the manuscript reports, which are Biopython 1.83, MAFFT v7.490,
IQ-TREE 2.0.7, HyPhy 2.5.103, BLAST+ 2.12.0 and exonerate 2.4.0. HyPhy and exonerate ran in their
own conda environments, and the scripts call them as `hyphy` and `exonerate`, so both must be on the
`PATH` when the steps that use them run. Steps 01, 07b and 14 query NCBI through Biopython, which
asks for a contact address. Replace `your.email@example.org` in those scripts with your own.

## Notes on what is not here

Genome assemblies and BLAST databases are left out for the reason above. IQ-TREE checkpoint files
and the bootstrap replicate file are also left out, since neither is needed to reproduce a result.
Primary sequence data are not duplicated beyond what the scripts downloaded, and Online Resource 4
in `tables/` lists the accession of every protein and transcript used, species by species, with the
reason each excluded species was dropped.

## License

Code is released under the MIT license. Data files are released under CC BY 4.0.
