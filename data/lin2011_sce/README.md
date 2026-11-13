# Pistas de Lin et al. (2011) para MAPT

Lin MF, Kheradpour P, Washietl S, Parker BJ, Pedersen JS, Kellis M (2011) Locating protein-coding
sequences under selection for additional, overlapping functions in 29 mammalian genomes.
Genome Res 21:1916-1928.

## De dónde salen

El archivo `SynonymousConstraintElements.tar.gz` (30,6 MB, fechado el 14 de abril de 2011), que
distribuyó el Broad Institute en
`http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/`. Se recuperó de la copia del
Internet Archive del 14 de septiembre de 2015:

`https://web.archive.org/web/20150914091453id_/http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/SynonymousConstraintElements.tar.gz`

Se usó la copia archivada porque el sitio del MIT que alojó los archivos originalmente
(`compbio.mit.edu/SCE/`) ya no los sirve. El 21 de octubre de 2026 devolvía 404.

El archivo contiene `SCE9.hg18.bed`, `SCE15.hg18.bed`, `SCE30.hg18.bed`, que son los elementos a
resolución de 9, 15 y 30 codones, y `lambda_s_ORF.hg18.bedGraph`, la tasa sinónima local relativa
al promedio del marco de lectura. Todo en coordenadas hg18.

## Qué hay en esta carpeta

Solo la región de MAPT, extraída del archivo completo:

- `SCE9/15/30.MAPT.hg18.bed`: las filas del cromosoma 17 cuyo nombre de gen es MAPT.
- `lambda_s.MAPT.hg18.bedGraph`: las ventanas del cromosoma 17 entre 41 395 557 y 41 467 498.
- `hg18ToHg38.over.chain.gz`: cadena oficial de UCSC para pasar de hg18 a GRCh38, de
  `https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz`.

El cruce con los exones humanos lo hace `scripts/33_lin2011_overlap.py`.
