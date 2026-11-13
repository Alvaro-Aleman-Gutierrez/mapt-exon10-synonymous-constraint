# Contexto genómico del exón 10: datos de entrada

`SCE9/15/30.hg18.bed` y `lambda_s_ORF.hg18.bedGraph` salen del archivo
`SynonymousConstraintElements.tar.gz` de Lin et al. (2011), distribuido por el Broad Institute y
recuperado de la copia del Internet Archive del 14 de septiembre de 2015:

`https://web.archive.org/web/20150914091453id_/http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/SynonymousConstraintElements.tar.gz`

Son los archivos del genoma completo, no solo de MAPT. Todo en coordenadas hg18.

`refGene.hg18.txt.gz` es la tabla de RefSeq del UCSC Genome Browser para hg18, bajada el
22 de octubre de 2026 de
`https://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz`.

Estos archivos no van en el depósito de Zenodo por su tamaño (unos 70 MB), y son públicos. El
paso 34 los lee desde esta carpeta. Se descartaron el archivo comprimido original y su
`HEADER.html`, que no hacen falta.
