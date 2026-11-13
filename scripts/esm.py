"""Escritura de los suplementarios tabulares como ESM_n, con el bloque identificatorio.

Los suplementarios se numeran de forma consecutiva y cada uno lleva dentro el título del artículo,
los autores y la afiliación y el correo del autor de correspondencia, y una línea más si
existe journal_meta.py. En un TSV eso va como líneas de comentario al inicio, que pandas salta con
`comment="#"`. Lo usan los pasos 16, 21 y 22 para que ningún suplementario dependa de código suelto.
"""
import os

try:
    from journal_meta import JOURNAL
except ImportError:
    JOURNAL = None

TITLE = ("Synonymous constraint on MAPT exon 10 exceeds what its encoded protein requires and is "
         "shared across amniotes")
HDR = (f"# {TITLE}\n" + (f"# {JOURNAL}\n" if JOURNAL else "") +
       "# Alvaro Alemán-Gutiérrez, Instituto Tecnológico de Costa Rica, "
       "aaleman@itcr.ac.cr\n")


def write_esm(df, n, name, outdir="results/tables"):
    os.makedirs(outdir, exist_ok=True)
    path = f"{outdir}/ESM_{n}_{name}.tsv"
    with open(path, "w") as fh:
        fh.write(HDR + f"# Online Resource {n}\n")
        df.to_csv(fh, sep="\t", index=False)
    print(f"  ESM_{n}_{name}.tsv  {len(df)} filas")
    return path
