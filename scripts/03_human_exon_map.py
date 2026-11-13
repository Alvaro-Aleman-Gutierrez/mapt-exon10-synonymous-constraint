"""Paso 3: mapa de exones codificantes de MAPT humano (unión de todos los transcritos NM_)."""
import json
import pandas as pd
from Bio import SeqIO

reports = json.load(open("data/mammal_product_reports.json"))
human = next(r["product"] for r in reports if str(r["product"]["tax_id"]) == "9606")
prot = {r.id: str(r.seq) for r in SeqIO.parse("data/mammal_protein.faa", "fasta")}
cds = {r.id.split(":")[0]: str(r.seq) for r in SeqIO.parse("data/mammal_cds.fna", "fasta")}

transcripts = []   # accesiones, para listarlas en el suplementario
segments = {}   # (inicio, fin) genómico codificante -> secuencia nt
for t in human["transcripts"]:
    if not t["accession_version"].startswith("NM_"):
        continue
    exons = sorted((int(e["begin"]), int(e["end"])) for e in t["genomic_locations"][0]["exons"])
    tx_len = sum(e - b + 1 for b, e in exons)
    c0 = int(t["cds"]["range"][0]["begin"]); c1 = int(t["cds"]["range"][0]["end"])
    pos = 0
    for b, e in exons:
        s, f = pos + 1, pos + e - b + 1          # coordenadas en el transcrito
        lo, hi = max(s, c0), min(f, c1)
        if lo <= hi:
            segments[(b + lo - s, b + hi - s)] = None
        pos = f
    print(t["accession_version"], t["protein"]["accession_version"], t["protein"]["length"], len(exons))
    transcripts.append({"transcript": t["accession_version"], "protein": t["protein"]["accession_version"],
                        "protein_length": t["protein"]["length"], "exons": len(exons),
                        "assembly_chromosome": "NC_000017.11"})

for k in sorted(segments):
    print(k, k[1] - k[0] + 1, "nt")

pd.DataFrame(transcripts).sort_values("transcript").to_csv(
    "results/human_reference_transcripts.tsv", sep="\t", index=False)
