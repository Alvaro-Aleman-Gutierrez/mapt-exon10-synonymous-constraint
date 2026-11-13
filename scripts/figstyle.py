"""Estilo común de las figuras, según las especificaciones de arte del envío.

Se usa tipografía Helvetica o Arial, cuerpo de 2 a 3 mm (8 a 12 pt) en el tamaño final y
variación mínima dentro de la ilustración, líneas de al menos 0.3 pt, color RGB, y ningún título ni
pie dentro de la figura. Las figuras compuestas se rotulan (a), (b). El ancho de columna es 84 mm y
se admiten 1.5 o 2 columnas cuando el detalle lo exige.

En EPS no hay canal alfa, así que ninguna figura usa transparencia. Los puntos superpuestos se
separan con un anillo blanco.
"""
import glob
import os
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from PIL import Image

MM = 1 / 25.4
W1, W15, W2 = 84 * MM, 129 * MM, 174 * MM   # 1, 1.5 y 2 columnas

# paleta validada con el verificador de paletas, separación CVD mínima ΔE 10.9
FOCUS = "#B4700E"   # exón 10, la entidad de interés
REF = "#2A6F9E"     # las demás repeticiones
MARK = "#8C2D5A"    # posiciones de FTDP-17 que actúan por splicing
INK = "#1a1a1a"
MUTED = "#595959"
GRID = "#dcdcdc"
PAPER = "#ffffff"

def _arial():
    """Registra Arial si está en el sistema."""
    pats = ["/mnt/c/Windows/Fonts/arial*.ttf", "/usr/share/fonts/**/LiberationSans*.ttf",
            "/usr/share/fonts/**/NimbusSans*.otf"]
    found = [f for p in pats for f in glob.glob(p, recursive=True)]
    for f in found:
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass
    names = {f.name for f in fm.fontManager.ttflist}
    for want in ("Arial", "Liberation Sans", "Nimbus Sans", "Helvetica"):
        if want in names:
            return want
    return "DejaVu Sans"

FONT = _arial()

# cuerpo de 8 pt en todo, 9 pt solo en los rótulos de panel, dentro del rango pedido
mpl.rcParams.update({
    "font.family": FONT, "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "grid.linewidth": 0.5, "lines.linewidth": 0.8,
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "pdf.fonttype": 42, "ps.fonttype": 42,  # TrueType incrustada en PDF y EPS
    "svg.fonttype": "none",
})

def panel(ax, letter, dx=-0.085, dy=1.045):
    """Rótulo de panel, (a), (b) en las figuras compuestas."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)

def despine(ax, keep=("left", "bottom")):
    for s, sp in ax.spines.items():
        sp.set_visible(s in keep)

def save(fig, num, outdir="manuscript/figuras", tiff_dpi=600):
    """Guarda EPS vectorial para el envío, TIFF para el respaldo y PNG para revisar.

    El nombre empieza por Fig y el número.
    """
    os.makedirs(outdir, exist_ok=True)
    base = f"{outdir}/Fig{num}"
    fig.savefig(base + ".eps", format="eps")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=tiff_dpi)
    im = Image.open(base + ".png").convert("RGB")
    im.save(base + ".tif", format="TIFF", compression="tiff_lzw",
            dpi=(tiff_dpi, tiff_dpi))
    sizes = {e: os.path.getsize(f"{base}.{e}") // 1024 for e in ("eps", "tif", "png")}
    print(f"Fig{num}  {im.size[0]}x{im.size[1]} px a {tiff_dpi} dpi  "
          + "  ".join(f"{k} {v} kB" for k, v in sizes.items()))


try:
    from journal_meta import JOURNAL
except ImportError:
    JOURNAL = None

TITLE = ("Synonymous constraint on MAPT exon 10 exceeds what its encoded protein requires and is "
         "shared across amniotes")


def save_esm(fig, num, description, outdir="manuscript/figuras"):
    """Guarda una figura suplementaria como ESM_n.pdf con la portada identificatoria.

    Cada archivo suplementario lleva, dentro de él, el título del artículo, los autores y la
    afiliación y el correo del autor de correspondencia, y una línea más si existe journal_meta.py.
    En un PDF eso se resuelve con una primera página, que no toca las dimensiones de la figura.
    """
    import io
    import textwrap
    from matplotlib.backends.backend_pdf import PdfPages
    from pypdf import PdfReader, PdfWriter

    os.makedirs(outdir, exist_ok=True)
    body = io.BytesIO()
    fig.savefig(body, format="pdf")

    cover = io.BytesIO()
    with PdfPages(cover) as pdf:
        c = plt.figure(figsize=(210 * MM, 297 * MM))
        lines = [("Online Resource %d" % num, 15, "bold"), ("", 9, "normal"),
                 ("Article title", 8.5, "bold"), (textwrap.fill(TITLE, 62), 10.5, "normal"),
                 ("", 9, "normal"),
                 *([("Journal", 8.5, "bold"), (JOURNAL, 10.5, "normal"), ("", 9, "normal")]
                   if JOURNAL else []),
                 ("Author", 8.5, "bold"), ("Alvaro Alemán-Gutiérrez", 10.5, "normal"),
                 ("", 9, "normal"),
                 ("Corresponding author", 8.5, "bold"),
                 (textwrap.fill("Alvaro Alemán-Gutiérrez, Laboratorio de Biología Molecular, "
                                "Escuela de Ciencias Naturales y Exactas, Campus Tecnológico Local "
                                "San Carlos, Instituto Tecnológico de Costa Rica, Alajuela, Costa Rica",
                                62), 10.5, "normal"),
                 ("aaleman@itcr.ac.cr", 10.5, "normal"),
                 ("", 9, "normal"),
                 ("Contents", 8.5, "bold"), (textwrap.fill(description, 62), 10.5, "normal")]
        y = 0.90
        for text, size, weight in lines:
            n = max(1, text.count("\n") + 1)
            c.text(0.12, y, text, fontsize=size, fontweight=weight, va="top", ha="left")
            y -= 0.0195 * n + (0.012 if text else 0.004)
        pdf.savefig(c)
        plt.close(c)

    out = PdfWriter()
    for src in (cover, body):
        src.seek(0)
        for page in PdfReader(src).pages:
            out.add_page(page)
    path = f"{outdir}/ESM_{num}.pdf"
    with open(path, "wb") as f:
        out.write(f)
    print(f"ESM_{num}.pdf  {os.path.getsize(path) // 1024} kB, con portada")
