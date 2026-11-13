"""Copias en tándem del exón 10 en tetrápodos no mamíferos.

Las aves llevan una copia duplicada del exón 10 que se empalma de forma alternativa y produce
isoformas de cinco repeticiones (Yoshida y Goedert 2002), y algunos escamosos llevan más. El
alineamiento global con la 2N4R humana reparte entonces el exón 10 humano entre las copias y mide
una quimera. Este módulo encuentra las copias de forma explícita para que el exón 10 de cada especie
sea una copia definida, la 5', y la 3' quede disponible como exón aparte.

Una copia es una ventana de 31 residuos con al menos 0,55 de identidad con el exón 10 humano y más
parecida a él que al exón 11, buscada sin solapamiento de 5' a 3'.
"""

MIN_ID = 0.55


def ident(a, b):
    return sum(x == y for x, y in zip(a, b)) / len(b)


def find(aa, h10, h11):
    """Posiciones (0-based) de inicio de cada copia del exón 10 en la proteína `aa`."""
    hits, i, L = [], 0, len(h10)
    while i <= len(aa) - L:
        w = aa[i:i + L]
        if ident(w, h10) >= MIN_ID and ident(w[:len(h11)], h11) < ident(w, h10):
            hits.append(i)
            i += L
        else:
            i += 1
    return hits


def remap(m, aa, hp, bounds):
    """Corrige el mapa humano -> especie `m` cuando hay copias en tándem del exón 10.

    El exón 10 humano pasa a la copia 5' y el exón 11 humano a los residuos que siguen a la última
    copia, siempre que esa ventana se parezca al exón 11 humano. Devuelve las posiciones de las
    copias.
    """
    s, e = bounds["E10"]
    s11, e11 = bounds["E11"]
    hits = find(aa, hp[s:e], hp[s11:e11])
    if len(hits) >= 2:
        for j in range(e - s):
            m[s + j] = hits[0] + j
        after = hits[-1] + (e - s)
        if ident(aa[after:after + (e11 - s11)], hp[s11:e11]) >= MIN_ID:
            for j in range(e11 - s11):
                m[s11 + j] = after + j
    return hits
