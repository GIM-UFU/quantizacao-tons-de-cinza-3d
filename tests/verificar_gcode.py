"""Confere um G-code gerado por src/gcode.py, lendo o ARQUIVO (não as estruturas internas).

Verifica: temperaturas, limites da mesa, alturas de Z, filamento extrudido, retrações,
cobertura da seção da peça por linhas de 0,4 mm (preenche? extrapola?) e gera prévias.

Uso: python tests/verificar_gcode.py <arquivo.gcode> <placa.stl> [--leito 220 220]
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from shapely.affinity import translate
from shapely.geometry import LineString
from shapely.ops import unary_union

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from gcode import poligonos, secao  # noqa: E402


def ler(caminho):
    camadas, atual = {}, None
    x = y = z = 0.0
    e = 0.0
    temps, retr, vmax = [], 0, 0
    absoluto = True
    with open(caminho, encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith(";LAYER:"):
                atual = int(ln.split(":")[1])
                camadas[atual] = {"ext": [], "viag": [], "z": None}
                continue
            if ln.startswith(";End of Gcode") or (atual is not None and ln.startswith("G91")):
                break  # a partir daqui é o bloco final (movimentos relativos)
            if re.match(r"^M(104|109|140|190) ", ln):
                temps.append(ln.split(";")[0].strip())
            if ln.startswith("G91"):
                absoluto = False
            if atual is None or not re.match(r"^G[01] ", ln):
                continue
            cod = ln.split(";")[0]
            vals = {k: float(v) for k, v in re.findall(r"([XYZEF])(-?\d+\.?\d*)", cod)}
            nx, ny, nz = vals.get("X", x), vals.get("Y", y), vals.get("Z", z)
            ne = vals.get("E", e)
            if "F" in vals:
                vmax = max(vmax, vals["F"])
            if "X" in vals or "Y" in vals:
                seg = ((x, y), (nx, ny))
                if ne > e + 1e-9:
                    camadas[atual]["ext"].append(seg)
                else:
                    camadas[atual]["viag"].append(seg)
            elif ne < e - 1e-9:
                retr += 1
            if "Z" in vals:
                camadas[atual]["z"] = nz
            x, y, z, e = nx, ny, nz, ne
    return camadas, temps, retr, vmax, e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gcode")
    ap.add_argument("stl")
    ap.add_argument("--leito", type=float, nargs=2, default=[220, 220])
    a = ap.parse_args()

    cam, temps, retr, vmax, e_final = ler(a.gcode)
    ext = [s for c in cam.values() for s in c["ext"]]
    xs = [p[0] for s in ext for p in s]
    ys = [p[1] for s in ext for p in s]
    zs = sorted({c["z"] for c in cam.values() if c["z"] is not None})
    comp = sum(math.hypot(s[1][0] - s[0][0], s[1][1] - s[0][1]) for s in ext)
    fil_mm = comp * 0.4 * 0.2 / (math.pi * 0.875 ** 2)
    ok_mesa = min(xs) >= 0 and max(xs) <= a.leito[0] and min(ys) >= 0 and max(ys) <= a.leito[1]

    m = trimesh.load(a.stl)
    P = secao(m, m.extents[2] / 2)
    with open(a.gcode, encoding="utf-8") as fh:  # mesmo giro que o gerador aplicou
        cab = fh.read(300)
    g = re.search(r";ROTATION_Z_DEG:(-?[\d.]+)", cab)
    if g and float(g.group(1)):
        from shapely.affinity import rotate
        P = rotate(P, float(g.group(1)), origin=unary_union(poligonos(P)).centroid)
        print(f"giro aplicado pelo gerador: {float(g.group(1)):g}°")
    mnx, mny, mxx, mxy = P.bounds
    P = translate(P, a.leito[0] / 2 - (mnx + mxx) / 2, a.leito[1] / 2 - (mny + mxy) / 2)
    ilhas = unary_union(poligonos(P))

    print(f"arquivo: {a.gcode} ({os.path.getsize(a.gcode) / 1e6:.2f} MB)")
    print(f"temperaturas: {temps}")
    print(f"camadas: {len(cam)} · Z de {zs[0]:.2f} a {zs[-1]:.2f} mm · passo(s) de Z: {sorted({round(b - a_, 3) for a_, b in zip(zs, zs[1:])})}")
    print(f"área de extrusão X {min(xs):.1f}..{max(xs):.1f} · Y {min(ys):.1f}..{max(ys):.1f} (mesa {a.leito[0]:g}×{a.leito[1]:g}) → dentro da mesa: {ok_mesa}")
    print(f"percurso extrudado: {comp / 1000:.1f} m · filamento (conta independente): {fil_mm / 1000:.2f} m · "
          f"massa PLA 1,24 g/cm³: {fil_mm * math.pi * 0.875 ** 2 / 1000 * 1.24:.1f} g")
    print(f"retrações: {retr} · velocidade máx. de avanço: F{vmax:.0f} (= {vmax / 60:.0f} mm/s)")
    print(f"E final (deve voltar a ficar retraído, perto de 0 após o último avanço): {e_final:.2f}")

    # cobertura numa camada intermediária (sem brim): linhas de 0,4 mm contra a seção
    k = len(cam) // 2
    cob = unary_union([LineString(s).buffer(0.2, cap_style=2) for s in cam[k]["ext"]])
    dentro = cob.intersection(ilhas).area
    print(f"\ncamada {k} (z = {cam[k]['z']:.1f} mm): seção {ilhas.area:.1f} mm² · coberta pelas linhas {dentro / ilhas.area:.1%} · "
          f"material fora da seção {cob.difference(ilhas).area:.2f} mm² ({cob.difference(ilhas).area / ilhas.area:.2%})")
    cb0 = unary_union([LineString(s).buffer(0.2, cap_style=2) for s in cam[0]["ext"]])
    print(f"camada 0: cobertura da seção {cb0.intersection(ilhas).area / ilhas.area:.1%} · "
          f"brim (área fora da seção) {cb0.difference(ilhas).area / 100:.1f} cm²")

    fig, ax = plt.subplots(1, 3, figsize=(18, 7))
    bx0, by0, bx1, by1 = ilhas.buffer(8).bounds
    for a_, (kk, tit) in zip(ax[:2], [(0, "camada 0 (brim + 1ª camada)"), (k, f"camada {k}")]):
        for s in cam[kk]["ext"]:
            a_.plot(*zip(*s), color="#1f77b4", lw=0.6)
        for s in cam[kk]["viag"]:
            a_.plot(*zip(*s), color="#d62728", lw=0.3, alpha=0.5)
        a_.set_title(tit + " · azul = extrusão, vermelho = deslocamento", fontsize=9)
    for p in poligonos(ilhas):
        ax[2].fill(*p.exterior.xy, color="#999", alpha=0.7)
    ax[2].set_title("seção da peça (cinza) na mesa", fontsize=9)
    for a_ in ax:
        a_.set_xlim(bx0, bx1); a_.set_ylim(by0, by1); a_.set_aspect("equal")
    fig.tight_layout()
    saida = os.path.splitext(a.gcode)[0] + "_previa.png"
    fig.savefig(saida, dpi=80)
    print("prévia:", saida)


if __name__ == "__main__":
    main()
