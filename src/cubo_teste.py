"""Corpo de prova para testar a infiltração de cera de abelha no PLA antes da peça única.

Uso:
    python src/cubo_teste.py                 # saida/10_cubo_teste/
    python src/cubo_teste.py --camada 0.2 --saida saida/10_cubo_teste

O que o bloco testa (ver src/peca_unica.py):
  1. Térmico: a cera de abelha funde a ≈ 62–64 °C e o PLA amolece perto de 55–60 °C (Tg).
     As paredes finas (1,2 mm = traço mínimo da peça única) e as dimensões externas conhecidas
     mostram se o PLA deformou: meça com paquímetro antes e depois da cera.
  2. Preenchimento: fendas de 0,8 a 5 mm e furos de 1 a 3 mm (cotas múltiplas do pixel de 0,2 mm), todos abertos para cima (como
     os vazios da peça única), mostram a menor abertura que a cera preenche sem bolha.
  3. Contração: a cera encolhe ao solidificar; o rebaixo largo (≥ 5 mm) mostra o rechupe.
  4. Degrau (bônus): 5 degraus de PLA (2 a 10 mm) para radiografar e medir o μ do PLA impresso,
     que é a validação pendente do fator de material.

Imprime deitado, sem suporte (todos os vazios abrem para cima), 100 % de preenchimento.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh

sys.path.insert(0, os.path.dirname(__file__))
from mama3d import solido_entre  # noqa: E402

# geometria (mm)
COMP, LARG, ALT = 60.0, 40.0, 15.0     # bloco: X × Y × Z
FUNDO = 2.0                            # piso sob os vazios
PAREDE = 1.2                           # paredes entre fendas = traço mínimo da peça única
FENDAS = [0.8, 1.0, 1.6, 2.0, 3.0, 5.0]  # larguras das fendas (mm), comprimento 14 mm; múltiplos de 0,2 (pixel)
FUROS = [1.0, 1.6, 2.0, 3.0]           # diâmetros dos furos (mm)
DEGRAUS = [2.0, 4.0, 6.0, 8.0, 10.0]   # alturas de PLA do degrau (mm), 8 mm de largura cada


def mapa_altura(px: float) -> tuple[np.ndarray, list[dict]]:
    nx, ny = int(round(COMP / px)), int(round(LARG / px))
    x = (np.arange(nx) + 0.5) * px
    y = (np.arange(ny) + 0.5) * px
    X, Y = np.meshgrid(x, y)  # linha = y, coluna = x
    h = np.full((ny, nx), ALT)
    feitos = []

    # faixa 1 (y 3–17): fendas paralelas a Y, separadas por paredes de 1,2 mm
    x0, y0, y1 = 4.0, 3.0, 17.0
    for w in FENDAS:
        h[(X >= x0) & (X < x0 + w) & (Y >= y0) & (Y < y1)] = FUNDO
        feitos.append({"tipo": "fenda", "largura_mm": w, "x_mm": [x0, x0 + w], "y_mm": [y0, y1],
                       "profundidade_mm": ALT - FUNDO})
        x0 += w + PAREDE
    # rebaixo largo (contração da cera), à direita das fendas
    xr0 = x0 + 2.0
    h[(X >= xr0) & (X < xr0 + 12.0) & (Y >= y0) & (Y < y1)] = FUNDO
    feitos.append({"tipo": "rebaixo", "largura_mm": 12.0, "x_mm": [xr0, xr0 + 12.0], "y_mm": [y0, y1],
                   "profundidade_mm": ALT - FUNDO})

    # faixa 2 (y 20–26): furos redondos
    xc = 6.0
    for dmt in FUROS:
        h[(X - xc) ** 2 + (Y - 23.0) ** 2 <= (dmt / 2) ** 2] = FUNDO
        feitos.append({"tipo": "furo", "diametro_mm": dmt, "centro_mm": [xc, 23.0], "profundidade_mm": ALT - FUNDO})
        xc += dmt + 4.0

    # faixa 3 (y 29–37): degrau de PLA (o resto da faixa é rebaixado até o fundo)
    xd = 4.0
    for hd in DEGRAUS:
        h[(X >= xd) & (X < xd + 8.0) & (Y >= 29.0) & (Y < 37.0)] = hd
        feitos.append({"tipo": "degrau", "altura_PLA_mm": hd, "x_mm": [xd, xd + 8.0], "y_mm": [29.0, 37.0]})
        xd += 8.0 + 2.0
    return h, feitos


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camada", type=float, default=0.2, help="altura de camada em mm (padrão 0.2)")
    ap.add_argument("--pixel", type=float, default=0.2, help="resolução lateral da malha em mm (padrão 0.2)")
    ap.add_argument("--saida", default=os.path.join("saida", "10_cubo_teste"))
    a = ap.parse_args(argv)
    os.makedirs(a.saida, exist_ok=True)

    h, feitos = mapa_altura(a.pixel)
    h = np.round(h / a.camada) * a.camada
    mask = np.ones_like(h, dtype=bool)
    # solido_entre usa (linha → −y); inverter as linhas faz o y da malha coincidir com o do desenho
    hm = h[::-1]
    m = solido_entre(np.zeros_like(hm), hm - a.camada, mask, a.pixel, dz_mm=a.camada, suavizar=(0.5, 0.5, 0))
    m.apply_translation(-m.bounds[0])
    m.export(os.path.join(a.saida, "cubo_teste.stl"))
    yup = m.copy()  # documento Y-up do Fusion: altura em Y, base em y = 0, centrado em X e Z
    yup.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    c = yup.bounds.mean(axis=0)
    yup.apply_translation([-c[0], -yup.bounds[0][1], -c[2]])
    yup.export(os.path.join(a.saida, "cubo_teste_yup.stl"))

    fig, ax = plt.subplots(figsize=(9, 6.2))
    im = ax.imshow(h, cmap="viridis", origin="lower", extent=[0, COMP, 0, LARG])
    fig.colorbar(im, ax=ax, label="altura de PLA [mm]")
    for f in feitos:
        if f["tipo"] in ("fenda", "rebaixo"):
            ax.text(np.mean(f["x_mm"]), 18.0, f"{f['largura_mm']:g}", ha="center", fontsize=7, color="k")
        elif f["tipo"] == "furo":
            ax.text(f["centro_mm"][0], 25.8, f"⌀{f['diametro_mm']:g}", ha="center", fontsize=7, color="k")
        else:
            ax.text(np.mean(f["x_mm"]), 37.6, f"{f['altura_PLA_mm']:g}", ha="center", fontsize=7, color="k")
    ax.set_title(f"corpo de prova cera × PLA · {COMP:g} × {LARG:g} × {ALT:g} mm (vista de cima, cotas em mm)", fontsize=9)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    fig.tight_layout()
    fig.savefig(os.path.join(a.saida, "desenho.png"), dpi=130)
    plt.close(fig)

    vol = float(m.volume) / 1000
    info = {
        "bloco_mm": [COMP, LARG, ALT], "fundo_mm": FUNDO, "parede_entre_fendas_mm": PAREDE, "camada_mm": a.camada,
        "elementos": feitos,
        "malha": {"faces": int(len(m.faces)), "corpos": int(m.body_count), "volume_valido": bool(m.is_volume),
                  "dimensoes_mm": [round(float(x), 2) for x in m.extents]},
        "PLA_cm3": round(vol, 2), "PLA_g": round(vol * 1.24, 1),
        "cera_para_encher_mL": round(COMP * LARG * ALT / 1000 - vol, 1),
    }
    with open(os.path.join(a.saida, "info.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in info.items() if k != "elementos"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
