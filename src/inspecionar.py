"""Inspeciona DICOMs de mamografia: metadados relevantes + miniatura PNG.

Uso:  python src/inspecionar.py dados/  [--saida saida/inspecao]

O que olhar no relatório:
- Quadros > 1  → é tomossíntese (DBT): há volume real, use o modo 'volume'.
- Intencao     → FOR PROCESSING (cru, ~linear com a atenuação) é o ideal para
                 qualquer estimativa física; FOR PRESENTATION já passou por
                 processamento do fabricante (não-linear).
- EspessuraMama_mm → necessário para dar escala física ao eixo z.
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

sys.path.insert(0, os.path.dirname(__file__))
from dicom_io import carregar, listar_dicoms  # noqa: E402

CAMPOS = ["Formato", "Vista", "Lateralidade", "Quadros", "Linhas", "Colunas", "Intencao",
          "Fotometria", "BitsArmazenados", "EspessuraMama_mm", "kVp", "Fabricante", "FontePixel"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entradas", nargs="+")
    ap.add_argument("--saida", default="saida/inspecao")
    ap.add_argument("--pixel-mm", dest="pixel_mm", type=float)
    a = ap.parse_args()

    arquivos = listar_dicoms(a.entradas)
    if not arquivos:
        sys.exit("Nenhum DICOM/TIFF encontrado.")
    os.makedirs(a.saida, exist_ok=True)

    relatorio = {}
    for arq in arquivos:
        nome = os.path.splitext(os.path.basename(arq))[0]
        try:
            m = carregar(arq, a.pixel_mm)
        except Exception as e:  # ex.: sintaxe de transferência sem decodificador
            print(f"[ERRO] {nome}: {e}")
            continue
        relatorio[nome] = m.meta | {"pixel_mm": m.pixel_mm}

        print(f"\n▸ {nome}")
        for c in CAMPOS:
            print(f"   {c:<17}: {m.meta.get(c)}")
        print(f"   {'pixel_mm':<17}: {m.pixel_mm}")

        img = m.imagem[m.imagem.shape[0] // 2] if m.eh_volume else m.imagem
        lo, hi = np.percentile(img, [1, 99])
        plt.figure(figsize=(4, 5))
        plt.imshow(np.clip((img - lo) / (hi - lo + 1e-9), 0, 1), cmap="gray")
        plt.title(f"{nome}  {m.meta['Vista']} {m.meta['Lateralidade']}", fontsize=8)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(a.saida, f"{nome}.png"), dpi=120)
        plt.close()

    with open(os.path.join(a.saida, "metadados.json"), "w", encoding="utf-8") as f:
        json.dump(relatorio, f, indent=2, ensure_ascii=False)
    print(f"\nMiniaturas e metadados.json em: {a.saida}")


if __name__ == "__main__":
    main()
