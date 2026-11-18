"""Refaz só o visualizador (modelo.glb e visualizar.html) de pastas já geradas, com o contorno do tecido
adiposo em prisma de parede reta (igual ao recipiente impresso). Não recalcula placas nem STL.

Motivo: nas versões antigas a "casca" translúcida usava um perfil semicircular que afinava perto da pele, e as
placas (empilhadas na vertical) pareciam passar do limite da mama no visualizador, embora estivessem dentro do
recipiente impresso. Rodadas novas de mama3d.py já saem corretas.

Uso: python src/regerar_visualizador.py <pasta-da-imagem> [<pasta> ...]   (a pasta que contém info.json)
     python src/regerar_visualizador.py --tudo saida                       (todas as pastas com info.json)
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import trimesh
from skimage import transform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mama3d  # noqa: E402


def refazer(pasta: str) -> bool:
    info_p, npz_p, glb_p = (os.path.join(pasta, n) for n in ("info.json", "intermediarios.npz", "modelo.glb"))
    if not all(os.path.exists(p) for p in (info_p, npz_p, glb_p)):
        return False
    info = json.load(open(info_p, encoding="utf-8"))
    mask = np.load(npz_p)["mask"]
    v, H, s = float(info["voxel_mm"]), float(info["espessura_H_mm"]), float(info.get("escala_impressao", 1.0))
    vc = max(v, 1.0)
    mask_c = transform.rescale(mask.astype(np.float32), v / vc, order=1) > 0.5 if vc != v else mask
    casca = mama3d.reduzir(mama3d.malha_solida(np.where(mask_c, H, 0.0), mask_c, vc), 30000)
    if s != 1.0:
        casca.apply_scale(s)
    casca.visual.face_colors = mama3d.cor_parte("adiposo", 2)
    cena = trimesh.load(glb_p, force="scene")
    if "adiposo" in cena.geometry:
        cena.delete_geometry("adiposo")
    cena.add_geometry(casca, geom_name="adiposo", node_name="adiposo")
    mama3d.exportar(cena, pasta, os.path.basename(pasta))
    return True


def main():
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    pastas = (sorted(os.path.dirname(p) for p in glob.glob(os.path.join(a[1], "**", "info.json"), recursive=True))
              if a[0] == "--tudo" else a)
    feitos = 0
    for p in pastas:
        try:
            ok = refazer(p)
        except Exception as ex:  # noqa: BLE001
            print("ERRO", p, ex, flush=True)
            continue
        if ok:
            feitos += 1
            print("ok", p, flush=True)
    print(f"{feitos} visualizadores refeitos")


if __name__ == "__main__":
    main()
