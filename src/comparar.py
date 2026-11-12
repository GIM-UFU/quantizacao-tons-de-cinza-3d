"""Compara o modelo 3D (o que está em visualizar.html) com a mamografia original.

Ideia: a mamografia é uma PROJEÇÃO. Então projetamos o modelo de volta em 2D — a
soma das espessuras das placas ao longo de z dá um mapa "espessura de material
denso por pixel" — e medimos o quanto esse mapa concorda com a imagem original.

Este script:
  1. extrai o modelo (GLB em base64) de dentro do visualizar.html — ou lê os STL
     das placas com --stl (resolução total);
  2. rasteriza cada placa e soma as espessuras → comparacao/projecao_modelo.tif
     (32 bits, em mm, mesma orientação da imagem original, calibrado para o ImageJ);
  3. calcula as métricas em Python (referência) → comparacao/resultado_python.json;
  4. deixa pronto o comando do ImageJ, que repete a comparação de forma independente
     com a macro src/comparar_imagej.ijm.

Uso: python src/comparar.py saida/placas/lcc1 [--stl]
"""
from __future__ import annotations

import argparse
import base64
import glob
import io
import json
import os
import re
import sys

import numpy as np
import tifffile
import trimesh
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from skimage import filters, transform

sys.path.insert(0, os.path.dirname(__file__))
from dicom_io import carregar  # noqa: E402
from mama3d import linearizar  # noqa: E402


def placas_do_html(caminho_html: str) -> dict[str, trimesh.Trimesh]:
    with open(caminho_html, encoding="utf-8") as f:
        html = f.read()
    m = re.search(r'const b64="([A-Za-z0-9+/=]+)"', html)
    if not m:
        raise RuntimeError("GLB embutido não encontrado no HTML.")
    cena = trimesh.load(io.BytesIO(base64.b64decode(m.group(1))), file_type="glb")
    return {k: g for k, g in cena.geometry.items() if k.startswith("placa_")}


def placas_dos_stl(pasta: str) -> dict[str, trimesh.Trimesh]:
    return {os.path.splitext(os.path.basename(f))[0]: trimesh.load(f)
            for f in sorted(glob.glob(os.path.join(pasta, "placas", "placa_*.stl")))}


def rasterizar(placa: trimesh.Trimesh, forma, voxel_mm: float, sub: int = 2) -> np.ndarray:
    """Cobertura (0–1) da placa em cada pixel. Uma placa é uma extrusão: sua projeção
    em z é a união dos triângulos da face de cima (normal apontando para +z).
    `sub` = superamostragem, para a borda não ficar serrilhada."""
    h, w = forma
    cima = placa.face_normals[:, 2] > 0.5
    tri = placa.triangles[cima][:, :, :2]
    # malha: x = coluna·voxel, y = −linha·voxel  (convenção de mama3d.superficie)
    px = np.stack([tri[:, :, 0] / voxel_mm + 0.5, -tri[:, :, 1] / voxel_mm + 0.5], axis=2) * sub
    im = Image.new("L", (w * sub, h * sub), 0)
    d = ImageDraw.Draw(im)
    for t in px:
        d.polygon([tuple(p) for p in t], fill=255)
    a = np.asarray(im, dtype=np.float32) / 255.0
    return a.reshape(h, sub, w, sub).mean((1, 3))


def metricas(orig: np.ndarray, esp: np.ndarray, mask: np.ndarray, interior: np.ndarray) -> dict:
    """Mesmas definições da macro do ImageJ (comparar_imagej.ijm)."""
    def bloco(reg):
        o, e = orig[reg], esp[reg]
        lim = float(filters.threshold_otsu(orig[mask]))  # gordura × denso, na mama inteira
        img_denso, mod_denso = o > lim, e > 0
        tp = int((img_denso & mod_denso).sum()); fp = int((~img_denso & mod_denso).sum())
        fn = int((img_denso & ~mod_denso).sum()); tn = int((~img_denso & ~mod_denso).sum())
        r = float(np.corrcoef(o, e)[0, 1])
        b, a0 = np.polyfit(e, o, 1)
        return {
            "pixels": int(reg.sum()), "limiar_otsu": round(lim, 2),
            "pearson_r": round(r, 4), "r2": round(r * r, 4),
            "ajuste_cinza_por_mm": round(float(b), 3), "ajuste_intercepto": round(float(a0), 2),
            "rmse_cinza": round(float(np.sqrt(np.mean((a0 + b * e - o) ** 2))), 2),
            "acerto_pixel_a_pixel": round((tp + tn) / max(tp + tn + fp + fn, 1), 4),
            "dice": round(2 * tp / max(2 * tp + fp + fn, 1), 4),
            "iou": round(tp / max(tp + fp + fn, 1), 4),
            "sensibilidade": round(tp / max(tp + fn, 1), 4),
            "precisao": round(tp / max(tp + fp, 1), 4),
            "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        }
    return {"mama_inteira": bloco(mask), "interior_sem_faixa_da_pele": bloco(interior)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta", help="pasta de saída do modo placas (contém visualizar.html e info.json)")
    ap.add_argument("--stl", action="store_true", help="usar os STL das placas em vez do GLB do HTML")
    a = ap.parse_args()

    with open(os.path.join(a.pasta, "info.json"), encoding="utf-8") as f:
        info = json.load(f)
    v, e = info["voxel_mm"], info["espessura_placa_mm"]
    m = carregar(info["arquivo"], info["pixel_mm"][0])
    forma = tuple(int(round(s * p / v)) for s, p in zip(m.imagem.shape, m.pixel_mm))

    placas = placas_dos_stl(a.pasta) if a.stl else placas_do_html(os.path.join(a.pasta, "visualizar.html"))
    if not placas:
        sys.exit("Nenhuma placa encontrada.")
    esp = np.zeros(forma, np.float32)
    for nome, p in sorted(placas.items()):
        esp += rasterizar(p, forma, v) * float(p.extents[2])
    if info.get("espelhado"):
        esp = esp[:, ::-1]  # volta para a orientação da imagem original

    out = os.path.join(a.pasta, "comparacao")
    os.makedirs(out, exist_ok=True)
    proj = os.path.join(out, "projecao_modelo.tif")
    tifffile.imwrite(proj, esp, imagej=True, resolution=(1 / v, 1 / v), metadata={"unit": "mm"})

    # referência em Python: mesma grade, mesmas definições da macro
    orig = transform.resize(m.imagem, forma, order=1, anti_aliasing=True, preserve_range=True).astype(np.float32)
    orig, _ = linearizar(orig)
    mask = ndi.binary_fill_holes(orig > filters.threshold_triangle(orig))
    lab, n = ndi.label(mask)
    if n > 1:
        mask = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
    npz = os.path.join(a.pasta, "intermediarios.npz")
    if os.path.exists(npz):  # usa a MESMA máscara do gerador (MLO: sem peitoral nem abdome)
        mk = np.load(npz)["mask"]
        mk = mk[:, ::-1] if info.get("espelhado") else mk
        if mk.shape == mask.shape:
            mask = mk
    margem_px = int(round(info.get("margem_pele_mm", 0) / v))
    interior = ndi.distance_transform_edt(mask) > margem_px if margem_px else mask
    res = {"fonte_do_modelo": "STL das placas" if a.stl else "GLB embutido no visualizar.html",
           "placas": len(placas), "espessura_placa_mm": e, "grade": list(forma), "voxel_mm": v,
           "espessura_max_projetada_mm": round(float(esp.max()), 2),
           "volume_projetado_cm3": round(float(esp.sum() * v * v / 1000), 1),
           "margem_pele_px": margem_px} | metricas(orig, esp, mask, interior)
    with open(os.path.join(out, "resultado_python.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    # mapa de acerto (mesmas cores da macro): verde = denso nos dois, vermelho = só no
    # modelo, azul = só na imagem
    lim = filters.threshold_otsu(orig[mask])
    g = np.clip(140 * (orig - orig[mask].min()) / (np.ptp(orig[mask]) + 1e-9), 0, 140).astype(np.uint8)
    rgb = np.stack([g, g, g], axis=2)
    di, dm = (orig > lim) & mask, (esp > 0) & mask
    rgb[di & dm] = (40, 200, 60)
    rgb[dm & ~di] = (230, 50, 50)
    rgb[di & ~dm] = (60, 110, 255)
    Image.fromarray(rgb).save(os.path.join(out, "mapa_acerto_python.png"))

    macro = os.path.join(os.path.dirname(os.path.abspath(__file__)), "comparar_imagej.ijm")
    arg = "|".join([os.path.abspath(info["arquivo"]), os.path.abspath(proj), os.path.abspath(out) + os.sep,
                    str(margem_px)])
    with open(os.path.join(out, "comando_imagej.txt"), "w", encoding="utf-8") as f:
        # ImageJ 1.x em lote. Não usar -Djava.awt.headless=true: o open() falha sem AWT.
        ij = os.path.join(os.path.dirname(os.path.dirname(macro)), "externo", "imagej", "ImageJ")
        f.write(f'"{ij}\\jre\\bin\\java.exe" -Xmx3g -jar "{ij}\\ij.jar" -ijpath "{ij}" -batch "{macro}" "{arg}"\n')

    for k in ("mama_inteira", "interior_sem_faixa_da_pele"):
        r = res[k]
        print(f"▸ {k}: r={r['pearson_r']} · R²={r['r2']} · acerto={r['acerto_pixel_a_pixel']:.1%} · "
              f"Dice={r['dice']} · sens={r['sensibilidade']:.1%} · prec={r['precisao']:.1%}")
    print(f"fonte: {res['fonte_do_modelo']} · volume projetado {res['volume_projetado_cm3']} cm³")
    print("saídas em:", out)


if __name__ == "__main__":
    main()
