"""Peça única vazada: o tecido fibroglandular inteiro numa só impressão; a gordura são os vazios.

Uso:
    python src/peca_unica.py saida/07_cera_abelha_1.0mm/lcc1/lcc1_pf06/lcc1
    python src/peca_unica.py <pasta com intermediarios.npz> --camada 0.2 --saida saida/09_peca_unica

Ideia (proposta da orientadora, 2026-11-18)
-------------------------------------------
Para o feixe primário só importa a espessura de material atravessada por cada raio
(Beer–Lambert: I = I0·exp(−Σ μ_i t_i), independente da ordem ao longo do raio). Em vez de
fatiar o mapa d(x,y) em ~22 placas de 1 mm, cada coluna (x,y) recebe uma coluna de PLA de
altura d(x,y), arredondada para a altura de camada da impressora (0,2 mm), numa peça só:

  * o erro de quantização cai de ≈ e/√12 (placas de espessura e) para ≈ camada/√12;
  * não há interfaces entre placas (lâminas de ar), nem empilhamento ou registro;
  * os vazios (gordura) ficam todos abertos para cima: a cera derretida entra e o ar sai
    sem canal de respiro — verificado aqui por rotulagem 3D dos vazios;
  * imprime deitada, sem suporte: cada camada está sobre a anterior (degraus só sobem).

Ilhas desconectadas (fibras soltas) podem ser unidas por pontes finas de PLA (--pontes),
ao custo de um excesso local de material, reportado no info.json.

Limitações: o mamograma 2D não informa a profundidade das fibras; a distribuição em z
aqui é uma escolha (tudo apoiado no fundo), que preserva a projeção primária, não a anatomia.
O espalhamento depende um pouco dessa distribuição (hipótese a validar com degrau radiografado).
Risco térmico: cera de abelha funde a ≈ 62–64 °C e o PLA amolece perto de 55–60 °C (Tg).
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
from scipy import ndimage as ndi
from skimage import draw, morphology

sys.path.insert(0, os.path.dirname(__file__))
from mama3d import limpar_regiao, reduzir, solido_entre  # noqa: E402

DENS_PLA = 1.24  # g/cm³ (ficha técnica do filamento do usuário)


def reamostrar_media(a: np.ndarray, f: int) -> np.ndarray:
    """Média em blocos f×f (preserva a espessura média projetada)."""
    if f == 1:
        return a.astype(np.float64)
    h, w = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    return a[:h, :w].reshape(h // f, f, w // f, f).mean(axis=(1, 3))


def niveis_relevo(d: np.ndarray, mask: np.ndarray, v: float, camada: float,
                  traco_min: float, area_min: float) -> tuple[np.ndarray, int]:
    """Altura impressa h(x,y), múltiplo da camada. Cada nível k é a região d ≥ (k−½)·camada,
    limpa do que a impressora não reproduz e contida no nível de baixo (sem balanço)."""
    n = int(np.ceil((d.max() - 0.5 * camada) / camada)) if d.max() > 0.5 * camada else 0
    h = np.zeros_like(d)
    anterior = mask
    for k in range(1, n + 1):
        # só o 1º nível exige área mínima (ilha manuseável); acima, basta o traço mínimo
        R = limpar_regiao(d >= (k - 0.5) * camada, v, traco_min, area_min if k == 1 else 1.0) & anterior
        if not R.any():
            break
        h[R] += camada
        anterior = R
    return h, int(round(h.max() / camada))


def pontes(base: np.ndarray, mask: np.ndarray, v: float, largura: float) -> tuple[np.ndarray, int]:
    """Une todas as ilhas de `base` à maior por segmentos retos (o mais curto a cada passo)."""
    lab, n = ndi.label(base)
    if n <= 1:
        return np.zeros_like(base), 0
    area = np.bincount(lab.ravel())
    area[0] = 0
    principal = lab == int(area.argmax())
    soltas = set(range(1, n + 1)) - {int(area.argmax())}
    ponte = np.zeros_like(base)
    r = max(1, int(round(largura / v / 2)))
    while soltas:
        dist, (iy, ix) = ndi.distance_transform_edt(~principal, return_indices=True)
        melhor, dmin = None, np.inf
        for i in soltas:
            ys, xs = np.nonzero(lab == i)
            j = int(dist[ys, xs].argmin())
            if dist[ys[j], xs[j]] < dmin:
                dmin, melhor = dist[ys[j], xs[j]], (i, ys[j], xs[j])
        i, y0, x0 = melhor
        y1, x1 = iy[y0, x0], ix[y0, x0]
        rr, cc = draw.line(int(y0), int(x0), int(y1), int(x1))
        seg = np.zeros_like(base)
        seg[rr, cc] = True
        seg = ndi.binary_dilation(seg, morphology.disk(r)) & mask
        ponte |= seg & ~base
        principal |= (lab == i) | seg
        soltas.discard(i)
    return ponte, n - 1


def vazios_presos(h: np.ndarray, mask: np.ndarray, camada: float) -> tuple[int, float]:
    """Rotula em 3D o espaço de cera (dentro do prisma da mama, acima e entre o PLA) e conta
    os componentes que não tocam o topo nem o lado do prisma: esses prenderiam ar."""
    nz = int(round(h.max() / camada)) + 1
    z = (np.arange(nz) + 0.5) * camada
    solido = h[..., None] > z[None, None, :]
    vazio = mask[..., None] & ~solido
    lab, n = ndi.label(vazio)  # conectividade 6 (faces)
    if n == 0:
        return 0, 0.0
    abertos = set(np.unique(lab[..., -1])) - {0}  # topo
    fora = ~ndi.binary_erosion(mask, border_value=0)  # coluna da borda do prisma
    abertos |= set(np.unique(lab[fora])) - {0}
    presos = [i for i in range(1, n + 1) if i not in abertos]
    vol = float(np.isin(lab, presos).sum()) if presos else 0.0
    return len(presos), vol


def processar(pasta_in: str, a) -> dict:
    z = np.load(os.path.join(pasta_in, "intermediarios.npz"))
    v0 = float(z["voxel"])
    f = max(1, int(round(a.voxel_malha / v0)))
    v = v0 * f
    d = reamostrar_media(z["d"], f)
    mask = reamostrar_media(z["mask"].astype(np.float32), f) > 0.5
    d[~mask] = 0
    info_in = {}
    p_info = os.path.join(pasta_in, "info.json")
    if os.path.exists(p_info):
        with open(p_info, encoding="utf-8") as fh:
            info_in = json.load(fh)
    fator = float(info_in.get("fator_material", 1.0))

    h, n_camadas = niveis_relevo(d, mask, v, a.camada, a.traco_min, a.area_min)
    base = h > 0
    ilhas = int(ndi.label(base)[1])
    ponte, n_pontes = (pontes(base, mask, v, a.largura_ponte) if a.pontes else (np.zeros_like(base), 0))
    h_pontes = np.round(a.altura_ponte / a.camada) * a.camada
    h = np.where(ponte, np.maximum(h, h_pontes), h)

    presos_n, presos_vox = vazios_presos(h, mask, a.camada)

    nome = os.path.basename(os.path.normpath(pasta_in))
    pasta = os.path.join(a.saida, os.path.basename(os.path.dirname(os.path.normpath(pasta_in))))  # p.ex. lcc1_pf06
    os.makedirs(pasta, exist_ok=True)

    reg = h > 0
    m = solido_entre(np.zeros_like(h), h - a.camada, reg, v, dz_mm=a.camada, suavizar=(0.7, 0.7, 0))
    m.apply_translation([0, 0, -m.bounds[0][2]])
    m = reduzir(m, a.max_faces)
    m.export(os.path.join(pasta, "peca_unica.stl"))
    # cópia para o Fusion (documento Y-up): espessura em Y, base em y = 0, centrada em X e Z
    yup = m.copy()
    yup.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    c = yup.bounds.mean(axis=0)
    yup.apply_translation([-c[0], -yup.bounds[0][1], -c[2]])
    yup.export(os.path.join(pasta, "peca_unica_yup.stl"))

    # erro da projeção (mm de PLA e mm de tecido equivalente) contra as placas da etapa de origem
    err = np.abs(h - d)[mask]
    err_sem_pontes = np.abs(np.where(ponte, 0, h) - d)[mask & ~ponte]
    dq = reamostrar_media(z["dq"], f)
    err_placas = np.abs(dq - d)[mask]
    vol_cm3 = float(m.volume) / 1000
    area_mask_cm2 = float(mask.sum()) * v**2 / 100
    H = float(z["t"].max())

    fig, ax = plt.subplots(1, 4, figsize=(16, 4.8))
    vmax = max(d.max(), 1e-6)
    i0 = ax[0].imshow(d, cmap="viridis", vmin=0, vmax=vmax); ax[0].set_title("d(x,y) desejado [mm PLA]", fontsize=9)
    fig.colorbar(i0, ax=ax[0], fraction=0.046)
    i1 = ax[1].imshow(np.where(mask, h, np.nan), cmap="viridis", vmin=0, vmax=vmax)
    ax[1].set_title(f"peça única: {n_camadas} camadas × {a.camada:g} mm", fontsize=9)
    fig.colorbar(i1, ax=ax[1], fraction=0.046)
    ax[1].contour(ponte, levels=[0.5], colors="r", linewidths=0.5)
    lim = max(float(np.percentile(np.abs(dq - d)[mask], 99)), 0.2)
    i2 = ax[2].imshow(np.where(mask, h - d, np.nan), cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax[2].set_title(f"erro peça única: média |·| = {err.mean():.3f} mm", fontsize=9)
    fig.colorbar(i2, ax=ax[2], fraction=0.046)
    i3 = ax[3].imshow(np.where(mask, dq - d, np.nan), cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax[3].set_title(f"erro placas de {float(z['e']):g} mm: média |·| = {err_placas.mean():.3f} mm", fontsize=9)
    fig.colorbar(i3, ax=ax[3], fraction=0.046)
    for x in ax:
        x.axis("off")
    fig.suptitle(f"{nome} · peça única vazada (vermelho = pontes)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "painel.png"), dpi=110)
    plt.close(fig)

    info = {
        "origem": pasta_in.replace("\\", "/"),
        "voxel_malha_mm": v, "camada_mm": a.camada, "fator_material": fator,
        "traco_min_mm": a.traco_min, "area_min_mm2": a.area_min,
        "n_camadas": n_camadas, "altura_peca_mm": round(float(h.max()), 2),
        "ilhas_antes_pontes": ilhas, "pontes": n_pontes,
        "pontes_altura_mm": float(h_pontes) if a.pontes else 0.0, "pontes_largura_mm": a.largura_ponte,
        "pontes_area_cm2": round(float(ponte.sum()) * v**2 / 100, 2),
        "vazios_presos": presos_n, "vazios_presos_mm3": round(presos_vox * v * v * a.camada, 2),
        "erro_medio_mm_PLA": round(float(err.mean()), 3),
        "erro_medio_mm_PLA_fora_pontes": round(float(err_sem_pontes.mean()), 3),
        "erro_medio_mm_tecido": round(float(err.mean() / fator), 3),
        "erro_medio_placas_origem_mm_PLA": round(float(err_placas.mean()), 3),
        "espessura_placas_origem_mm": float(z["e"]),
        "fracao_glandular_impressa": round(float((h / fator).sum() / z["t"].sum() * f * f), 3) if fator else None,
        "malha": {"faces": int(len(m.faces)), "corpos": int(m.body_count), "volume_valido": bool(m.is_volume),
                  "dimensoes_mm": [round(float(x), 2) for x in m.extents]},
        "PLA_cm3": round(vol_cm3, 2), "PLA_g": round(vol_cm3 * DENS_PLA, 1),
        "cera_estimada_mL": round(area_mask_cm2 * H / 10 - vol_cm3, 1), "altura_compressao_mm": H,
        "avisos": [
            "distribuição em z escolhida (fibroso apoiado no fundo): preserva a projeção primária, não a anatomia",
            "cera de abelha funde a ~62–64 °C e o PLA amolece perto de 55–60 °C: testar a infiltração num cubo antes",
            "a contração da cera ao solidificar pode abrir frestas: verter em etapas e resfriar devagar",
        ],
    }
    with open(os.path.join(pasta, "info.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, ensure_ascii=False)
    info["pasta"] = pasta
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entradas", nargs="+", help="pastas com intermediarios.npz (saída do modo placas)")
    ap.add_argument("--camada", type=float, default=0.2, help="altura de camada da impressora em mm (padrão 0.2)")
    ap.add_argument("--voxel-malha", dest="voxel_malha", type=float, default=0.5,
                    help="resolução lateral da malha em mm (padrão 0.5; bico de 0,4 mm)")
    ap.add_argument("--traco-min", dest="traco_min", type=float, default=1.2, help="traço mínimo impresso em mm")
    ap.add_argument("--area-min", dest="area_min", type=float, default=25.0, help="área mínima de ilha em mm²")
    ap.add_argument("--sem-pontes", dest="pontes", action="store_false", help="não unir as ilhas")
    ap.add_argument("--largura-ponte", dest="largura_ponte", type=float, default=1.5, help="mm (padrão 1.5)")
    ap.add_argument("--altura-ponte", dest="altura_ponte", type=float, default=0.6, help="mm (padrão 0.6)")
    ap.add_argument("--max-faces", dest="max_faces", type=int, default=None)
    ap.add_argument("--saida", default=os.path.join("saida", "09_peca_unica"))
    a = ap.parse_args(argv)
    for p in a.entradas:
        info = processar(p, a)
        print(json.dumps({k: info[k] for k in info if k != "avisos"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
