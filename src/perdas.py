"""Quanto da variação de cinza convertida em relevo se perde em cada etapa até o STL.

Cadeia de conversão (modo placas), do cinza ao modelo impresso:

   cinza ──(corte gordura × fibroso)──► d_ref   relevo "ideal": espessura de tecido denso
                                                que a imagem pede, pixel a pixel (mm)
   d_ref ─(faixa da pele)────────────► d_pele   pele = parede do recipiente, sem placa
         ─(suavização 0,75 mm)───────► d        o que o gerador de fato quantiza
         ─(quantização em placas)────► dq_q     degraus de `e` mm
         ─(limpeza de manufatura)────► dq       traço mínimo, ilhas e furos pequenos
         ─(malha: voxel + STL)───────► modelo   projeção das placas do STL (o que será impresso)

Para cada etapa, contra d_ref, dentro da mama:
   volume retido  = Σ etapa / Σ d_ref
   variação preservada (R²) = 1 − Σ(etapa − d_ref)² / Σ(d_ref − média)²
        (100 % = relevo idêntico ao ideal; 0 % = tão ruim quanto um relevo constante na média)
   RMSE e erro absoluto médio, em mm de tecido denso.

Antes da cadeia, o corte gordura × fibroso é uma DECISÃO de projeto (gordura = parafina,
não é impressa), e não perda de conversão; mesmo assim reportamos quanto do sinal de
cinza fica abaixo do corte, para a decisão ficar visível.

Uso: python src/perdas.py saida/placas/lcc1
Requer a pasta gerada por mama3d.py no modo placas (intermediarios.npz + placas/*.stl).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(__file__))
from comparar import rasterizar  # noqa: E402


def projetar_stl(pasta: str, forma, voxel: float, escala: float) -> np.ndarray:
    """Soma, em mm de MODELO, a espessura das placas do STL projetadas em z."""
    arqs = sorted(glob.glob(os.path.join(pasta, "placas", "placa_*.stl")))
    if not arqs:
        sys.exit("Nenhuma placa em placas/.")
    esp = np.zeros(forma, np.float32)
    for f in arqs:
        p = trimesh.load(f)
        # x_impresso = x_modelo · escala → coordenadas em unidades de (voxel·escala)
        esp += rasterizar(p, forma, voxel * escala) * float(p.extents[2]) / escala
    return esp


def etapa(nome: str, X: np.ndarray, ref: np.ndarray, mask: np.ndarray, v: float) -> dict:
    x, r = X[mask].astype(np.float64), ref[mask].astype(np.float64)
    sst = float(((r - r.mean()) ** 2).sum())
    sse = float(((x - r) ** 2).sum())
    return {"etapa": nome, "volume_cm3": round(float(x.sum() * v * v / 1000), 2),
            "volume_retido": round(float(x.sum() / r.sum()), 4),
            "variacao_preservada_R2": round(1 - sse / sst, 4),
            "rmse_mm": round(float(np.sqrt(np.mean((x - r) ** 2))), 3),
            "erro_abs_medio_mm": round(float(np.mean(np.abs(x - r))), 3),
            "correlacao": round(float(np.corrcoef(x, r)[0, 1]), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta")
    a = ap.parse_args()

    with open(os.path.join(a.pasta, "info.json"), encoding="utf-8") as f:
        info = json.load(f)
    Z = np.load(os.path.join(a.pasta, "intermediarios.npz"))
    img, mask, t = Z["img"], Z["mask"], Z["t"]
    d_ref, d, dq_q, dq = Z["d_ref"], Z["d"], Z["dq_sem_limpeza"], Z["dq"]
    v, e, s = float(Z["voxel"]), float(Z["e"]), float(Z["escala"])
    piso, zero, um = float(Z["piso"]), float(Z["zero"]), float(Z["um"])

    # ── antes da cadeia: quanto do sinal de cinza fica abaixo do corte (gordura) ──
    g = np.clip((img - piso) / (um - piso + 1e-9), 0, 1)
    gm, abaixo = g[mask], (img[mask] <= zero)
    sst_g = float(((gm - gm.mean()) ** 2).sum())
    corte = {"sinal_de_cinza_abaixo_do_corte": round(float((gm * t[mask])[abaixo].sum() / (gm * t[mask]).sum()), 4),
             "variancia_de_cinza_abaixo_do_corte": round(float(((gm[abaixo] - gm.mean()) ** 2).sum() / sst_g), 4),
             "pixels_da_mama_abaixo_do_corte": round(float(abaixo.mean()), 4)}

    # ── cadeia ──
    margem = float(info.get("margem_pele_mm", 0))
    d_pele = d_ref * (ndi.distance_transform_edt(mask) * v >= margem) if margem > 0 else d_ref
    modelo = projetar_stl(a.pasta, img.shape, v, s)
    etapas = [etapa("d_ref (relevo ideal)", d_ref, d_ref, mask, v),
              etapa("faixa da pele", d_pele, d_ref, mask, v),
              etapa("suavização 0,75 mm", d, d_ref, mask, v),
              etapa("quantização em placas", dq_q, d_ref, mask, v),
              etapa("limpeza de manufatura", dq, d_ref, mask, v),
              etapa("malha / STL final", modelo, d_ref, mask, v)]
    for i in range(1, len(etapas)):
        etapas[i]["perda_nesta_etapa_R2"] = round(etapas[i - 1]["variacao_preservada_R2"] - etapas[i]["variacao_preservada_R2"], 4)
        etapas[i]["delta_volume_nesta_etapa"] = round(etapas[i]["volume_retido"] - etapas[i - 1]["volume_retido"], 4)

    # custo de usar placas de espessura e: erro de quantização esperado ≈ e/√12 (uniforme)
    res = {"pasta": a.pasta, "escala_impressao": s, "espessura_placa_modelo_mm": e,
           "erro_quantizacao_teorico_mm": round(e / np.sqrt(12), 3),
           "corte_gordura_fibroso": corte, "etapas": etapas,
           "resumo": {"variacao_preservada_no_STL": etapas[-1]["variacao_preservada_R2"],
                      "volume_retido_no_STL": etapas[-1]["volume_retido"],
                      "variacao_perdida_no_STL": round(1 - etapas[-1]["variacao_preservada_R2"], 4)}}
    with open(os.path.join(a.pasta, "resultado_perdas.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    print(f"escala {s:g} · placas de {e:g} mm (modelo) · erro de quantização teórico ±{res['erro_quantizacao_teorico_mm']} mm RMS\n")
    print(f"Corte gordura×fibroso (decisão de projeto): {corte['pixels_da_mama_abaixo_do_corte']:.1%} dos pixels da mama, "
          f"{corte['sinal_de_cinza_abaixo_do_corte']:.1%} do sinal de cinza e "
          f"{corte['variancia_de_cinza_abaixo_do_corte']:.1%} da variância do cinza ficam abaixo do corte.\n")
    print(f"{'etapa':<26}{'volume':>9}{'retido':>9}{'preserv.':>10}{'perda':>8}{'RMSE':>9}{'|erro|':>9}")
    for x in etapas:
        print(f"{x['etapa']:<26}{x['volume_cm3']:>8.1f}c{x['volume_retido']:>9.1%}{x['variacao_preservada_R2']:>10.1%}"
              f"{x.get('perda_nesta_etapa_R2', 0):>8.1%}{x['rmse_mm']:>8.2f}m{x['erro_abs_medio_mm']:>8.2f}m")

    # ── figura ──
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))
    vmax = float(d_ref.max())
    ax[0].imshow(d_ref, cmap="viridis", vmin=0, vmax=vmax); ax[0].set_title("d_ref: relevo ideal [mm]", fontsize=9)
    ax[1].imshow(modelo, cmap="viridis", vmin=0, vmax=vmax); ax[1].set_title("modelo (projeção do STL) [mm]", fontsize=9)
    err = (modelo - d_ref) * mask
    lim = float(np.percentile(np.abs(err[mask]), 99))
    h = ax[2].imshow(err, cmap="RdBu_r", vmin=-lim, vmax=lim); ax[2].set_title("erro = modelo − ideal [mm]", fontsize=9)
    fig.colorbar(h, ax=ax[2], fraction=0.046)
    for x in ax[:3]:
        x.axis("off")
    nomes = [x["etapa"].replace(" (relevo ideal)", "") for x in etapas]
    r2 = [100 * x["variacao_preservada_R2"] for x in etapas]
    vr = [100 * x["volume_retido"] for x in etapas]
    xs = np.arange(len(etapas))
    ax[3].bar(xs - 0.2, r2, 0.4, label="variação preservada (R²)", color="#2a9d8f")
    ax[3].bar(xs + 0.2, vr, 0.4, label="volume retido", color="#e9c46a")
    for i in xs:
        ax[3].text(i - 0.2, r2[i] + 1, f"{r2[i]:.0f}", ha="center", fontsize=7)
        ax[3].text(i + 0.2, vr[i] + 1, f"{vr[i]:.0f}", ha="center", fontsize=7)
    ax[3].set_xticks(xs); ax[3].set_xticklabels(nomes, rotation=35, ha="right", fontsize=8)
    ax[3].set_ylim(0, 112); ax[3].set_ylabel("% do relevo ideal"); ax[3].legend(fontsize=8, loc="lower left")
    ax[3].set_title("o que sobra após cada etapa", fontsize=9)
    fig.suptitle(os.path.basename(os.path.dirname(a.pasta + os.sep)) + f" · escala {s:g}", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(a.pasta, "perdas.png"), dpi=110)
    plt.close(fig)
    print("\nsaídas:", os.path.join(a.pasta, "resultado_perdas.json"), "·", os.path.join(a.pasta, "perdas.png"))


if __name__ == "__main__":
    main()
