"""Gráficos do erro de quantização: peça única (camadas de 0,2 mm) × placas de 1 mm e 2 mm.

Erro = espessura impressa − espessura desejada d(x,y) de PLA, por pixel dentro da mama, na mesma
grade de 0,5 mm da peça única. Em mm de tecido equivalente = erro de PLA / fator de material.

Uso:
    python src/graficos_erro.py                       # lcc1_pf06 → saida/ENCAMINHAR EXEMPLO/graficos_erro
    python src/graficos_erro.py --img lcc1 --pf pf06 --saida <pasta>
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from peca_unica import niveis_relevo, pontes, reamostrar_media  # noqa: E402

CORES = {"peça única 0,2 mm": "#2a7fb8", "placas 1 mm": "#e08a2c", "placas 2 mm": "#c43c3c"}


def carregar(pasta: str, f: int):
    z = np.load(os.path.join(pasta, "intermediarios.npz"))
    mask = reamostrar_media(z["mask"].astype(np.float32), f) > 0.5
    d = reamostrar_media(z["d"], f)
    dq = reamostrar_media(z["dq"], f)
    d[~mask] = 0
    return z, mask, d, dq


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--img", default="lcc1")
    ap.add_argument("--pf", default="pf06")
    ap.add_argument("--voxel-malha", dest="voxel_malha", type=float, default=0.5)
    ap.add_argument("--saida", default=os.path.join("saida", "ENCAMINHAR EXEMPLO", "graficos_erro"))
    a = ap.parse_args(argv)
    os.makedirs(a.saida, exist_ok=True)

    p1 = os.path.join("saida", "07_cera_abelha_1.0mm", a.img, f"{a.img}_{a.pf}", a.img)
    p2 = os.path.join("saida", "08_cera_abelha_2.0mm", a.img, f"{a.img}_{a.pf}", a.img)
    z = np.load(os.path.join(p1, "intermediarios.npz"))
    v0 = float(z["voxel"])
    f = max(1, int(round(a.voxel_malha / v0)))
    v = v0 * f
    _, mask, d, dq1 = carregar(p1, f)
    _, mask2, d2, dq2 = carregar(p2, f)
    assert np.allclose(d, d2, atol=1e-6), "d(x,y) difere entre 1 mm e 2 mm"
    fator = float(json.load(open(os.path.join(p1, "info.json"), encoding="utf-8")).get("fator_material", 1.0))

    # peça única com os parâmetros padrão de src/peca_unica.py
    h, ncam = niveis_relevo(d, mask, v, 0.2, 1.2, 25.0)
    ponte, _ = pontes(h > 0, mask, v, 1.5)
    h = np.where(ponte, np.maximum(h, 0.6), h)

    imp = {"peça única 0,2 mm": h, "placas 1 mm": dq1, "placas 2 mm": dq2}
    err = {k: (x - d)[mask] for k, x in imp.items()}
    res = {k: {"medio_abs_mm_PLA": float(np.abs(e).mean()), "rms_mm_PLA": float(np.sqrt((e ** 2).mean())),
               "p95_abs_mm_PLA": float(np.percentile(np.abs(e), 95)), "vies_mm_PLA": float(e.mean())}
           for k, e in err.items()}
    for r in res.values():
        r["medio_abs_mm_tecido"] = r["medio_abs_mm_PLA"] / fator
    sub = f"{a.img} · limiar {a.pf[2]},{a.pf[3]} · fator de material {fator:.3f} (PLA × cera, 20 keV)"
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    # 1) mapas do erro, mesma escala
    lim = 1.0
    fig, ax = plt.subplots(1, 4, figsize=(16, 5.4))
    i0 = ax[0].imshow(np.where(mask, d, np.nan), cmap="viridis")
    ax[0].set_title("espessura desejada de PLA\nd(x,y) [mm]")
    fig.colorbar(i0, ax=ax[0], fraction=0.046)
    for x, (k, hi) in zip(ax[1:], imp.items()):
        im = x.imshow(np.where(mask, hi - d, np.nan), cmap="RdBu_r", vmin=-lim, vmax=lim)
        x.set_title(f"{k}\nerro médio |·| = {res[k]['medio_abs_mm_PLA']:.3f} mm")
    fig.colorbar(im, ax=ax[1:], fraction=0.015, label="impresso − desejado [mm de PLA]")
    for x in ax:
        x.axis("off")
    fig.suptitle(f"Mapa do erro de quantização · {sub}", fontsize=11)
    fig.savefig(os.path.join(a.saida, "1_mapas_erro.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 2) distribuição do erro
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bins = np.linspace(-1.1, 1.1, 111)
    for k, e in err.items():
        ax.hist(e, bins=bins, density=True, histtype="stepfilled", alpha=0.35, color=CORES[k])
        ax.hist(e, bins=bins, density=True, histtype="step", lw=1.6, color=CORES[k],
                label=f"{k}: RMS {res[k]['rms_mm_PLA']:.3f} mm")
    ax.set_xlabel("erro (impresso − desejado) [mm de PLA]")
    ax.set_ylabel("densidade de pixels")
    ax.set_title(f"Distribuição do erro por pixel da mama\n{sub}", fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.saida, "2_histograma_erro.png"), dpi=140)
    plt.close(fig)

    # 3) erro × passo de quantização (arredondamento ideal) + teoria + pipelines reais
    passos = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
    dm = d[mask]
    mae = [np.abs(np.round(dm / s) * s - dm).mean() for s in passos]
    rms = [np.sqrt(((np.round(dm / s) * s - dm) ** 2).mean()) for s in passos]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    xs = np.linspace(0.05, 5.2, 200)
    ax.plot(xs, xs / np.sqrt(12), "k--", lw=1, label="teoria RMS = passo/√12")
    ax.plot(xs, xs / 4, "k:", lw=1, label="teoria média |·| = passo/4")
    ax.plot(passos, rms, "o-", color="#555", ms=4, label="arredondamento ideal: RMS")
    ax.plot(passos, mae, "s-", color="#999", ms=4, label="arredondamento ideal: média |·|")
    for k, s in zip(imp, [0.2, 1.0, 2.0]):
        ax.plot(s, res[k]["rms_mm_PLA"], "*", ms=15, color=CORES[k], label=f"{k} (gerado): RMS")
    ax.axvspan(4, 5, color="#c43c3c", alpha=0.07)
    ax.text(4.5, 0.05, "placas de\n4–5 mm", ha="center", fontsize=8, color="#c43c3c")
    ax.set_xlabel("passo de quantização (espessura da placa ou da camada) [mm]")
    ax.set_ylabel("erro [mm de PLA]")
    ax.set_title(f"Erro de quantização × passo\n{sub}", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(0, 5.3)
    ax.set_ylim(0, None)
    fig.tight_layout()
    fig.savefig(os.path.join(a.saida, "3_erro_vs_passo.png"), dpi=140)
    plt.close(fig)

    # 4) resumo em barras: mm de PLA e mm de tecido equivalente
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ks = list(imp)
    for x, campo, tit in ((ax[0], "medio_abs_mm_PLA", "mm de PLA"), (ax[1], "medio_abs_mm_tecido", "mm de tecido equivalente")):
        vals = [res[k][campo] for k in ks]
        b = x.bar(ks, vals, color=[CORES[k] for k in ks])
        x.bar_label(b, labels=[f"{y:.3f}" for y in vals], padding=3)
        x.set_title(f"erro médio |·| [{tit}]")
        x.set_ylim(0, max(vals) * 1.2)
    fig.suptitle(f"Resumo do erro · {sub}", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(a.saida, "4_resumo_erro.png"), dpi=140)
    plt.close(fig)

    out = {"imagem": a.img, "limiar": a.pf, "voxel_mm": v, "fator_material": fator, "camadas_peca": ncam,
           "pixels_mama": int(mask.sum()),
           "resultados": {k: {c: round(x, 4) for c, x in r.items()} for k, r in res.items()},
           "arredondamento_ideal": {f"{s:g}": {"medio_abs": round(m, 4), "rms": round(r, 4)}
                                    for s, m, r in zip(passos, mae, rms)},
           "nota": "fator de material monoenergético (20 keV); erro só de quantização da projeção, sem impressão real"}
    with open(os.path.join(a.saida, "erros.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print(json.dumps(out["resultados"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
