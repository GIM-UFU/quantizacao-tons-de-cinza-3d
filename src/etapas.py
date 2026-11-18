"""Gera e compara as ETAPAS da pesquisa (antes × depois da atenuação, várias espessuras de placa).

Organização de saida/ (cada etapa: <etapa>/<imagem>/<imagem>_pf06|pf08|pf10/<imagem>/ com STL, visualizador, info.json):
  01_antes_sem_atenuacao_3.4mm      fator de material 1,0 (sem atenuação), placas de 3,4 mm
  02_depois_com_atenuacao_3.4mm     fator por atenuação linear (PLA, 20 keV), 3,4 mm
  03_depois_com_atenuacao_1.0mm     idem, 1,0 mm
  04_depois_com_atenuacao_1.4mm     idem, 1,4 mm
  06_depois_com_atenuacao_2.0mm     idem, 2,0 mm (ressalva: placas mais robustas)
  07_cera_abelha_1.0mm              fator com CERA DE ABELHA no lugar da parafina (0,478), 1,0 mm
  08_cera_abelha_2.0mm              idem, 2,0 mm
  05_comparacoes/                   tabelas, gráficos e mostradores
  exploratorios/                    testes antigos e avulsos (nada foi apagado)
O sufixo pf06/pf08/pf10 é o limiar fibroso 0,6/0,8/1,0.

Subcomandos:
  python src/etapas.py gerar      cria as etapas 03, 04, 06, 07 e 08 que faltarem (8 imagens × 3 limiares)
  python src/etapas.py comparar   mede o erro de quantização das 4 etapas e escreve 05_comparacoes/quantizacao_etapas/

Erro de quantização (igual para todas as etapas, em mm de TECIDO equivalente):
  erro = | dq_sem_limpeza − d | / fator_material, na região da mama
(d = relevo suavizado em mm de placa; dq_sem_limpeza = relevo após virar degraus de espessura e;
 dividir pelo fator converte mm de placa em mm de tecido glandular).
"""
from __future__ import annotations

import csv
import glob
import json
import os
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAIDA = os.path.join(RAIZ, "saida")
ETAPAS = [("01_antes_sem_atenuacao_3.4mm", "antes · fator 1,0 · 3,4 mm", None),
          ("02_depois_com_atenuacao_3.4mm", "depois · atenuação · 3,4 mm", 3.4),
          ("03_depois_com_atenuacao_1.0mm", "depois · atenuação · 1,0 mm", 1.0),
          ("04_depois_com_atenuacao_1.4mm", "depois · atenuação · 1,4 mm", 1.4),
          ("06_depois_com_atenuacao_2.0mm", "depois · atenuação · 2,0 mm", 2.0),
          ("07_cera_abelha_1.0mm", "cera de abelha · 1,0 mm", 1.0),
          ("08_cera_abelha_2.0mm", "cera de abelha · 2,0 mm", 2.0)]
CERA = ("07", "08")  # etapas com cera de abelha no lugar da parafina (fator 0,478 em vez de 0,387)
LIMIARES = [("pf06", "0.6", "0,6"), ("pf08", "0.8", "0,8"), ("pf10", "1.0", "1,0")]
CMP = os.path.join(SAIDA, "05_comparacoes")


def imagens() -> list[dict]:
    out = {}
    for f in sorted(glob.glob(os.path.join(SAIDA, ETAPAS[0][0], "*", "*_pf10", "*", "info.json"))):
        d = json.load(open(f, encoding="utf-8"))
        nome = os.path.basename(os.path.dirname(f))
        out[nome] = {"nome": nome, "arquivo": os.path.join(RAIZ, d["arquivo"]), "pixel": d["pixel_mm"][0], "voxel": d["voxel_mm"]}
    return list(out.values())


def gerar():
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    for etapa, _, esp in ETAPAS:
        if esp is None or etapa.startswith("02"):
            continue
        for img in imagens():
            for suf, lim, _ in LIMIARES:
                saida = os.path.join(SAIDA, etapa, img["nome"], f"{img['nome']}_{suf}")
                if os.path.exists(os.path.join(saida, img["nome"], "info.json")):
                    continue
                cmd = [sys.executable, os.path.join(RAIZ, "src", "mama3d.py"), img["arquivo"], "--pixel-mm", str(img["pixel"]),
                       "--voxel", str(img["voxel"]), "--espessura-placa", str(esp), "--camada", "0.2", "--max-placas", "40",
                       "--limiar-fibroso", lim, "--energia-kev", "20", "--material-placa", "PLA", "--saida", saida]
                if etapa[:2] in CERA:
                    cmd += ["--material-fundo", "cera_abelha"]
                r = subprocess.run(cmd, capture_output=True, text=True, env=env, encoding="utf-8", cwd=RAIZ)
                linha = [x for x in r.stdout.splitlines() if x.startswith("[")] or [r.stderr[-200:]]
                print(etapa, img["nome"], suf, linha[-1].encode("ascii", "replace").decode()[:140], flush=True)


def medir(pasta: str) -> dict:
    info = json.load(open(os.path.join(pasta, "info.json"), encoding="utf-8"))
    Z = np.load(os.path.join(pasta, "intermediarios.npz"))
    f = float(info["fator_material"])
    reg = Z["mask"] & (Z["t"] > 0)
    err = np.abs(Z["dq_sem_limpeza"] - Z["d"]) / f
    cad = np.abs(Z["dq"] - Z["d_ref"]) / f
    e = float(info["espessura_placa_modelo_mm"])
    return {"fator": round(f, 4), "espessura_mm": info["espessura_placa_mm"], "n_placas": info["n_placas"],
            "pecas": info["impressao"]["corpos"], "filamento_cm3": info["impressao"]["volume_cm3"],
            "fracao_glandular": info["fracao_glandular_impressa"], "degrau_tecido_mm": round(e / f, 2),
            "quant_medio_mm": round(float(err[reg].mean()), 3),
            "quant_rmse_mm": round(float(np.sqrt((err[reg] ** 2).mean())), 3),
            "quant_p95_mm": round(float(np.percentile(err[reg], 95)), 3),
            "cadeia_medio_mm": round(float(cad[reg].mean()), 3), "_mapa": np.where(reg, err, np.nan)}


def comparar():
    saida = os.path.join(CMP, "quantizacao_todas_etapas")
    os.makedirs(saida, exist_ok=True)
    linhas, mapas = [], {}
    for img in imagens():
        for etapa, rot, _ in ETAPAS:
            for suf, _, lim in LIMIARES:
                p = os.path.join(SAIDA, etapa, img["nome"], f"{img['nome']}_{suf}", img["nome"])
                if not os.path.exists(os.path.join(p, "info.json")):
                    print("faltando:", p)
                    continue
                m = medir(p)
                mapas[(img["nome"], etapa, suf)] = m.pop("_mapa")
                linhas.append({"imagem": img["nome"], "etapa": etapa, "limiar": lim, **m})
    cols = [k for k in linhas[0] if k not in ("imagem", "etapa", "limiar")]
    with open(os.path.join(saida, "quantizacao.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["imagem", "etapa", "limiar", *cols])
        w.writeheader()
        w.writerows(linhas)

    def media(etapa, lim, k):
        v = [r[k] for r in linhas if r["etapa"] == etapa and r["limiar"] == lim]
        return float(np.mean(v)) if v else float("nan")

    with open(os.path.join(saida, "quantizacao.md"), "w", encoding="utf-8") as f:
        f.write("# Erro de quantização por etapa e limiar (média das imagens)\n\n"
                "Erros em **mm de tecido equivalente** (mm de placa ÷ fator), na região da mama. "
                "\"Depois\" = fator por atenuação linear (PLA, 20 keV).\n\n"
                "| limiar | etapa | placas | peças | filamento (cm³) | degrau (mm tecido) | erro médio | RMSE | p95 |\n"
                "|---|---|---|---|---|---|---|---|---|\n")
        for _, _, lim in LIMIARES:
            for etapa, rot, _ in ETAPAS:
                f.write(f"| {lim} | {rot} | {media(etapa, lim, 'n_placas'):.1f} | {media(etapa, lim, 'pecas'):.0f} | "
                        f"{media(etapa, lim, 'filamento_cm3'):.1f} | {media(etapa, lim, 'degrau_tecido_mm'):.2f} | "
                        f"{media(etapa, lim, 'quant_medio_mm'):.3f} | {media(etapa, lim, 'quant_rmse_mm'):.3f} | "
                        f"{media(etapa, lim, 'quant_p95_mm'):.3f} |\n")
        f.write("\n## Por imagem (erro médio, mm de tecido equivalente)\n\n| imagem | limiar | "
                + " | ".join(r for _, r, _ in ETAPAS) + " |\n|---|---|" + "---|" * len(ETAPAS) + "\n")
        for img in imagens():
            for _, _, lim in LIMIARES:
                v = [next((r["quant_medio_mm"] for r in linhas
                           if r["imagem"] == img["nome"] and r["etapa"] == e and r["limiar"] == lim), float("nan"))
                     for e, _, _ in ETAPAS]
                f.write(f"| {img['nome']} | {lim} | " + " | ".join(f"{x:.3f}" for x in v) + " |\n")

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    x = np.arange(len(LIMIARES))
    for j, (etapa, rot, _) in enumerate(ETAPAS):
        ax[0].bar(x + (j - (len(ETAPAS) - 1) / 2) * 0.8 / len(ETAPAS), [media(etapa, lim, "quant_medio_mm") for _, _, lim in LIMIARES], 0.8 / len(ETAPAS), label=rot)
        ax[1].bar(x + (j - (len(ETAPAS) - 1) / 2) * 0.8 / len(ETAPAS), [media(etapa, lim, "n_placas") for _, _, lim in LIMIARES], 0.8 / len(ETAPAS), label=rot)
    for a, tit, yl in ((ax[0], "erro médio de quantização (menor é melhor)", "mm de tecido equivalente"),
                       (ax[1], "número de placas", "placas")):
        a.set_xticks(x, [f"limiar {lim}" for _, _, lim in LIMIARES])
        a.set_title(tit)
        a.set_ylabel(yl)
        a.grid(axis="y", alpha=.3)
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "erro_por_etapa_e_limiar.png"), dpi=120)
    plt.close(fig)

    for img in imagens():  # mapas de erro: linhas = limiar, colunas = etapa (mesma escala de cor)
        vals = [m[~np.isnan(m)] for (n, _, _), m in mapas.items() if n == img["nome"]]
        vmax = max(np.percentile(v, 99) for v in vals)
        fig, ax = plt.subplots(len(LIMIARES), len(ETAPAS), figsize=(4.2 * len(ETAPAS), 3.6 * len(LIMIARES)))
        for i, (suf, _, lim) in enumerate(LIMIARES):
            for j, (etapa, rot, _) in enumerate(ETAPAS):
                m = mapas.get((img["nome"], etapa, suf))
                a = ax[i, j]
                a.axis("off")
                if m is None:
                    continue
                a.imshow(m, cmap="magma", vmin=0, vmax=vmax)
                r = next(r for r in linhas if r["imagem"] == img["nome"] and r["etapa"] == etapa and r["limiar"] == lim)
                a.set_title(f"limiar {lim} · {rot}\n{r['n_placas']} placas · erro médio {r['quant_medio_mm']} mm", fontsize=8)
        fig.suptitle(f"{img['nome']}: erro de quantização (mm de tecido equivalente, mesma escala)")
        fig.tight_layout()
        fig.savefig(os.path.join(saida, f"mapas_{img['nome']}.png"), dpi=90)
        plt.close(fig)

    for etapa, _, esp in ETAPAS[2:]:  # mostradores das novas etapas
        subprocess.run([sys.executable, os.path.join(RAIZ, "src", "mostrador_pf.py"), os.path.join(SAIDA, etapa)], cwd=RAIZ,
                       env=dict(os.environ, PYTHONUTF8="1"), capture_output=True)
        src, dst = os.path.join(SAIDA, etapa, "mostrador_pf"), os.path.join(CMP, f"mostrador_{etapa[:2]}_{etapa[3:]}")
        if os.path.isdir(src):
            shutil.rmtree(dst, ignore_errors=True)
            shutil.move(src, dst)
    print("em:", saida)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("gerar", "comparar"):
        sys.exit(__doc__)
    gerar() if sys.argv[1] == "gerar" else comparar()
