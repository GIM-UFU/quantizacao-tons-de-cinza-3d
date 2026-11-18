"""Recipiente (parte externa da mama) para a peça única vazada de src/peca_unica.py.

Uso:
    python src/recipiente_peca_unica.py saida/07_cera_abelha_1.0mm/lcc1/lcc1_pf06/lcc1
    python src/recipiente_peca_unica.py <pasta com intermediarios.npz> --parede 2 --fundo 1.5 --tampa 1.5

Gera, na mesma pasta da peça única (saida/09_peca_unica/<imagem>_pfXX/):
  recipiente.stl         bandeja com o contorno da mama: fundo + parede (≈ pele), aberta em cima
  tampa.stl              placa com o mesmo contorno (≈ pele superior), com aba de encaixe;
                         imprimir com a aba para cima (sem suporte)
  conjunto_integrado.stl recipiente + peça única num só corpo (a peça nasce do fundo): uma
                         impressão, registro perfeito, a peça não desliza ao verter a cera
  *_yup.stl              cópias para o documento Y-up do Fusion (altura em Y, base em y = 0)

Mesmo sistema de coordenadas da peça única (x = coluna·voxel, y = −linha·voxel): a peça
importada junto cai exatamente dentro da cavidade.

Espessuras padrão: pele da mama ≈ 1,45 mm em média (Huang et al., Med. Phys. 35:1199, 2008,
breast CT), por isso fundo e tampa de 1,5 mm. A parede lateral é atravessada de lado pelo feixe
(até H mm de PLA na borda), então fica tão fina quanto a impressão permite (2 mm = 5 perímetros
de bico 0,4); a borda do phantom vai aparecer mais clara que a linha de pele real.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import trimesh
from scipy import ndimage as ndi
from skimage import morphology

sys.path.insert(0, os.path.dirname(__file__))
from mama3d import reduzir, superficie  # noqa: E402
from peca_unica import DENS_PLA, niveis_relevo, pontes, reamostrar_media  # noqa: E402


def disco(mm: float, v: float):
    return morphology.disk(max(1, int(round(mm / v))))


def malha(vol: np.ndarray, v: float, dz: float, pad: int) -> trimesh.Trimesh:
    m = superficie(vol, v, dz, 0.0, suavizar=(0.7, 0.7, 0))
    m.apply_translation([-pad * v, pad * v, 0])  # volta ao sistema da máscara sem padding
    m.apply_translation([0, 0, -m.bounds[0][2]])
    return m


def yup(m: trimesh.Trimesh, centro_xz: np.ndarray) -> trimesh.Trimesh:
    """Rotaciona para Y-up e centra pelo centro do recipiente (todas as peças com o mesmo deslocamento)."""
    c = m.copy()
    c.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    c.apply_translation([-centro_xz[0], 0, -centro_xz[1]])
    return c


def processar(pasta_in: str, a) -> dict:
    z = np.load(os.path.join(pasta_in, "intermediarios.npz"))
    v0 = float(z["voxel"])
    f = max(1, int(round(a.voxel_malha / v0)))
    v = v0 * f
    d = reamostrar_media(z["d"], f)
    mask = reamostrar_media(z["mask"].astype(np.float32), f) > 0.5
    d[~mask] = 0
    H = float(a.altura) if a.altura else float(z["t"].max())

    # peça única, idêntica à de peca_unica.py (mesmos parâmetros padrão)
    h, _ = niveis_relevo(d, mask, v, a.camada, a.traco_min, a.area_min)
    ponte, _ = pontes(h > 0, mask, v, a.largura_ponte)
    h = np.where(ponte, np.maximum(h, np.round(a.altura_ponte / a.camada) * a.camada), h)

    pad = int(np.ceil((a.parede + a.folga) / v)) + 3
    M = np.pad(mask, pad)
    hP = np.pad(h, pad)
    interno = ndi.binary_dilation(M, disco(a.folga, v)) if a.folga > 0 else M
    externo = ndi.binary_dilation(interno, disco(a.parede, v))

    dz = a.camada
    fundo = np.round(a.fundo / dz) * dz
    nz = int(round((fundo + H) / dz))
    zc = ((np.arange(nz) + 0.5) * dz)[None, None, :]
    rec = externo[..., None] & ((zc < fundo) | ~interno[..., None])
    peca = (zc >= fundo) & (zc < fundo + hP[..., None])

    pasta = os.path.join(a.saida, os.path.basename(os.path.dirname(os.path.normpath(pasta_in))))
    os.makedirs(pasta, exist_ok=True)
    m_rec = malha(rec, v, dz, pad)
    if a.max_faces:  # paredes lisas: a decimação quase não muda o volume; só aceita se continuar estanque
        q = reduzir(m_rec.copy(), a.max_faces)
        m_rec = q if q.is_watertight else m_rec
    m_int = malha(rec | peca, v, dz, pad)
    del rec, peca

    # tampa: placa no contorno externo + aba de encaixe (anel 1,2 mm, 2 mm de altura) por dentro da cavidade
    et = np.round(a.tampa / dz) * dz
    aba_ext = ndi.binary_erosion(interno, disco(0.3, v))          # 0,3 mm de folga na parede
    aba = aba_ext & ~ndi.binary_erosion(aba_ext, disco(1.2, v))
    nzt = int(round((et + 2.0) / dz))
    zt = ((np.arange(nzt) + 0.5) * dz)[None, None, :]
    tam = (externo[..., None] & (zt < et)) | (aba[..., None] & (zt >= et))
    m_tam = malha(tam, v, dz, pad)

    centro = m_rec.bounds.mean(axis=0)
    centro_xz = np.array([centro[0], -centro[1]])  # após girar −90° em X, y vira −z
    saidas = {"recipiente": m_rec, "tampa": m_tam, "conjunto_integrado": m_int}
    for nome, m in saidas.items():
        m.export(os.path.join(pasta, f"{nome}.stl"))
        yup(m, centro_xz).export(os.path.join(pasta, f"{nome}_yup.stl"))

    # verificações
    folga_min = float(ndi.distance_transform_edt(interno)[hP > 0].min() * v)  # peça → parede
    area_cav = float(interno.sum()) * v**2
    vol_peca = float(h.sum()) * v**2  # mm³ (colunas exatas, sem a suavização da malha)
    mesa = np.asarray(a.mesa)
    info = {
        "origem": pasta_in.replace("\\", "/"), "voxel_malha_mm": v, "camada_mm": dz,
        "altura_cavidade_mm": H, "parede_mm": a.parede, "fundo_mm": float(fundo), "tampa_mm": float(et),
        "folga_mm": a.folga, "folga_minima_peca_parede_mm": round(folga_min, 2),
        "cavidade_cm3": round(area_cav * H / 1000, 1),
        "cera_mL": round((area_cav * H - vol_peca) / 1000, 1),
        "altura_peca_mm": round(float(h.max()), 2),
        "cera_acima_da_peca_mm": round(H - float(h.max()), 2),
    }
    for nome, m in saidas.items():
        info[nome] = {"faces": int(len(m.faces)), "corpos": int(m.body_count), "volume_valido": bool(m.is_volume),
                      "dimensoes_mm": [round(float(x), 1) for x in m.extents],
                      "PLA_cm3": round(float(m.volume) / 1000, 1), "PLA_g": round(float(m.volume) / 1000 * DENS_PLA, 1),
                      "cabe_na_mesa": bool(np.all(m.extents <= mesa))}
    info["avisos"] = [
        "a parede é atravessada de lado pelo feixe: a borda do phantom aparece mais clara que a pele real",
        "no recipiente separado a peça pode deslizar ao verter a cera; o conjunto integrado evita isso",
        "cera de abelha ≈ 62–64 °C × PLA Tg ≈ 55–60 °C: validar com o corpo de prova (saida/10_cubo_teste)",
    ]
    with open(os.path.join(pasta, "info_recipiente.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, ensure_ascii=False)
    info["pasta"] = pasta
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entradas", nargs="+", help="pastas com intermediarios.npz (as mesmas da peça única)")
    ap.add_argument("--parede", type=float, default=2.0, help="parede lateral em mm (padrão 2)")
    ap.add_argument("--fundo", type=float, default=1.5, help="fundo em mm (≈ pele, padrão 1.5)")
    ap.add_argument("--tampa", type=float, default=1.5, help="tampa em mm (≈ pele, padrão 1.5)")
    ap.add_argument("--folga", type=float, default=0.5, help="folga entre contorno da mama e parede (padrão 0.5)")
    ap.add_argument("--altura", type=float, default=None, help="altura da cavidade em mm (padrão: espessura H)")
    ap.add_argument("--camada", type=float, default=0.2)
    ap.add_argument("--voxel-malha", dest="voxel_malha", type=float, default=0.5)
    ap.add_argument("--traco-min", dest="traco_min", type=float, default=1.2)
    ap.add_argument("--area-min", dest="area_min", type=float, default=25.0)
    ap.add_argument("--largura-ponte", dest="largura_ponte", type=float, default=1.5)
    ap.add_argument("--altura-ponte", dest="altura_ponte", type=float, default=0.6)
    ap.add_argument("--max-faces", dest="max_faces", type=int, default=200000,
                    help="decimação do recipiente (padrão 200000; 0 desliga)")
    ap.add_argument("--mesa", type=float, nargs=3, default=[220, 220, 250], metavar=("X", "Y", "Z"))
    ap.add_argument("--saida", default=os.path.join("saida", "09_peca_unica"))
    a = ap.parse_args(argv)
    for p in a.entradas:
        info = processar(p, a)
        print(json.dumps({k: info[k] for k in info if k != "avisos"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
