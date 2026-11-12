"""Mamografia (DICOM ou TIFF) → modelo 3D (STL, GLB colorido e visualizador HTML).

Uso:
    python src/mama3d.py dados/ --modo solido   # padrão
    python src/mama3d.py dados/ --modo relevo
    python src/mama3d.py dados/ --modo volume   # apenas para DBT (multi-frame)

Modos
-----
solido  (modelo de forma, recomendado para 2D)
    A mama comprimida fica entre duas placas paralelas separadas por H
    (BodyPartThickness). No centro a espessura é H; na borda da pele ela cai com
    perfil semicircular de raio H/2. Isso dá um sólido fechado (estanque, imprimível)
    com escala em mm. A intensidade da mamografia vira a COR da superfície.
    → Forma = hipótese geométrica plausível; NÃO é reconstrução anatômica.

relevo  (mapa de altura)
    z = intensidade normalizada × escala. Visualização da distribuição de tecido
    denso como "terreno". Sem significado físico no eixo z.

volume  (DBT / tomossíntese)
    Dado já é 3D: isosuperfície por marching cubes (Lorensen & Cline, 1987).
    Único modo que representa anatomia real em z.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, transform

sys.path.insert(0, os.path.dirname(__file__))
from dicom_io import Mamografia, carregar, listar_dicoms  # noqa: E402


# ───────────────────────── pré-processamento 2D ─────────────────────────

def reamostrar(img: np.ndarray, pixel_mm, voxel_mm: float) -> np.ndarray:
    """Leva a imagem para pixels isotrópicos de `voxel_mm`.
    Mamografias têm ~0.05–0.1 mm/pixel (≈ 3000×4000); a 1 mm a malha fica leve."""
    escala = (pixel_mm[0] / voxel_mm, pixel_mm[1] / voxel_mm)
    return transform.rescale(img, escala, anti_aliasing=True, preserve_range=True).astype(np.float32)


def linearizar(img: np.ndarray) -> tuple[np.ndarray, bool]:
    """Garante 'tecido = claro' e, em imagens de transmissão (RAW), converte para atenuação.

    Em imagem RAW/"FOR PROCESSING" o pixel é ~ proporcional à intensidade transmitida
    I = I0·exp(−∫μ dz): o AR é o mais claro e o tecido denso é escuro. Em
    "FOR PRESENTATION" é o contrário. Detectamos pelo fundo: a borda da imagem é
    quase toda ar (a parede torácica ocupa no máximo um dos quatro lados). Se a
    maioria da borda cai na classe clara do Otsu, a imagem é de transmissão e
    aplicamos A = ln(I0/I) — a mesma linearização de Beer–Lambert do mammoreplicator
    (Badal et al., J. Med. Imaging 5:033501, 2018). A resultante é ∝ ∫μ dz."""
    k = max(2, min(img.shape) // 50)
    borda = np.concatenate([img[:k].ravel(), img[-k:].ravel(), img[:, :k].ravel(), img[:, -k:].ravel()])
    limiar = filters.threshold_otsu(img)
    if np.mean(borda > limiar) <= 0.5:
        return img, False
    I0 = float(np.median(borda[borda > limiar]))  # intensidade do ar ≈ feixe incidente
    return np.log(I0 / np.clip(img, 1.0, None)).astype(np.float32), True


def segmentar_mama(img: np.ndarray, fator_limiar: float = 1.0) -> np.ndarray:
    """Máscara da mama: limiar + maior componente conexa + preenchimento de buracos.
    A maior componente descarta marcadores/etiquetas no fundo.

    Limiar: em mamografia real o histograma tem um pico enorme no fundo (ar) e uma
    cauda longa (gordura escura → glândula clara). O Otsu, que busca duas classes
    equilibradas, cai ENTRE gordura e glândula e corta a gordura fora da máscara.
    O método do triângulo (Zack et al., 1977) é feito para esse caso: acha o 'pé'
    do pico do fundo. Usamos o menor dos dois; se o resultado cobrir quase a imagem
    toda (fundo ruidoso), voltamos ao Otsu."""
    otsu = filters.threshold_otsu(img)
    limiar = min(otsu, filters.threshold_triangle(img)) * fator_limiar
    m = img > limiar
    m = ndi.binary_opening(m, morphology.disk(2))
    if m.mean() > 0.9:
        m = ndi.binary_opening(img > otsu * fator_limiar, morphology.disk(2))
    lab = measure.label(m)
    if lab.max() == 0:
        raise RuntimeError("Segmentação vazia; ajuste --limiar.")
    maior = np.argmax(np.bincount(lab.ravel())[1:]) + 1
    m = ndi.binary_fill_holes(lab == maior)
    return ndi.binary_closing(m, morphology.disk(3)) | m


def orientar(img: np.ndarray, mask: np.ndarray):
    """Padroniza: parede torácica à esquerda (x = 0). Detecta pelo lado com mais mama."""
    w = mask.shape[1]
    faixa = max(1, w // 10)
    if mask[:, :faixa].sum() < mask[:, -faixa:].sum():
        return img[:, ::-1].copy(), mask[:, ::-1].copy(), True
    return img, mask, False


def remover_peitoral(img: np.ndarray, mask: np.ndarray, voxel_mm: float, margem_mm: float = 2.0):
    """MLO: devolve a máscara do músculo peitoral (para excluir da mama).

    Na MLO o peitoral é um triângulo muito claro no canto superior, junto à parede
    torácica (imagem já orientada, parede à esquerda). Sem removê-lo, o gerador
    imprimiria o músculo como se fosse tecido fibroglandular.
    Método: (1) região clara (> percentil 80 da mama, suavizada) ligada ao canto
    superior esquerdo; (2) sua borda direita, linha a linha, vira uma reta por
    mínimos quadrados aparados (a fronteira do peitoral é quase reta); (3) tudo à esquerda
    da reta, mais `margem_mm`, é peitoral. Heurístico: conferir o resultado na imagem."""
    h, w = img.shape
    s = ndi.gaussian_filter(img, 2.0 / voxel_mm)
    alto = (s > np.percentile(s[mask], 80)) & mask
    lab, _ = ndi.label(alto)
    semente = lab[: max(1, int(0.08 * h)), : max(1, int(0.25 * w))]
    if not (semente > 0).any():
        return np.zeros_like(mask)
    alvo = lab == np.bincount(semente[semente > 0]).argmax()
    ys = np.where(alvo.any(1))[0]
    if len(ys) < 10:
        return np.zeros_like(mask)
    xs = np.array([np.where(alvo[y])[0].max() for y in ys], dtype=float)
    keep = np.ones(len(ys), bool)
    for _ in range(3):  # mínimos quadrados aparados: descarta os 20 % piores
        a, b = np.polyfit(ys[keep], xs[keep], 1)
        res = np.abs(xs - (a * ys + b))
        keep = res <= np.percentile(res[keep], 80)
    Y, X = np.mgrid[0:h, 0:w]
    reta = a * Y + b + margem_mm / voxel_mm
    peit = (X <= reta) & mask & (a < 0)  # inclinação negativa = triângulo no canto superior
    # mantém só o que está ligado ao canto (descarta pontas soltas)
    lab2, _ = ndi.label(peit)
    c = lab2[: max(1, int(0.08 * h)), : max(1, int(0.25 * w))]
    return lab2 == np.bincount(c[c > 0]).argmax() if (c > 0).any() else np.zeros_like(mask)


def cortar_abdome(mask: np.ndarray, voxel_mm: float):
    """MLO: corta a faixa de tronco/abdome que a máscara inclui abaixo do sulco inframamário.

    Com a parede torácica à esquerda, o contorno inferior da mama sobe a partir do sulco;
    a faixa de tronco continua colada à borda esquerda até o fim da imagem. Acha a primeira
    coluna (da parede para fora) em que a máscara deixa de chegar à base da imagem e corta,
    nas colunas anteriores, tudo abaixo do contorno inferior dessa coluna.
    Devolve a máscara do que foi cortado."""
    h, w = mask.shape
    base = np.array([np.where(mask[:, x])[0].max() if mask[:, x].any() else -1 for x in range(w)])
    if base[0] < 0:
        return np.zeros_like(mask)
    # contorno inferior da mama de referência: mediana do fundo da máscara entre 15 e 30 mm da
    # parede (já fora da faixa de tronco). Colunas mais perto da parede que ainda descem mais
    # de 10 mm além dessa referência são tronco: corta abaixo da referência.
    x1, x2 = int(round(15 / voxel_mm)), int(round(30 / voxel_mm))
    if x2 >= w:
        return np.zeros_like(mask)
    ref = int(np.median(base[x1:x2]))
    excesso = np.where(base[:x1] > ref + int(round(10 / voxel_mm)))[0]
    if len(excesso) == 0:
        return np.zeros_like(mask)
    xf = int(excesso.max()) + 1
    corte = np.zeros_like(mask)
    corte[ref + int(round(2 / voxel_mm)):, :xf + int(round(2 / voxel_mm))] = True
    return corte & mask


def intensidade_normalizada(img: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Proxy de densidade em [0,1] dentro da mama (percentis 1–99).
    Não é densidade calibrada: para isso seria preciso imagem FOR PROCESSING e o
    modelo físico do equipamento (cf. Highnam & Brady, 1999)."""
    v = img[mask]
    lo, hi = np.percentile(v, [1, 99])
    p = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)
    p[~mask] = 0
    return p


def mapa_espessura(mask: np.ndarray, H_mm: float, voxel_mm: float) -> np.ndarray:
    """Espessura local t(x,y) em mm.
    d = distância até a pele (a borda da imagem do lado da parede torácica NÃO conta
    como pele, pois a EDT do scipy só mede distância até pixels zero do array).
    Perfil semicircular de raio R = H/2:  t = 2·sqrt(R² − (R − min(d,R))²)."""
    d = ndi.distance_transform_edt(mask) * voxel_mm
    R = H_mm / 2
    dc = np.minimum(d, R)
    t = 2 * np.sqrt(np.maximum(R**2 - (R - dc) ** 2, 0))
    t[~mask] = 0
    return t


# ───────────────────────────── malhas ─────────────────────────────

def colorir(malha: trimesh.Trimesh, campo: np.ndarray, voxel_mm: float, cmap="bone"):
    """Cor por vértice amostrando o campo 2D na posição (x, y) do vértice."""
    h, w = campo.shape
    xs = np.clip((malha.vertices[:, 0] / voxel_mm).round().astype(int), 0, w - 1)
    ys = np.clip((-malha.vertices[:, 1] / voxel_mm).round().astype(int), 0, h - 1)
    rgba = (plt.get_cmap(cmap)(campo[ys, xs]) * 255).astype(np.uint8)
    malha.visual.vertex_colors = rgba
    return malha


def solido_entre(zmin: np.ndarray, zmax: np.ndarray, mask: np.ndarray, voxel_mm: float,
                 dz_mm: float | None = None, suavizar: float = 1.0) -> trimesh.Trimesh:
    """Voxeliza {(x,y,z) : mask(x,y) e zmin(x,y) ≤ z ≤ zmax(x,y)} e extrai a superfície
    por marching cubes. Voxelizar (em vez de costurar topo+fundo+laterais à mão)
    garante, por construção, uma malha fechada (estanque) e 2-manifold — o que o
    fatiador (slicer) precisa para saber o que é 'dentro' e 'fora'."""
    dz = dz_mm or voxel_mm
    z0, z1 = zmin[mask].min(), zmax[mask].max()
    nz = int(np.ceil((z1 - z0) / dz)) + 1
    z = z0 + np.arange(nz) * dz
    vol = (mask[:, :, None] & (z >= zmin[:, :, None] - 1e-6) & (z <= zmax[:, :, None] + 1e-6))
    return superficie(vol, voxel_mm, dz, z0, suavizar)


def superficie(vol: np.ndarray, voxel_mm: float, dz: float, z0: float,
               suavizar: float = 1.0) -> trimesh.Trimesh | None:
    """Volume binário (linha, coluna, z) → malha fechada por marching cubes.

    Propriedade usada no modo 'composto': se A e B são rótulos complementares
    dentro de um sólido, G*A + G*B = G*(A∪B) (filtro gaussiano é linear). Onde a
    união vale 1, a isosuperfície 0.5 de A coincide com a de B → as partes
    compartilham exatamente a mesma interface, sem frestas nem sobreposição."""
    if not vol.any():
        return None
    vol = np.pad(vol.astype(np.float32), 2)  # borda vazia → superfície fechada
    if isinstance(suavizar, tuple):   # (σ_linha, σ_coluna, σ_z) em voxels — p.ex. σ_z = 0
        vol = ndi.gaussian_filter(vol, suavizar)  # preserva a espessura exata de placas finas
    elif suavizar > 0:
        vol = ndi.gaussian_filter(vol, (suavizar, suavizar, suavizar * voxel_mm / dz))
    if vol.max() < 0.5:
        return None  # estrutura menor que a resolução: some na suavização
    verts, faces, _, _ = measure.marching_cubes(vol, level=0.5, spacing=(voxel_mm, voxel_mm, dz))
    verts -= (2 * voxel_mm, 2 * voxel_mm, 2 * dz)  # desfaz o pad
    # (linha, coluna, z) → (x, y, z) com y para cima
    v = np.c_[verts[:, 1], -verts[:, 0], verts[:, 2] + z0]
    m = trimesh.Trimesh(v, faces[:, ::-1], process=True)
    trimesh.repair.fix_normals(m)
    return limpar_degeneradas(m)


def limpar_degeneradas(m: trimesh.Trimesh) -> trimesh.Trimesh:
    """Marching cubes ocasionalmente gera triângulos de área ~0 (vértices colineares).
    Removemos e, se a malha abrir, fundimos os vértices próximos; só aceitamos o
    resultado se continuar estanque — senão mantém a original."""
    deg = m.area_faces < 1e-10
    if not deg.any():
        return m
    c = m.copy()
    c.update_faces(~deg)
    c.remove_unreferenced_vertices()
    if not c.is_watertight:
        c.merge_vertices(digits_vertex=4)
        c.update_faces(c.nondegenerate_faces())
        c.remove_unreferenced_vertices()
    if c.is_watertight and c.is_winding_consistent:
        trimesh.repair.fix_normals(c)
        return c
    return m


# ─────────────────── modo composto: interior por tons de cinza ───────────────────

def volume_rotulado(t: np.ndarray, p: np.ndarray, mask: np.ndarray, niveis: int,
                    relevo_mm: float, dz: float):
    """Volume de rótulos: 0 = fora, 1 = adiposo, 2..niveis = tecido denso por faixa de cinza.

    Modelo de duas componentes (adiposo + fibroglandular), padrão na análise de
    densidade mamográfica: cada pixel é a projeção de uma coluna de altura t(x,y)
    com uma fração p(x,y) de tecido denso. Colocamos esse tecido denso como uma
    camada centrada em z = 0 com espessura p·t:
        denso(x,y,z)  ⇔  |z| ≤ p(x,y)·t(x,y)/2
    → a QUANTIDADE de tecido denso em cada coluna reproduz o tom de cinza (o que a
      imagem de fato mede); a POSIÇÃO em z (centrada) é hipótese — a projeção não
      informa a profundidade.
    Superfície superior: z_topo = t/2 + relevo_mm·p  → os tons também viram relevo externo.
    Com niveis > 2 a camada densa é subdividida por faixas de p (cinzas mais claros
    → rótulos mais altos), p.ex. para infills/materiais diferentes por faixa."""
    zb = -t / 2
    zt = t / 2 + relevo_mm * p
    z0, z1 = zb[mask].min(), zt[mask].max()
    nz = int(np.ceil((z1 - z0) / dz)) + 1
    z = (z0 + np.arange(nz) * dz)[None, None, :]

    dentro = mask[..., None] & (z >= zb[..., None]) & (z <= zt[..., None])
    denso = dentro & (np.abs(z) <= (p * t / 2)[..., None])
    faixa = 2 + np.minimum((p * (niveis - 1)).astype(np.int16), niveis - 2)
    rot = np.where(denso, faixa[..., None], dentro).astype(np.uint8)
    return rot, z0


def nomes_rotulos(niveis: int) -> dict[int, str]:
    return {1: "adiposo"} | {k: (f"denso_{k - 1}" if niveis > 2 else "denso")
                             for k in range(2, niveis + 1)}


def partes_compostas(rot: np.ndarray, niveis: int, voxel_mm: float, dz: float, z0: float):
    partes = {}
    for k, nome in nomes_rotulos(niveis).items():
        m = superficie(rot == k, voxel_mm, dz, z0)
        if m is not None:
            partes[nome] = m
    return partes


def cor_parte(nome: str, niveis: int) -> np.ndarray:
    if nome == "adiposo":
        return np.array([235, 205, 170, 255], np.uint8)
    k = int(nome.split("_")[1]) if "_" in nome else 1
    f = (k - 1) / max(niveis - 2, 1)
    return (np.array([215, 95, 95, 255]) * (1 - f) + np.array([110, 15, 45, 255]) * f).astype(np.uint8)


def malha_solida(t: np.ndarray, mask: np.ndarray, voxel_mm: float) -> trimesh.Trimesh:
    """Mama comprimida: sólido simétrico {|z| ≤ t(x,y)/2}."""
    return solido_entre(-t / 2, t / 2, mask, voxel_mm)


def malha_relevo(campo: np.ndarray, mask: np.ndarray, voxel_mm: float,
                 escala_mm: float, base_mm: float) -> trimesh.Trimesh:
    """Relevo com base plana: 0 ≤ z ≤ base + campo·escala. Sólido fechado, imprimível
    deitado (base na mesa). dz menor que o voxel preserva a resolução vertical."""
    zmin = np.zeros_like(campo)
    zmax = base_mm + campo * escala_mm
    return solido_entre(zmin, zmax, mask, voxel_mm, dz_mm=min(voxel_mm, 0.5))


# ───────────── modo placas: só o tecido fibroglandular, em placas empilháveis ─────────────
#
# Segue a montagem de Bliznakova et al., Front. Phys. 12:1441740 (2024): imprime-se
# apenas a estrutura densa (PLA/resina, 100 % de infill) e o tecido adiposo é
# substituído por um líquido (parafina líquida, ρ = 0.82 g/cm³) dentro de um
# recipiente impresso com o formato da mama comprimida, cuja parede faz o papel da pele.
#
# Por que placas: numa mamografia 2D o que importa é a espessura TOTAL de material
# denso atravessada por cada raio (Beer–Lambert: I = I0·exp(−∫μ dz)); a posição em z
# não altera a projeção. Logo o mapa de espessura d(x,y) pode ser fatiado em N placas
# planas de espessura e — a placa k cobre a região onde d ≥ (k−½)·e. Empilhadas, somam
# d quantizado em passos de e. Placas planas imprimem deitadas, sem suporte.

def mapa_fibroso(img: np.ndarray, mask: np.ndarray, voxel_mm: float, margem_pele_mm: float,
                 fator_zero: float = 1.0, suavizacao_mm: float = 0.75, margem_borda_mm: float = 3.0):
    """Fração fibroglandular pf(x,y) ∈ [0,1]; gordura → 0 (não é impressa).

    Dentro da mama, o Otsu separa as duas populações (adiposo × fibroglandular); esse
    limiar vira o zero da escala e o percentil 99.5 vira 1. A faixa junto à pele é
    zerada: em imagens 'FOR PRESENTATION' a borda é realçada pelo processamento
    (equalização periférica) e a pele já é representada pela parede do recipiente."""
    v = img[mask]
    um = float(np.percentile(v, 99.5))
    piso = float(np.percentile(v, 1))
    # fator_zero < 1 baixa o limiar (inclui feixes fibrosos mais tênues); > 1 só o mais denso
    zero = piso + (float(filters.threshold_otsu(v)) - piso) * fator_zero
    pf = np.clip((img - zero) / (um - zero + 1e-9), 0, 1)
    pf[~mask] = 0
    # Borda da imagem: o filtro/borda do detector deixa uma faixa clara de ~2–3 mm que não
    # é tecido (cinza 98 contra 32 ao lado, medido no lcc1). Ela virava lascas de placa
    # coladas ao limite da imagem, fora do contorno real da mama.
    kb = int(round(margem_borda_mm / voxel_mm))
    borda = np.zeros_like(mask, dtype=bool)
    if kb > 0:
        borda[:, :kb] = borda[:, -kb:] = borda[:kb, :] = borda[-kb:, :] = True
    pf[borda] = 0
    pf_bruto = pf.copy()  # conversão cinza → fração, antes de faixa da pele e suavização
    if margem_pele_mm > 0:
        pf[ndi.distance_transform_edt(mask) * voxel_mm < margem_pele_mm] = 0
    if suavizacao_mm > 0:
        pf = ndi.gaussian_filter(pf, suavizacao_mm / voxel_mm)  # tira ruído sem apagar os feixes
    pf[~mask] = 0
    pf[borda] = 0  # a suavização não pode devolver tecido à faixa da borda
    return pf, zero, um, piso, pf_bruto


def limpar_regiao(R: np.ndarray, voxel_mm: float, traco_min_mm: float, area_min_mm2: float):
    """Remove o que a impressora não reproduz ou não dá para manusear: traços mais
    finos que `traco_min_mm`, ilhas menores que `area_min_mm2` e furos < 4 mm²."""
    # disk(r) tem largura 2r+1 pixels; a abertura remove traços mais estreitos que isso.
    # (antes: r = round(traco/2/voxel), que superestimava a largura removida em voxels grossos)
    r = int(round((traco_min_mm / voxel_mm - 1) / 2))
    if r >= 1:
        R = ndi.binary_opening(R, morphology.disk(r))
    lab, n = ndi.label(R)
    if n:
        area = np.bincount(lab.ravel()) * voxel_mm**2
        area[0] = 0
        R = (area >= area_min_mm2)[lab]
    furos, n = ndi.label(~R)
    if n:
        area = np.bincount(furos.ravel()) * voxel_mm**2
        area[0] = np.inf
        R = R | (area < 4.0)[furos]
    return R


def malha_placa(R: np.ndarray, mask: np.ndarray, voxel_mm: float, e: float, base: float):
    """Placa plana de espessura e sobre a região R. Com base > 0, uma lâmina de
    espessura `base` cobre toda a mama e une as ilhas numa peça só (ao custo de uma
    camada uniforme extra de material)."""
    dz = e / max(4, int(round(e / 0.25)))  # dz divide e exatamente → espessura sem arredondar
    base = round(base / dz) * dz
    reg = mask if base > 0 else R
    topo = np.where(R, base + e, base).astype(np.float32)
    # amostras em z = 0, dz, …: (topo − dz) faz a peça ter exatamente a espessura pedida
    m = solido_entre(np.zeros_like(topo), topo - dz, reg, voxel_mm, dz_mm=dz, suavizar=(0.7, 0.7, 0))
    if m is not None:
        m.apply_translation([0, 0, -m.bounds[0][2]])
    return m


def malha_recipiente(mask: np.ndarray, voxel_mm: float, H: float, parede: float, fundo: float,
                     folga: float):
    """Bandeja com o contorno da mama: parede (≈ pele) + fundo, aberta em cima.
    Cavidade = contorno da mama + folga, altura H (espessura de compressão)."""
    pad = int(np.ceil((parede + folga) / voxel_mm)) + 3
    M = np.pad(mask, pad)
    interno = ndi.binary_dilation(M, morphology.disk(max(1, int(round(folga / voxel_mm))))) if folga > 0 else M
    externo = ndi.binary_dilation(interno, morphology.disk(max(2, int(round(parede / voxel_mm)))))
    z = (np.arange(int(np.ceil((fundo + H) / voxel_mm)) + 1) * voxel_mm)[None, None, :]
    vol = externo[..., None] & ((z < fundo) | ~interno[..., None])
    m = superficie(vol, voxel_mm, voxel_mm, 0.0, suavizar=0.6)
    m.apply_translation([-pad * voxel_mm, pad * voxel_mm, -m.bounds[0][2]])
    return m, float(interno.sum() * voxel_mm**2)


def guia_montagem(caminho: str, regioes, mask: np.ndarray, voxel_mm: float, e: float):
    """PDF em escala 1:1 (imprimir a 100 %): uma página por placa, com o contorno da
    mama, para posicionar as peças de cada camada."""
    from matplotlib.backends.backend_pdf import PdfPages
    h, w = mask.shape
    ext = [0, w * voxel_mm, h * voxel_mm, 0]
    with PdfPages(caminho) as pdf:
        for k, R in enumerate(regioes, 1):
            fig = plt.figure(figsize=(w * voxel_mm / 25.4, h * voxel_mm / 25.4))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.imshow(np.where(R, 0.45, 1.0), cmap="gray", vmin=0, vmax=1, extent=ext, interpolation="nearest")
            ax.contour(mask, levels=[0.5], colors="k", linewidths=0.6, extent=[0, w * voxel_mm, 0, h * voxel_mm], origin="upper")
            n = ndi.label(R)[1]
            ax.text(0.98, 0.02, f"placa {k}/{len(regioes)} · {e:g} mm · {n} peça(s) · escala 1:1",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
            ax.plot([w * voxel_mm - 60, w * voxel_mm - 10], [8, 8], "k-", lw=1.5)
            ax.text(w * voxel_mm - 35, 6, "50 mm", ha="center", va="bottom", fontsize=7)
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)


def processar_placas(a, info: dict, pasta: str, nome: str, img, mask, t, H: float) -> dict:
    os.makedirs(os.path.join(pasta, "placas"), exist_ok=True)
    v = a.voxel
    s = a.escala_impressao  # fator linear do protótipo (1 = tamanho real)
    pf, zero, um, piso, pf_bruto = mapa_fibroso(img, mask, v, a.margem_pele, a.limiar_fibroso, a.suavizacao,
                                                a.margem_borda)

    # espessura de material denso por pixel (mm)
    d = pf * t
    fracao_bruta = float(d.sum() / t.sum())
    k = 1.0
    if a.fracao_alvo:
        # reescala para a fração glandular volumétrica desejada (literatura: 13.7–25.6 %
        # em média — Yaffe et al., Med. Phys. 36:5437, 2009, 'The myth of the 50-50 breast')
        for _ in range(20):  # iterativo porque d é limitado a t
            atual = np.minimum(k * d, t).sum() / t.sum()
            if abs(atual - a.fracao_alvo) < 1e-4 or atual == 0:
                break
            k *= a.fracao_alvo / atual
    d = np.minimum(k * d, t) * a.fator_material
    d_ref = np.minimum(k * pf_bruto * t, t) * a.fator_material  # relevo "ideal", sem perdas

    # número de placas
    def _n(esp):
        return int(np.ceil((d.max() - 0.5 * esp) / esp)) if d.max() > 0.5 * esp else 0

    def _camada(esp):  # fatiador: a placa vira nº inteiro de camadas (arredonda p/ cima)
        return float(np.ceil(esp / a.camada - 1e-9) * a.camada) if a.camada > 0 else esp

    e = _camada(a.espessura_placa)
    n = _n(e)
    if n > a.max_placas:
        e = _camada(float(d.max() / (a.max_placas + 0.5)))
        n = _n(e)
        info["aviso_placas"] = f"espessura da placa aumentada para {e:.2f} mm para caber em {n} placas"
    if n == 0:
        raise RuntimeError("Nenhum tecido denso acima de meia placa; reduza --espessura-placa.")

    # restrições de manufatura valem no tamanho IMPRESSO: com escala s, o traço mínimo
    # de traco_min mm impressos corresponde a traco_min/s mm no modelo (área: /s²)
    regioes, anterior = [], mask
    for i in range(1, n + 1):
        R = limpar_regiao(d >= (i - 0.5) * e, v, a.traco_min / s, a.area_min / s**2) & anterior
        if not R.any():
            break
        regioes.append(R)
        anterior = R
    n = len(regioes)
    dq = np.sum(regioes, axis=0) * e  # espessura efetivamente impressa (quantizada)
    dq_sem_limpeza = sum((d >= (i - 0.5) * e) & mask for i in range(1, n + 1)) * e  # só a quantização
    np.savez_compressed(os.path.join(pasta, "intermediarios.npz"), img=img, mask=mask, t=t,
                        pf_bruto=pf_bruto, d_ref=d_ref, d=d, dq_sem_limpeza=dq_sem_limpeza, dq=dq,
                        piso=piso, zero=zero, um=um, voxel=v, e=e, escala=s)

    cmap = plt.get_cmap("magma")
    cena = trimesh.Scene()
    tabela, ok, vol_total, ext_max, pecas_total = [], True, 0.0, np.zeros(3), 0
    fora_n, x_min, mask_d = 0, np.inf, ndi.binary_dilation(mask)  # tolerância de 1 voxel
    H_px, W_px = mask.shape
    for i, R in enumerate(regioes, 1):
        m = reduzir(malha_placa(R, mask, v, e, a.base_placa), a.max_faces)
        # verificação: nenhum vértice de placa pode cair fora do contorno da mama
        cols = np.round(m.vertices[:, 0] / v).astype(int)
        lins = np.round(-m.vertices[:, 1] / v).astype(int)
        dentro_img = (cols >= 0) & (cols < W_px) & (lins >= 0) & (lins < H_px)
        ok_v = dentro_img.copy()
        ok_v[dentro_img] = mask_d[lins[dentro_img], cols[dentro_img]]
        fora_n += int((~ok_v).sum())
        x_min = min(x_min, float(m.vertices[:, 0].min()))
        if s != 1.0:
            m.apply_scale(s)
        m.export(os.path.join(pasta, "placas", f"placa_{i:02d}.stl"))
        valido = bool(m.is_volume)
        ok &= valido
        vol_total += float(m.volume) / 1000
        ext_max = np.maximum(ext_max, m.extents)
        pecas_total += int(m.body_count)
        tabela.append({"placa": i, "pecas": int(m.body_count),
                       "area_cm2": round(float(R.sum() * (v * s) ** 2 / 100), 1),
                       "volume_cm3": round(float(m.volume) / 1000, 2), "faces": int(len(m.faces)),
                       "volume_valido": valido})
        leve = reduzir(m.copy(), 30000)
        leve.apply_translation([0, 0, s * ((i - 1) * (e + a.base_placa) - H / 2)])  # pilha apoiada no fundo
        leve.visual.face_colors = (np.array(cmap(0.25 + 0.65 * (i - 1) / max(n - 1, 1))) * 255).astype(np.uint8)
        cena.add_geometry(leve, geom_name=f"placa_{i:02d}", node_name=f"placa_{i:02d}")

    # recipiente e contorno da mama em resolução ≥ 1 mm (peças grandes e lisas)
    vc = max(v, 1.0)
    mask_c = transform.rescale(mask.astype(np.float32), v / vc, order=1) > 0.5 if vc != v else mask
    # parede/fundo/folga são medidas de impressão: no modelo valem /s, para que, depois
    # de escalar, a parede impressa continue com `a.parede` mm
    rec, area_cav = malha_recipiente(mask_c, vc, H, a.parede / s, a.fundo / s, a.folga / s)
    rec = reduzir(rec, a.max_faces)
    if s != 1.0:
        rec.apply_scale(s)
    rec.export(os.path.join(pasta, "recipiente.stl"))
    cav_cm3 = area_cav * H * s**3 / 1000
    casca = reduzir(malha_solida(mapa_espessura(mask_c, H, vc), mask_c, vc), 30000)
    if s != 1.0:
        casca.apply_scale(s)
    casca.visual.face_colors = cor_parte("adiposo", 2)
    cena.add_geometry(casca, geom_name="adiposo", node_name="adiposo")
    exportar(cena, pasta, nome)
    guia_montagem(os.path.join(pasta, "guia_montagem_1para1.pdf"), regioes, mask, v * s, e * s)

    # painel
    fig, ax = plt.subplots(1, 4, figsize=(15, 4.6))
    lo, hi = np.percentile(img, [1, 99])
    ax[0].imshow(np.clip((img - lo) / (hi - lo + 1e-9), 0, 1), cmap="gray"); ax[0].set_title("original + contorno da mama", fontsize=9)
    ax[0].contour(mask, levels=[0.5], colors="c", linewidths=0.6)
    ax[1].imshow(pf, cmap="bone"); ax[1].set_title("fração fibroglandular pf (gordura = 0)", fontsize=9)
    h1 = ax[2].imshow(d, cmap="viridis"); ax[2].set_title("espessura densa d(x,y) [mm]", fontsize=9)
    fig.colorbar(h1, ax=ax[2], fraction=0.046)
    h2 = ax[3].imshow(dq, cmap="viridis", vmin=0, vmax=max(d.max(), 1e-6)); ax[3].set_title(f"impresso: {n} placas × {e:g} mm", fontsize=9)
    fig.colorbar(h2, ax=ax[3], fraction=0.046)
    for x in ax:
        x.axis("off")
    fig.suptitle(nome, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "painel.png"), dpi=110)
    plt.close(fig)

    cabe = bool(np.all(ext_max <= np.asarray(a.mesa)) and np.all(rec.extents <= np.asarray(a.mesa)))
    altura_pilha = n * (e + a.base_placa) * s
    info |= {
        "escala_impressao": s, "traco_min_impresso_mm": a.traco_min, "area_min_impressa_mm2": a.area_min,
        "margem_borda_mm": a.margem_borda, "vertices_das_placas_fora_da_mama": fora_n,
        "x_min_das_placas_mm": round(x_min, 2),
        "limiar_gordura_glandula": zero, "nivel_100pct": um, "margem_pele_mm": a.margem_pele,
        "fracao_glandular_bruta": round(fracao_bruta, 3),
        "fracao_glandular_final": round(float((d / a.fator_material).sum() / t.sum()), 3),
        "fracao_glandular_impressa": round(float((dq / a.fator_material).sum() / t.sum()), 3),
        "fator_escala_k": round(k, 3), "fator_material": a.fator_material,
        "n_placas": n, "espessura_placa_mm": round(e * s, 3), "base_placa_mm": a.base_placa * s,
        "espessura_placa_modelo_mm": round(e, 3),
        "altura_pilha_mm": round(altura_pilha, 1), "pilha_cabe_no_recipiente": bool(altura_pilha <= H * s + 1e-6),
        "erro_quantizacao_medio_mm": round(float(np.abs(dq - d)[mask].mean()), 3),
        "placas": tabela,
        "recipiente": {"volume_valido": bool(rec.is_volume), "corpos": int(rec.body_count),
                       "dimensoes_mm": [round(float(x), 1) for x in rec.extents],
                       "material_cm3": round(float(rec.volume) / 1000, 1),
                       "cavidade_cm3": round(cav_cm3, 1),
                       "parede_mm": a.parede, "fundo_mm": a.fundo, "folga_mm": a.folga},
        "parafina_estimada_mL": round(cav_cm3 - vol_total, 1),
        "faces": int(sum(p["faces"] for p in tabela)),
        "impressao": {"volume_valido": bool(ok and rec.is_volume), "corpos": pecas_total,
                      "dimensoes_mm": [round(float(x), 1) for x in ext_max],
                      "volume_cm3": round(vol_total, 1), "mesa_mm": list(a.mesa),
                      "cabe_na_mesa": cabe, "faces_degeneradas": 0,
                      "pronto": bool(ok and rec.is_volume and cabe)},
    }
    with open(os.path.join(pasta, "info.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return info


def reduzir(malha: trimesh.Trimesh, max_faces: int | None) -> trimesh.Trimesh:
    """Decimação por quádricas de erro (Garland & Heckbert, SIGGRAPH 1997).
    Cada colapso de aresta escolhe o vértice que minimiza a soma das distâncias²
    aos planos das faces originais vizinhas → remove faces primeiro onde a
    superfície é plana (topo/fundo do sólido) e preserva curvatura (bordas).
    Útil para o 'Convert Mesh' do Fusion, que fica muito lento com 10⁵ faces."""
    if not max_faces or len(malha.faces) <= max_faces:
        return malha
    r = malha.simplify_quadric_decimation(face_count=int(max_faces))
    r.merge_vertices()
    r.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(r)
    if not r.is_watertight:
        trimesh.repair.fill_holes(r)
    return r


# ─────────────────────────── impressão 3D ───────────────────────────

def preparar_impressao(malha: trimesh.Trimesh, modo: str, escala: float) -> trimesh.Trimesh:
    """Orienta, escala e apoia o modelo na mesa (z = 0), em mm.

    solido: a face plana da parede torácica (x mínimo) vai para baixo. Nessa posição
            o contorno da mama só 'afunila' para cima → praticamente sem balanços
            (overhangs), dispensando suportes, e a área de contato é grande.
    relevo: já tem base plana em z = 0.
    volume: apoia o z mínimo (face do compressor na DBT)."""
    tf = trimesh.transformations
    T = np.eye(4)
    if modo in ("solido", "composto"):
        # Ry(−90°): x → z, logo x mínimo (parede torácica) vira z mínimo
        T = tf.rotation_matrix(-np.pi / 2, [0, 1, 0]) @ T
    if escala != 1.0:
        T = tf.scale_matrix(escala) @ T
    m = malha.copy()
    m.apply_transform(T)
    lo, hi = m.bounds
    Tr = tf.translation_matrix([-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]])
    m.apply_transform(Tr)
    return m, Tr @ T  # a matriz é reaplicada às partes internas (mesmo referencial)


def verificar_impressao(m: trimesh.Trimesh, mesa_mm) -> dict:
    """Checagens que fatiadores (Cura, PrusaSlicer, Bambu Studio) exigem."""
    ext = m.extents
    r = {
        "estanque": bool(m.is_watertight),                   # sem buracos
        "normais_consistentes": bool(m.is_winding_consistent),
        "volume_valido": bool(m.is_volume),                  # fechado + orientado p/ fora
        "corpos": int(m.body_count),                         # 1 = peça única
        "faces_degeneradas": int((m.area_faces < 1e-10).sum()),
        "dimensoes_mm": [round(float(e), 1) for e in ext],
        "volume_cm3": round(float(m.volume) / 1000, 1),
        "mesa_mm": list(mesa_mm),
        "cabe_na_mesa": bool(np.all(ext <= np.asarray(mesa_mm))),
    }
    if not r["cabe_na_mesa"]:
        r["escala_max_sugerida"] = round(float(np.min(np.asarray(mesa_mm) / ext)) * 0.95, 3)
    r["pronto"] = (r["volume_valido"] and r["corpos"] == 1 and r["faces_degeneradas"] == 0
                   and r["cabe_na_mesa"])
    return r


def malha_volume(m: Mamografia, voxel_mm: float):
    """DBT: reamostra para voxels ~isotrópicos e extrai isosuperfície da mama."""
    vol = m.imagem
    esc = (m.fatia_mm / voxel_mm, m.pixel_mm[0] / voxel_mm, m.pixel_mm[1] / voxel_mm)
    vol = transform.rescale(vol, esc, anti_aliasing=True, preserve_range=True).astype(np.float32)
    nivel = filters.threshold_otsu(vol)
    vol = np.pad(vol, 2, constant_values=vol.min())  # fecha a superfície nas bordas
    verts, faces, _, _ = measure.marching_cubes(vol, level=nivel, spacing=(voxel_mm,) * 3)
    verts -= 2 * voxel_mm
    v = np.c_[verts[:, 2], -verts[:, 1], verts[:, 0]]
    return trimesh.Trimesh(v, faces[:, ::-1], process=True)


# ─────────────────────────── exportação ───────────────────────────

HTML = """<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>__TITULO__</title>
<style>html,body{margin:0;height:100%;background:#111;color:#ddd;font:13px system-ui}
#info{position:absolute;left:12px;top:10px;opacity:.8;pointer-events:none}</style>
<script type="importmap">{"imports":{
 "three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
 "three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}</script></head>
<body><div id="info">__TITULO__ · arraste para girar, roda para zoom</div><script type="module">
import * as THREE from "three";
import {OrbitControls} from "three/addons/controls/OrbitControls.js";
import {GLTFLoader} from "three/addons/loaders/GLTFLoader.js";
const r=new THREE.WebGLRenderer({antialias:true});r.setPixelRatio(devicePixelRatio);
r.setSize(innerWidth,innerHeight);document.body.appendChild(r.domElement);
const s=new THREE.Scene();s.background=new THREE.Color(0x111111);
const c=new THREE.PerspectiveCamera(40,innerWidth/innerHeight,0.1,5000);
s.add(new THREE.HemisphereLight(0xffffff,0x333344,1.6));
const dl=new THREE.DirectionalLight(0xffffff,1.4);dl.position.set(1,1,2);s.add(dl);
const ctl=new OrbitControls(c,r.domElement);ctl.enableDamping=true;
const b64="__GLB__";const bin=Uint8Array.from(atob(b64),ch=>ch.charCodeAt(0)).buffer;
new GLTFLoader().parse(bin,"",g=>{const o=g.scene;o.traverse(n=>{if(n.isMesh){
 n.geometry.computeVertexNormals();
 const nm=(n.name+" "+(n.parent?n.parent.name:"")).toLowerCase(),adip=nm.includes("adiposo");
 n.material=new THREE.MeshStandardMaterial({vertexColors:!!n.geometry.attributes.color,
 color:n.geometry.attributes.color?0xffffff:0xd9b8a8,roughness:.6,side:THREE.DoubleSide,
 transparent:adip,opacity:adip?0.22:1,depthWrite:!adip});if(adip)n.renderOrder=1;}});
 const bx=new THREE.Box3().setFromObject(o),ctr=bx.getCenter(new THREE.Vector3()),
 sz=bx.getSize(new THREE.Vector3()).length();o.position.sub(ctr);s.add(o);
 c.position.set(sz*.35,-sz*.9,sz*.7);c.up.set(0,0,1);ctl.update();});
addEventListener("resize",()=>{c.aspect=innerWidth/innerHeight;c.updateProjectionMatrix();
 r.setSize(innerWidth,innerHeight);});
(function loop(){requestAnimationFrame(loop);ctl.update();r.render(s,c);})();
</script></body></html>"""


def exportar(malha: trimesh.Trimesh | trimesh.Scene, pasta: str, titulo: str):
    os.makedirs(pasta, exist_ok=True)
    glb = malha.export(file_type="glb")              # com cores por vértice
    with open(os.path.join(pasta, "modelo.glb"), "wb") as f:
        f.write(glb)
    html = HTML.replace("__TITULO__", titulo).replace("__GLB__", base64.b64encode(glb).decode())
    with open(os.path.join(pasta, "visualizar.html"), "w", encoding="utf-8") as f:
        f.write(html)


def painel(pasta, img, mask, p, t, titulo, rot=None, voxel_mm=1.0, dz=1.0):
    n = 5 if rot is not None else 4
    fig, ax = plt.subplots(1, n, figsize=(3.5 * n, 4.5))
    lo, hi = np.percentile(img, [1, 99])
    for a, (im, nome, cm) in zip(ax, [
        (np.clip((img - lo) / (hi - lo + 1e-9), 0, 1), "original (reamostrada)", "gray"),
        (mask, "máscara da mama", "gray"),
        (p, "intensidade normalizada p(x,y)", "bone"),
        (t, "espessura t(x,y) [mm]", "viridis"),
    ]):
        h = a.imshow(im, cmap=cm)
        a.set_title(nome, fontsize=9)
        a.axis("off")
        if nome.startswith("espessura"):
            fig.colorbar(h, ax=a, fraction=0.046)
    if rot is not None:
        # corte transversal na linha central: mostra o interior (x horizontal, z vertical)
        linha = int(np.argmax(mask.sum(1)))
        ax[3].axhline(linha, color="w", lw=0.8, ls="--")
        c = rot[linha].T[::-1]  # (z, x), z para cima
        ax[4].imshow(c, cmap="magma", interpolation="nearest", aspect=dz / voxel_mm)
        ax[4].set_title("corte na linha tracejada (z × x)\n0 fora · 1 adiposo · ≥2 denso", fontsize=9)
        ax[4].axis("off")
    fig.suptitle(titulo, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(pasta, "painel.png"), dpi=110)
    plt.close(fig)


# ───────────────────────────── principal ─────────────────────────────

def processar(arq: str, a) -> dict:
    nome = os.path.splitext(os.path.basename(arq))[0]
    pasta = os.path.join(a.saida, nome)
    m = carregar(arq, a.pixel_mm)
    info = {"arquivo": arq, "formato": m.meta["Formato"], "modo": a.modo, "voxel_mm": a.voxel,
            "pixel_mm": list(m.pixel_mm), "fonte_pixel": m.meta["FontePixel"],
            "vista": m.meta["Vista"], "lateralidade": m.meta["Lateralidade"]}
    partes: dict[str, trimesh.Trimesh] = {}

    if m.eh_volume or a.modo == "volume":
        if not m.eh_volume:
            raise RuntimeError("modo 'volume' exige DICOM multi-frame (DBT).")
        malha = reduzir(malha_volume(m, a.voxel), a.max_faces)
        info["modo"] = "volume"
    else:
        img = reamostrar(m.imagem, m.pixel_mm, a.voxel)
        img, transmissao = linearizar(img)
        info["convertida_para_atenuacao"] = transmissao
        mask = segmentar_mama(img, a.limiar)
        img, mask, espelhou = orientar(img, mask)
        mlo = {"sim": True, "nao": False}.get(a.peitoral, m.meta["Vista"] == "MLO" or "mlo" in nome.lower())
        if mlo:
            peit = remover_peitoral(img, mask, a.voxel)
            mask = mask & ~peit
            abd = cortar_abdome(mask, a.voxel)
            mask = mask & ~abd
            info["abdome_cortado_cm2"] = round(float(abd.sum() * a.voxel**2 / 100), 1)
            lab_m, _n = ndi.label(mask)
            if _n > 1:  # sobras desconectadas depois do corte
                mask = lab_m == (np.argmax(np.bincount(lab_m.ravel())[1:]) + 1)
            info["peitoral_removido_cm2"] = round(float(peit.sum() * a.voxel**2 / 100), 1)
            info["mlo"] = True
        p = intensidade_normalizada(img, mask)

        H = a.espessura or m.espessura_mm
        if H is None:
            H = 50.0
            info["aviso"] = "BodyPartThickness ausente; H = 50 mm assumido"
        t = mapa_espessura(mask, H, a.voxel)
        info |= {"espessura_H_mm": H, "espelhado": espelhou,
                 "area_mama_cm2": float(mask.sum() * a.voxel**2 / 100)}

        if a.modo == "placas":
            os.makedirs(pasta, exist_ok=True)
            return processar_placas(a, info, pasta, nome, img, mask, t, H)

        rot, dz = None, a.voxel
        if a.modo == "solido":
            malha = malha_solida(t, mask, a.voxel)
            info["volume_estimado_cm3"] = float(t.sum() * a.voxel**2 / 1000)
        elif a.modo == "relevo":
            malha = malha_relevo(ndi.gaussian_filter(p, 1.0), mask, a.voxel, a.escala_relevo, a.base)
        else:  # composto
            # suaviza p: evita ilhas de 1 voxel que viram fragmentos não imprimíveis
            ps = ndi.gaussian_filter(p, 1.5) * mask
            dz = min(a.voxel, 0.5)
            rot, z0 = volume_rotulado(t, ps, mask, a.niveis, a.relevo_superficie, dz)
            malha = superficie(rot > 0, a.voxel, dz, z0)
            partes = partes_compostas(rot, a.niveis, a.voxel, dz, z0)
            vox_cm3 = a.voxel * a.voxel * dz / 1000
            vols = {nome: round(float((rot == i).sum() * vox_cm3), 1)
                    for i, nome in nomes_rotulos(a.niveis).items()}
            total = sum(vols.values())
            info |= {"niveis": a.niveis, "relevo_superficie_mm": a.relevo_superficie,
                     "volume_por_parte_cm3": vols,
                     "fracao_densa_estimada": round(1 - vols["adiposo"] / total, 3) if total else None}
        info["faces_originais"] = int(len(malha.faces))
        malha = reduzir(malha, a.max_faces)
        colorir(malha, p, a.voxel)  # depois da redução: cor amostrada nos vértices finais
        os.makedirs(pasta, exist_ok=True)
        painel(pasta, img, mask, p, t, f"{nome} · {m.meta['Vista']} {m.meta['Lateralidade']}",
               rot, a.voxel, dz)

    if partes:
        for k in partes:
            partes[k] = reduzir(partes[k], a.max_faces)
            partes[k].visual.face_colors = cor_parte(k, a.niveis)
        # visualizador: cópias reduzidas (o navegador não precisa da resolução de impressão)
        leve = {k: reduzir(v.copy(), 60000) for k, v in partes.items()}
        for k, v in leve.items():
            v.visual.face_colors = cor_parte(k, a.niveis)
        exportar(trimesh.Scene(leve), pasta, nome)
    else:
        exportar(malha, pasta, nome)

    imp, T = preparar_impressao(malha, info["modo"], a.escala_impressao)
    imp.export(os.path.join(pasta, "impressao.stl"))  # mm, orientado, apoiado em z = 0
    info |= {"vertices": int(len(malha.vertices)), "faces": int(len(malha.faces)),
             "estanque": bool(malha.is_watertight),
             "escala_impressao": a.escala_impressao,
             "impressao": verificar_impressao(imp, a.mesa)}

    if partes:
        # Partes no MESMO referencial da peça inteira: importadas juntas no fatiador
        # (PrusaSlicer/Bambu/Cura: "objeto com várias partes") encaixam sem deslocar.
        pp = os.path.join(pasta, "partes")
        os.makedirs(pp, exist_ok=True)
        cena = trimesh.Scene()
        info["partes"] = {}
        for k, v in partes.items():
            w = v.copy()
            w.apply_transform(T)
            w.export(os.path.join(pp, f"{k}.stl"))
            cena.add_geometry(w, geom_name=k, node_name=k)
            info["partes"][k] = {"faces": int(len(w.faces)), "volume_valido": bool(w.is_volume),
                                 "corpos": int(w.body_count),
                                 "volume_cm3": round(float(w.volume) / 1000, 1)}
        cena.export(os.path.join(pasta, "impressao_multimaterial.3mf"))
    with open(os.path.join(pasta, "info.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entradas", nargs="+", help="arquivos, pastas ou curingas")
    ap.add_argument("--modo", choices=["placas", "composto", "solido", "relevo", "volume"], default="placas")
    g = ap.add_argument_group("modo placas (só tecido fibroglandular; gordura = parafina no recipiente)")
    g.add_argument("--espessura-placa", dest="espessura_placa", type=float, default=2.0,
                   help="espessura de cada placa em mm (padrão 2.0)")
    g.add_argument("--max-placas", dest="max_placas", type=int, default=15,
                   help="limite de placas; se exceder, a espessura aumenta (padrão 15)")
    g.add_argument("--fracao-alvo", dest="fracao_alvo", type=float,
                   help="reescala para esta fração glandular volumétrica (ex.: 0.20)")
    g.add_argument("--peitoral", choices=["auto", "sim", "nao"], default="auto",
                   help="remove o músculo peitoral (MLO): auto detecta por vista/nome do arquivo")
    g.add_argument("--margem-borda", dest="margem_borda", type=float, default=3.0,
                   help="faixa junto à borda da imagem sem placa, em mm (artefato de borda; padrão 3)")
    g.add_argument("--camada", type=float, default=0.0,
                   help="altura de camada do fatiador em mm (ex.: 0.2): a espessura da placa vira múltiplo dela")
    g.add_argument("--suavizacao", type=float, default=0.75,
                   help="desvio-padrão da suavização do relevo, em mm (padrão 0.75; menor = mais detalhe e mais ruído)")
    g.add_argument("--limiar-fibroso", dest="limiar_fibroso", type=float, default=1.0,
                   help="ajusta o corte gordura×fibroso: <1 inclui feixes mais tênues, >1 só o mais denso")
    g.add_argument("--fator-material", dest="fator_material", type=float, default=1.0,
                   help="mm de placa por mm de tecido glandular: (μ_glandular − μ_gordura)/(μ_placa − μ_parafina). "
                        "Padrão 1 (não calibrado); ver docs/CALIBRACAO_RADIOGRAFICA.md")
    g.add_argument("--margem-pele", dest="margem_pele", type=float, default=5.0,
                   help="faixa junto à pele sem tecido denso, em mm (padrão 5)")
    g.add_argument("--traco-min", dest="traco_min", type=float, default=1.2,
                   help="menor largura de traço imprimível, em mm (padrão 1.2)")
    g.add_argument("--area-min", dest="area_min", type=float, default=25.0,
                   help="ilhas menores que isto (mm²) são descartadas (padrão 25)")
    g.add_argument("--base-placa", dest="base_placa", type=float, default=0.0,
                   help="lâmina de base que une as ilhas de cada placa, em mm (padrão 0 = sem base)")
    g.add_argument("--parede", type=float, default=3.0, help="parede do recipiente em mm (padrão 3)")
    g.add_argument("--fundo", type=float, default=2.0, help="fundo do recipiente em mm (padrão 2)")
    g.add_argument("--folga", type=float, default=0.5, help="folga entre placas e recipiente em mm (padrão 0.5)")
    ap.add_argument("--niveis", type=int, default=2,
                    help="composto: nº de materiais (2 = adiposo + denso; 3+ subdivide o denso por cinza)")
    ap.add_argument("--relevo-superficie", dest="relevo_superficie", type=float, default=3.0,
                    help="composto: altura máx. do relevo externo pelos tons, em mm (0 desliga)")
    ap.add_argument("--voxel", type=float, default=1.0, help="resolução da malha em mm (padrão 1.0)")
    ap.add_argument("--espessura", type=float, help="força H em mm (sobrepõe BodyPartThickness; "
                    "em TIFF não há essa informação e o padrão é 50 mm)")
    ap.add_argument("--pixel-mm", dest="pixel_mm", type=float,
                    help="tamanho do pixel em mm (necessário em TIFF sem calibração)")
    ap.add_argument("--limiar", type=float, default=1.0, help="multiplicador do limiar de Otsu")
    ap.add_argument("--escala-relevo", dest="escala_relevo", type=float, default=20.0,
                    help="altura máxima do relevo em mm")
    ap.add_argument("--base", type=float, default=3.0, help="espessura da base do relevo em mm")
    ap.add_argument("--escala-impressao", dest="escala_impressao", type=float, default=1.0,
                    help="fator de escala do STL (1.0 = tamanho real)")
    ap.add_argument("--mesa", type=float, nargs=3, default=[220, 220, 250], metavar=("X", "Y", "Z"),
                    help="volume útil da impressora em mm (padrão 220 220 250)")
    ap.add_argument("--max-faces", dest="max_faces", type=int,
                    help="reduz a malha para no máximo N faces (ex.: 20000 para editar no Fusion)")
    ap.add_argument("--saida", default="saida")
    a = ap.parse_args()

    arquivos = listar_dicoms(a.entradas)
    if not arquivos:
        sys.exit("Nenhum DICOM/TIFF encontrado.")
    for arq in arquivos:
        try:
            info = processar(arq, a)
            r = info["impressao"]
            dx, dy, dz = r["dimensoes_mm"]
            if info["modo"] == "placas":
                print(f"[{'OK' if r['pronto'] else 'ATENÇÃO'}] {os.path.basename(arq)}: "
                      f"{info['n_placas']} placas × {info['espessura_placa_mm']} mm · {r['corpos']} peças · "
                      f"{r['volume_cm3']} cm³ impressos · fração glandular {info['fracao_glandular_impressa']:.1%} · "
                      f"parafina ≈ {info['parafina_estimada_mL']} mL · válido={r['volume_valido']} · "
                      f"cabe_na_mesa={r['cabe_na_mesa']}")
                continue
            print(f"[{'OK' if r['pronto'] else 'ATENÇÃO'}] {os.path.basename(arq)}: "
                  f"{info['faces']} faces · {dx}×{dy}×{dz} mm · {r['volume_cm3']} cm³ · "
                  f"volume_valido={r['volume_valido']} · corpos={r['corpos']} · "
                  f"cabe_na_mesa={r['cabe_na_mesa']}"
                  + (f" (use --escala-impressao {r['escala_max_sugerida']})"
                     if not r["cabe_na_mesa"] else ""))
        except Exception as e:
            print(f"[ERRO] {os.path.basename(arq)}: {e}")
    print(f"Resultados em: {os.path.abspath(a.saida)}")


if __name__ == "__main__":
    main()
