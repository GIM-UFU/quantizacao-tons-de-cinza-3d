"""Estimativa de filamento (FDM) para uma placa em STL.

Uso: python src/filamento.py saida/placas_ender3/lcc1/placas/placa_01.stl [opções]

O que entra na conta (valores em mm, mm³ e g):
  1. Modelo: com 100 % de preenchimento, o volume extrudado é o volume do STL
     (× fluxo). Paredes e preenchimento somam o mesmo sólido.
  2. Brim: faixa de `--brim` mm ao redor da base, impressa só na 1ª camada. Calculada
     pela pegada real da peça (seção em z = meia camada), não por aproximação de caixa.
  3. Skirt: `--skirt` voltas ao redor de tudo, a `--skirt-dist` mm de distância.
  4. Purga: filamento gasto na linha de preparo do início do G-code (valor típico,
     depende do perfil do fatiador).
  5. Reserva: margem para perdas e reimpressão parcial.

Limites: o fatiador (Cura, Creality Print) é quem dá o número final; esta estimativa
costuma ficar próxima dele (cerca de ±5 %), e não inclui refugo de falhas.
"""
from __future__ import annotations

import argparse
import math

import numpy as np
import trimesh
from shapely.geometry import Polygon
from shapely.ops import unary_union


def pegada(m: trimesh.Trimesh, h_camada: float):
    """União das seções da peça em z = meia camada (a pegada na mesa)."""
    sec = m.section(plane_origin=[0, 0, 0.5 * h_camada], plane_normal=[0, 0, 1])
    if sec is None:
        raise RuntimeError("Seção vazia: a peça está apoiada em z = 0?")
    planar, _ = sec.to_2D()
    # cada laço fechado vira um polígono; furos entram por "ou exclusivo" (regra par-ímpar),
    # o que dispensa a árvore de inclusão do trimesh (que exigiria o pacote rtree)
    laços = [Polygon(l).buffer(0) for l in planar.discrete if len(l) >= 4]
    forma = Polygon()
    for p in sorted(laços, key=lambda q: -q.area):
        forma = forma.symmetric_difference(p)
    return forma


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stl")
    ap.add_argument("--diametro", type=float, default=1.75, help="diâmetro do filamento, mm (Ender 3: 1.75)")
    ap.add_argument("--densidade", type=float, default=1.24, help="g/cm³ (PLA = 1,24, Bliznakova 2024, Tab. 2)")
    ap.add_argument("--fluxo", type=float, default=100.0, help="fluxo do fatiador, %%")
    ap.add_argument("--camada", type=float, default=0.2, help="altura de camada, mm")
    ap.add_argument("--linha", type=float, default=0.4, help="largura da linha, mm")
    ap.add_argument("--brim", type=float, default=4.0, help="largura do brim, mm (0 = sem)")
    ap.add_argument("--skirt", type=int, default=2, help="voltas de skirt (0 = sem)")
    ap.add_argument("--skirt-dist", dest="skirt_dist", type=float, default=3.0)
    ap.add_argument("--purga", type=float, default=0.5, help="g gastos na linha de preparo (valor típico)")
    ap.add_argument("--reserva", type=float, default=10.0, help="margem de segurança, %%")
    ap.add_argument("--preco-kg", dest="preco", type=float, help="preço do filamento por kg (opcional)")
    ap.add_argument("--vel", type=float, default=50.0, help="velocidade de impressão, mm/s (só para o tempo)")
    a = ap.parse_args()

    m = trimesh.load(a.stl)
    V = float(m.volume)                                  # mm³
    area_fil = math.pi * (a.diametro / 2) ** 2           # mm²
    g_por_mm3 = a.densidade / 1000.0

    base = pegada(m, a.camada)
    area_base = base.area
    v_brim = 0.0
    if a.brim > 0:
        brim = base.buffer(a.brim, join_style=1).difference(base)
        v_brim = brim.area * a.camada
    v_skirt = 0.0
    if a.skirt > 0:
        ctorno = (base.buffer(a.brim) if a.brim > 0 else base).convex_hull.buffer(a.skirt_dist)
        v_skirt = a.skirt * ctorno.length * a.linha * a.camada

    v_modelo = V * a.fluxo / 100.0
    v_total = v_modelo + v_brim + v_skirt
    g_purga = a.purga

    def g(v):
        return v * g_por_mm3

    modelo_g = g(v_modelo)
    provavel_g = g(v_total) + g_purga
    reserva_g = provavel_g * (1 + a.reserva / 100.0)

    def comp_m(grama):  # comprimento de filamento, em metros
        return grama / g_por_mm3 / area_fil / 1000.0

    # tempo: fluxo volumétrico = vel · linha · camada; +35 % de deslocamentos, retrações e 1ª camada lenta
    t_h = (v_total / (a.vel * a.linha * a.camada)) * 1.35 / 3600.0

    print(f"Peça: {a.stl}")
    print(f"  dimensões      : {' × '.join(f'{x:.2f}' for x in m.extents)} mm · {m.body_count} peça(s) · "
          f"pegada {area_base / 100:.1f} cm²")
    print(f"  volume (STL)   : {V / 1000:.2f} cm³  →  fluxo {a.fluxo:g} %: {v_modelo / 1000:.2f} cm³")
    print(f"  brim {a.brim:g} mm     : {v_brim / 1000:.2f} cm³ · skirt {a.skirt} volta(s): {v_skirt / 1000:.3f} cm³")
    print(f"  filamento {a.diametro:g} mm, ρ = {a.densidade:g} g/cm³\n")
    print(f"{'cenário':<34}{'massa':>9}{'comprimento':>14}{'% de 1 kg':>11}" + (f"{'custo':>10}" if a.preco else ""))
    for nome, grama in (("só o modelo (100 % de preench.)", modelo_g),
                        ("modelo + brim + skirt + purga", provavel_g),
                        (f"com reserva de {a.reserva:g} %", reserva_g)):
        linha = f"{nome:<34}{grama:>7.1f} g{comp_m(grama):>12.2f} m{grama / 10:>10.1f}%"
        if a.preco:
            linha += f"{grama / 1000 * a.preco:>9.2f}"
        print(linha)
    print(f"\ntempo de impressão (ordem de grandeza, {a.vel:g} mm/s): ≈ {t_h:.1f} h")


if __name__ == "__main__":
    main()
