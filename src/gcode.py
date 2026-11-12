"""Gerador de G-code (Marlin, Ender 3) para as placas: sólidos prismáticos de preenchimento 100 %.

Por que existe um gerador próprio: cada placa do modo `placas` é uma EXTRUSÃO (mesma
seção em todas as camadas), então o fatiamento é simples: a cada camada, as mesmas
paredes e o mesmo preenchimento. Não substitui o Cura/Creality Print: o G-code do
fatiador continua sendo a referência. Este gerador serve para entregar o arquivo no mesmo
formato dos que já estão no pendrive (Marlin, extrusão absoluta, start/end do Ender 3).

Estrutura de cada camada (estilo Cura):
  1. brim (só na camada 0): `--brim` mm de anéis ao redor das peças, de fora para dentro;
  2. por peça: parede interna, parede externa e preenchimento retilíneo a 45°/135°
     alternado, com espaçamento igual à largura da linha (100 %).
  Peças estreitas (< 1,5 mm) saem só de paredes, que se sobrepõem e fecham o miolo.

Extrusão: E = largura · altura · comprimento / (π d²/4), com largura 0,4 e altura 0,2
(0,03326 mm de filamento por mm, o mesmo valor que se confere nos arquivos do pendrive).
Retração de 5 mm a 45 mm/s em deslocamentos maiores que 1,5 mm.

Uso: python src/gcode.py saida/placas_ender3/lcc1/placas/placa_01.stl saida/gcode/CE3E3V2_placa01.gcode
"""
from __future__ import annotations

import argparse
import io
import math
import os

import numpy as np
import trimesh
from shapely.affinity import rotate, translate
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

START = """M140 S{bed}
M105
M190 S{bed}
M104 S{nozzle}
M105
M109 S{nozzle}
; Ender 3 Custom Start G-code
G92 E0 ; Reset Extruder
G28 ; Home all axes
G1 Z2.0 F3000 ; Move Z Axis up little to prevent scratching of Heat Bed
G1 X0.1 Y20 Z0.3 F5000.0 ; Move to start position
G1 X0.1 Y200.0 Z0.3 F1500.0 E15 ; Draw the first line
G1 X0.4 Y200.0 Z0.3 F5000.0 ; Move to side a little
G1 X0.4 Y20 Z0.3 F1500.0 E30 ; Draw the second line
G92 E0 ; Reset Extruder
G1 Z2.0 F3000 ; Move Z Axis up little to prevent scratching of Heat Bed
G1 X5 Y20 Z0.3 F5000.0 ; Move over to prevent blob squish
M82 ;absolute extrusion mode
G92 E0
G92 E0
G1 F2700 E-5
"""

END = """G91 ;Relative positioning
G1 E-2 F2700 ;Retract a bit
G1 E-2 Z0.2 F2400 ;Retract and raise Z
G1 X5 Y5 F3000 ;Wipe out
G1 Z10 ;Raise Z more
G90 ;Absolute positioning

G1 X0 Y235 ;Present print
M106 S0 ;Turn-off fan
M104 S0 ;Turn-off hotend
M140 S0 ;Turn-off bed

M84 X Y E ;Disable all steppers but Z

M82 ;absolute extrusion mode
M104 S0
;End of Gcode
"""


def secao(m: trimesh.Trimesh, z: float):
    """Seção horizontal em coordenadas reais (xy); furos pela regra par-ímpar."""
    sec = m.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
    if sec is None:
        raise RuntimeError("Seção vazia.")
    forma = Polygon()
    for p in sorted((Polygon(l[:, :2]).buffer(0) for l in sec.discrete if len(l) >= 4), key=lambda q: -q.area):
        forma = forma.symmetric_difference(p)
    return forma.simplify(0.02)


def poligonos(geom):
    return [p for p in getattr(geom, "geoms", [geom]) if p.geom_type == "Polygon" and not p.is_empty]


def aneis(geom, so_exterior=False):
    """Contornos (exterior e furos) como arrays Nx2."""
    out = []
    for p in poligonos(geom):
        out.append(np.asarray(p.exterior.coords))
        if not so_exterior:
            out += [np.asarray(i.coords) for i in p.interiors]
    return out


def linhas_de_preenchimento(regiao, ang_graus: float, esp: float):
    """Segmentos paralelos (espaçamento `esp`) dentro da região, em zigue-zague."""
    if regiao.is_empty:
        return []
    c = np.array([regiao.centroid.x, regiao.centroid.y])
    minx, miny, maxx, maxy = regiao.bounds
    R = math.hypot(maxx - minx, maxy - miny) / 2 + 1
    a = math.radians(ang_graus)
    d, n = np.array([math.cos(a), math.sin(a)]), np.array([-math.sin(a), math.cos(a)])
    segs = []
    for o in np.arange(-R, R, esp):
        inter = LineString([c + n * o - d * R, c + n * o + d * R]).intersection(regiao)
        for g in getattr(inter, "geoms", [inter]):
            if g.geom_type == "LineString" and g.length > 0.3:
                segs.append(np.asarray(g.coords))
    return segs


def ordenar_vizinho_mais_proximo(segs, x, y):
    """Ordem gulosa: sempre o segmento cuja ponta está mais perto (inverte se preciso).
    Encurta os deslocamentos e, com isso, as retrações em formas ramificadas."""
    if not segs:
        return []
    p0 = np.array([s[0] for s in segs])
    p1 = np.array([s[-1] for s in segs])
    usados = np.zeros(len(segs), dtype=bool)
    out, cur = [], np.array([x, y])
    for _ in range(len(segs)):
        d0 = ((p0 - cur) ** 2).sum(1)
        d1 = ((p1 - cur) ** 2).sum(1)
        d0[usados] = d1[usados] = np.inf
        i0, i1 = int(np.argmin(d0)), int(np.argmin(d1))
        if d0[i0] <= d1[i1]:
            i, s = i0, segs[i0]
        else:
            i, s = i1, segs[i1][::-1]
        usados[i] = True
        out.append(s)
        cur = s[-1]
    return out


class Escritor:
    def __init__(self, f, largura, altura, diametro, retracao, vel_retracao):
        self.f = f
        self.K = largura * altura / (math.pi * (diametro / 2) ** 2)  # mm de filamento por mm de percurso
        self.E, self.x, self.y, self.z = -5.0, 0.0, 0.0, 0.0        # a partida termina retraída (E-5)
        self.retr, self.vel_retr = retracao, vel_retracao
        self.retraido, self.filamento_mm, self.t = True, 0.0, 0.0
        self.livre = None  # região onde se pode viajar sem retrair (dentro da peça já impressa)

    def viaja(self, x, y, z=None):
        dist = math.hypot(x - self.x, y - self.y)
        # "combing": deslocamento inteiro dentro da própria peça não precisa de retração
        dentro = self.livre is not None and dist > 0 and self.livre.contains(LineString([(self.x, self.y), (x, y)]))
        if dist > 1.5 and not dentro and not self.retraido:
            self.E -= self.retr
            self.f.write(f"G1 F{self.vel_retr} E{self.E:.5f}\n")
            self.retraido = True
        zs = f" Z{z:.3f}" if z is not None and abs(z - self.z) > 1e-9 else ""
        self.f.write(f"G0 F6000 X{x:.3f} Y{y:.3f}{zs}\n")
        self.t += dist / 100.0
        self.x, self.y = x, y
        if z is not None:
            self.z = z

    def desretrai(self):
        if self.retraido:
            self.E += self.retr
            self.f.write(f"G1 F{self.vel_retr} E{self.E:.5f}\n")
            self.retraido = False

    def percurso(self, xy: np.ndarray, vel_mm_s: float, z: float):
        self.viaja(xy[0][0], xy[0][1], z)
        self.desretrai()
        F, primeiro = int(round(vel_mm_s * 60)), True
        for x, y in xy[1:]:
            L = math.hypot(x - self.x, y - self.y)
            if L < 1e-6:
                continue
            self.E += L * self.K
            self.filamento_mm += L * self.K
            self.f.write(f"G1 {'F%d ' % F if primeiro else ''}X{x:.3f} Y{y:.3f} E{self.E:.5f}\n")
            primeiro = False
            self.t += L / vel_mm_s
            self.x, self.y = x, y


def a_partir_do_mais_proximo(ring: np.ndarray, x, y):
    pts = ring[:-1]
    i = int(np.argmin((pts[:, 0] - x) ** 2 + (pts[:, 1] - y) ** 2))
    pts = np.roll(pts, -i, axis=0)
    return np.vstack([pts, pts[:1]])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stl")
    ap.add_argument("saida")
    ap.add_argument("--camada", type=float, default=0.2)
    ap.add_argument("--linha", type=float, default=0.4, help="largura da linha, mm")
    ap.add_argument("--diametro", type=float, default=1.75)
    ap.add_argument("--bico", type=int, default=200, help="°C (PLA: 200, como nos arquivos do pendrive)")
    ap.add_argument("--mesa", type=int, default=60, help="°C")
    ap.add_argument("--brim", type=float, default=4.0, help="mm (0 = sem)")
    ap.add_argument("--leito", type=float, nargs=2, default=[220, 220])
    ap.add_argument("--rotacao", default="auto",
                    help="giro da peça no plano da mesa, em graus; 'auto' só gira se não couber (padrão)")
    a = ap.parse_args()

    m = trimesh.load(a.stl)
    esp = float(m.extents[2])
    n = int(round(esp / a.camada))
    if abs(esp - n * a.camada) > 0.05:
        raise SystemExit(f"Espessura {esp:.3f} mm não é múltiplo de {a.camada} mm (use --camada no gerador de placas).")
    P = secao(m, esp / 2)
    ang = 0.0
    if a.rotacao != "auto":
        ang = float(a.rotacao)
    else:  # peça maior que a mesa: procura o ângulo (no plano) que a faz caber, com folga para o brim
        lim = (a.leito[0] - 2 * (a.brim + 1), a.leito[1] - 2 * (a.brim + 1))
        c0 = unary_union(poligonos(P)).centroid

        def lados(ang_):
            b = rotate(P, ang_, origin=c0).bounds
            return b[2] - b[0], b[3] - b[1]
        if lados(0)[0] > lim[0] or lados(0)[1] > lim[1]:
            cabem = [(max(lados(g)), g) for g in range(1, 180) if lados(g)[0] <= lim[0] and lados(g)[1] <= lim[1]]
            if not cabem:
                raise SystemExit(f"A peça não cabe na mesa {a.leito} mesmo girada; divida ou reduza a placa.")
            ang = float(min(cabem)[1])
    if ang:
        P = rotate(P, ang, origin=unary_union(poligonos(P)).centroid)
    minx, miny, maxx, maxy = P.bounds   # centraliza a peça na mesa
    P = translate(P, a.leito[0] / 2 - (minx + maxx) / 2, a.leito[1] / 2 - (miny + maxy) / 2)
    ilhas = sorted(poligonos(P), key=lambda g: (round(g.centroid.y / 20), g.centroid.x))
    W, H = a.linha, a.camada
    uniao = unary_union(ilhas)

    buf = io.StringIO()
    w = Escritor(buf, W, H, a.diametro, 5.0, 2700)
    for camada in range(n):
        z = round((camada + 1) * H, 3)
        buf.write(f";LAYER:{camada}\n")
        if camada == 0:
            buf.write("M107\n")
        elif camada <= 3:
            buf.write(f"M106 S{[85, 170, 255][camada - 1]}\n")
        v_int, v_ext, v_inf = (20, 20, 20) if camada == 0 else (40, 30, 40)  # mm/s; 1ª camada lenta
        if camada == 0 and a.brim > 0:
            buf.write(";TYPE:SKIRT\n")
            for k in range(int(round(a.brim / W)) - 1, -1, -1):  # de fora para dentro
                anel = uniao.buffer((k + 0.5) * W, join_style=1, resolution=6).simplify(0.02)
                for ring in aneis(anel, so_exterior=True):
                    w.percurso(a_partir_do_mais_proximo(ring, w.x, w.y), v_int, z)
        for ilha in ilhas:
            externo = ilha.buffer(-0.5 * W, join_style=1, resolution=6)
            if externo.is_empty:
                continue
            interno = ilha.buffer(-1.5 * W, join_style=1, resolution=6)
            w.livre = ilha.buffer(-0.25, join_style=1, resolution=4)
            for tipo, geom, vel in (("WALL-INNER", interno, v_int), ("WALL-OUTER", externo, v_ext)):
                if geom.is_empty:
                    continue
                buf.write(f";TYPE:{tipo}\n")
                for ring in aneis(geom.simplify(0.02)):
                    w.percurso(a_partir_do_mais_proximo(ring, w.x, w.y), vel, z)
            regiao = ilha.buffer(-(2 * W - 0.1 * W), join_style=1, resolution=6)  # 10 % de sobreposição nas paredes
            if not regiao.is_empty:
                buf.write(";TYPE:FILL\n")
                segs = linhas_de_preenchimento(regiao, 45 if camada % 2 == 0 else 135, W)
                for seg in ordenar_vizinho_mais_proximo(segs, w.x, w.y):
                    w.percurso(seg, v_inf, z)
            w.livre = None

    bx0, by0, bx1, by1 = (uniao.buffer(a.brim) if a.brim > 0 else uniao).bounds
    os.makedirs(os.path.dirname(os.path.abspath(a.saida)), exist_ok=True)
    with open(a.saida, "w", newline="\n") as f:
        f.write(f";ROTATION_Z_DEG:{ang:g}\n;FLAVOR:Marlin\n;TIME:{int(w.t + 90)}\n;Filament used: {w.filamento_mm / 1000:.5f}m\n")
        f.write(f";Layer height: {H}\n;MINX:{bx0:.3f}\n;MINY:{by0:.3f}\n;MINZ:{H}\n")
        f.write(f";MAXX:{bx1:.3f}\n;MAXY:{by1:.3f}\n;MAXZ:{n * H:.3f}\n;TARGET_MACHINE.NAME:Creality Ender-3\n")
        f.write(";Generated by src/gcode.py (gerador proprio, NAO e o Cura). Conferir a previa antes de imprimir.\n")
        f.write(START.format(bed=a.mesa, nozzle=a.bico))
        f.write(f";LAYER_COUNT:{n}\n")
        f.write(buf.getvalue())
        f.write(f"G1 F2700 E{w.E - 5.0:.5f}\nM140 S0\nM107\n")
        f.write(END)
    print(f"{a.saida}: {n} camadas de {H} mm · {len(ilhas)} peça(s) · centralizada em "
          f"({a.leito[0] / 2:g}, {a.leito[1] / 2:g}) · filamento {w.filamento_mm / 1000:.2f} m · ~{w.t / 3600:.1f} h")


if __name__ == "__main__":
    main()
