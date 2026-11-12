"""Valida a saída do mammoreplicator contra o gabarito do tests/sintetico_raw.py.

Uso: python tests/validar_mammoreplicator.py <arquivo.ply> <gabarito.npz> <params.json>
"""
import json
import math
import sys

import numpy as np
import trimesh

ply, gab, prm = sys.argv[1:4]
P = json.load(open(prm))
G = np.load(gab)
b = P["output_binning"]

# ── 1. integridade do arquivo PLY ──
with open(ply) as f:
    linhas = f.read().splitlines()
fim = linhas.index("end_header")
nv = int(next(l for l in linhas if l.startswith("element vertex")).split()[-1])
nf = int(next(l for l in linhas if l.startswith("element face")).split()[-1])
corpo = linhas[fim + 1:]
faces_escritas = sum(1 for l in corpo[nv:] if l.startswith("3 "))
extras = len(corpo) - nv - faces_escritas
print("▸ PLY")
print(f"   vértices declarados {nv} · escritos {sum(1 for l in corpo[:nv] if len(l.split()) == 3)}")
print(f"   faces declaradas {nf} · escritas {faces_escritas} · linhas extras após as faces: {extras}")

V = np.array([[float(x) for x in l.split()] for l in corpo[:nv]])
F = np.array([[int(x) for x in l.split()[1:]] for l in corpo[nv:nv + faces_escritas]])
print(f"   índices de face fora do intervalo: {int((F >= nv).sum() + (F < 0).sum())}")
m = trimesh.Trimesh(V, F[(F < nv).all(1)], process=True)
print(f"   estanque={m.is_watertight} · normais_consistentes={m.is_winding_consistent} · "
      f"corpos={m.body_count} · faces degeneradas={int((m.area_faces < 1e-9).sum())}")
print(f"   dimensões (mm) = {np.round(m.extents, 1).tolist()}")

# ── 2. alturas das colunas × atenuação física ──
nr, nc = (P["rows"][1] - P["rows"][0]) // b, (P["columns"][1] - P["columns"][0]) // b
C = V.reshape(-1, 8, 3)
h = np.linalg.norm(C[:, 4] - C[:, 0], axis=1).reshape(nr, nc)  # mm (vértices de topo − base)

def binar(a):
    return a[:nr * b, :nc * b].reshape(nr, b, nc, b).mean((1, 3))

t = binar(G["t"]) / 10  # cm
f = binar(G["f"])
mask = binar(G["mask"].astype(float)) > 0.999
A = t * ((1 - f) * 0.46 + f * 0.80)  # atenuação verdadeira ∫μ dz (adimensional)

slope_teo = 10 * P["compression_thickness"] / (-math.log(P["average_breast_intensity"] / P["air_intensity"])) / 1.284
coef = np.polyfit(A[mask], h[mask], 1)
r2 = np.corrcoef(A[mask], h[mask])[0, 1] ** 2
print("\n▸ ALTURA × ATENUAÇÃO VERDADEIRA (pixels de mama)")
print(f"   h = {coef[0]:.3f}·A + {coef[1]:.3f} mm   (inclinação teórica {slope_teo:.3f}, intercepto 0)")
print(f"   R² = {r2:.5f}")
print(f"   altura: média {h[mask].mean():.1f} mm · máx {h.max():.1f} mm · ar (fora) média {h[~mask].mean():.2f} mm")

# ── 3. recuperação da fração glandular onde t = espessura plena ──
cheio = mask & (t >= P["compression_thickness"] - 1e-3)
A_est = h / slope_teo
f_est = (A_est / t - 0.46) / (0.80 - 0.46)
err = f_est[cheio] - f[cheio]
print("\n▸ FRAÇÃO GLANDULAR RECUPERADA (região de espessura plena)")
print(f"   correlação r = {np.corrcoef(f[cheio], f_est[cheio])[0, 1]:.4f} · "
      f"erro médio {err.mean():+.4f} · desvio {err.std():.4f}")
print(f"   f média verdadeira {f[cheio].mean():.3f} · estimada {f_est[cheio].mean():.3f}")
