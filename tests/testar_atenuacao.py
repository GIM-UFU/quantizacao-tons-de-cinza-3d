"""Confere src/atenuacao.py contra a tabela NIST XCOM (ICRU-44) e checa propriedades físicas.

Valores de referência: NIST, Tissue Adipose / Breast (ICRU-44), mu/rho total em cm²/g
(physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/adipose.html e breast.html).
Uso: python tests/testar_atenuacao.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import atenuacao as at  # noqa: E402

NIST = {  # keV: (adiposo, mama)
    15: (1.083, 1.378), 20: (0.5677, 0.6889), 30: (0.3063, 0.3403), 40: (0.2396, 0.2530),
}
falhas = 0
for e, (ra, rb) in NIST.items():
    ma, mb = at.mu_massico(at.MATERIAIS["gordura"], e), at.mu_massico(at.MATERIAIS["glandular"], e)
    erro = max(abs(ma / ra - 1), abs(mb / rb - 1))
    ok = erro < 0.005
    falhas += not ok
    print(f"{e:3d} keV  adiposo {ma:.4f} (NIST {ra})  mama {mb:.4f} (NIST {rb})  erro máx {erro:.2%}  {'OK' if ok else 'FALHOU'}")

# propriedades: mu cai com a energia; fator positivo; fator = razão dos contrastes
r = at.fator_material(20)
ok = abs(r["fator_material"] - r["contraste_tecido_cm-1"] / r["contraste_placa_cm-1"]) < 1e-12
falhas += not ok
print("fator = razão dos contrastes:", "OK" if ok else "FALHOU")
ok = all(at.mu_massico(at.MATERIAIS["PLA"], e1) > at.mu_massico(at.MATERIAIS["PLA"], e2) for e1, e2 in ((15, 20), (20, 30), (30, 40)))
falhas += not ok
print("mu/rho decresce com a energia:", "OK" if ok else "FALHOU")
ok = 0.2 < r["fator_material"] < 1.0
falhas += not ok
print(f"fator PLA a 20 keV = {r['fator_material']:.3f} (faixa plausível 0,2–1,0):", "OK" if ok else "FALHOU")
# densidade menor da placa (vazios) → fator maior
ok = at.fator_material(20, dens_placa=1.0)["fator_material"] > r["fator_material"]
falhas += not ok
print("placa menos densa exige mais espessura:", "OK" if ok else "FALHOU")
sys.exit(1 if falhas else 0)
