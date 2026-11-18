"""Cálculo de atenuação linear e do fator de material (mm de placa por mm de tecido glandular).

Por que isto existe
-------------------
Na mamografia, o contraste entre tecido glandular e gordura vem da diferença de atenuação:
    contraste_tecido = (mu_gland - mu_gordura) * t_glandular
No phantom, a "gordura" é parafina e o "glandular" é a placa impressa:
    contraste_placa  = (mu_placa - mu_parafina) * e_placa
Para a placa reproduzir o contraste do tecido, e_placa/t_glandular deve ser
    fator_material = (mu_gland - mu_gordura) / (mu_placa - mu_parafina)        [mm de placa / mm de tecido]
onde mu = coeficiente de atenuação LINEAR (cm^-1) = rho (g/cm³) · (mu/rho) (cm²/g), na energia de interesse.

O coeficiente mássico (mu/rho) de cada material é obtido por mistura (regra das frações de massa):
    mu/rho = soma_i  w_i · (mu/rho)_i
com os (mu/rho)_i elementares da base Elam (pacote xraydb). Conferido contra a tabela NIST XCOM/ICRU-44
(tecido adiposo e mama, 15–40 keV): concorda até o 4º dígito (ver tests/testar_atenuacao.py).

Limites (importante para a apresentação)
----------------------------------------
* Cálculo MONOENERGÉTICO: o feixe real é polienergético e endurece ao atravessar o material; use a
  energia efetiva do seu espectro (tipicamente ~18–25 keV) ou veja a varredura em energia.
* As composições dos tecidos são as do ICRU-44 (NIST). "Tecido mamário ICRU-44" é a composição de
  referência da mama; o fibroglandular 100 % de outros autores pode diferir. As composições são editáveis.
* PETG/ABS e a resina são aproximações (composição de repetição); para o filamento exato, meça (ver
  docs/CALIBRACAO_RADIOGRAFICA.md) ou informe os mu diretamente (--mu-placa).

Uso:  python src/atenuacao.py --energia 20 --placa PLA
      python src/atenuacao.py --varredura 15 40 --placa PLA --png saida/atenuacao.png
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass

import xraydb

# número atômico → símbolo
_SIMBOLO = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 15: "P", 16: "S", 17: "Cl"}


@dataclass(frozen=True)
class Material:
    nome: str
    massa: dict          # {Z: fração de massa}
    densidade: float     # g/cm³
    fonte: str


def _de_formula(formula: dict) -> dict:
    """{'C': 3, 'H': 4, 'O': 2} (átomos por unidade) → frações de massa por Z."""
    pesos = {el: n * xraydb.atomic_mass(el) for el, n in formula.items()}
    tot = sum(pesos.values())
    return {xraydb.atomic_number(el): p / tot for el, p in pesos.items()}


MATERIAIS: dict[str, Material] = {
    # ICRU-44 via NIST (physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html)
    "gordura": Material("Tecido adiposo (ICRU-44)",
                        {1: .114, 6: .598, 7: .007, 8: .278, 11: .001, 16: .001, 17: .001}, 0.95,
                        "ICRU-44 / NIST"),
    "glandular": Material("Tecido mamário (ICRU-44)",
                          {1: .106, 6: .332, 7: .030, 8: .527, 11: .001, 15: .001, 16: .002, 17: .001}, 1.02,
                          "ICRU-44 / NIST"),
    # filamentos: unidade de repetição do polímero (aproximação)
    "PLA": Material("PLA (C3H4O2)n", _de_formula({"C": 3, "H": 4, "O": 2}), 1.24, "unidade de repetição; 1,24 g/cm³"),
    "PETG": Material("PETG ≈ PET (C10H8O4)n", _de_formula({"C": 10, "H": 8, "O": 4}), 1.27,
                     "aproximação: PET; glicol-modificado difere pouco"),
    "ABS": Material("ABS ≈ C15H17N", _de_formula({"C": 15, "H": 17, "N": 1}), 1.05,
                    "aproximação: copolímero acrilonitrila-butadieno-estireno"),
    # "gordura" do phantom
    "parafina": Material("Parafina (C25H52)", _de_formula({"C": 25, "H": 52}), 0.82,
                         "alcano; densidade 0,82 g/cm³ (literatura do projeto) — ajuste com --dens-parafina"),
    # cera de abelha: componente principal palmitato de miricila C46H92O2 (aproximação: ignora ~10–14 % de
    # hidrocarbonetos e ácidos livres); densidade 0,95–0,96 g/cm³ (tabelas de fornecedores/catálogos)
    "cera_abelha": Material("Cera de abelha ≈ C46H92O2", _de_formula({"C": 46, "H": 92, "O": 2}), 0.96,
                            "aproximação: palmitato de miricila; 0,96 g/cm³"),
}


def mu_massico(mat: Material, energia_kev: float) -> float:
    """mu/rho total (com espalhamento coerente), cm²/g."""
    return sum(w * xraydb.mu_elam(_SIMBOLO[z], energia_kev * 1000.0, kind="total") for z, w in mat.massa.items())


def mu_linear(mat: Material, energia_kev: float, densidade: float | None = None) -> float:
    """mu = rho · mu/rho, em cm^-1."""
    return (densidade if densidade is not None else mat.densidade) * mu_massico(mat, energia_kev)


def fator_material(energia_kev: float, placa: str = "PLA", dens_placa: float | None = None,
                   dens_parafina: float | None = None, mu_placa: float | None = None,
                   glandular: str = "glandular", gordura: str = "gordura", fundo: str = "parafina") -> dict:
    """Devolve o fator e todas as parcelas do cálculo (para registrar na execução).
    `fundo` = material que imita a gordura ao redor das placas ("parafina" ou "cera_abelha");
    `dens_parafina` sobrescreve a densidade do material de fundo, seja qual for."""
    mg = mu_linear(MATERIAIS[glandular], energia_kev)
    mf = mu_linear(MATERIAIS[gordura], energia_kev)
    mpar = mu_linear(MATERIAIS[fundo], energia_kev, dens_parafina)
    mp = mu_placa if mu_placa is not None else mu_linear(MATERIAIS[placa], energia_kev, dens_placa)
    if mp <= mpar:
        raise ValueError("mu da placa deve ser maior que o da parafina nessa energia.")
    return {
        "energia_keV": energia_kev, "placa": placa if mu_placa is None else "mu informado", "fundo": fundo,
        "mu_glandular_cm-1": mg, "mu_gordura_cm-1": mf, "mu_placa_cm-1": mp, "mu_parafina_cm-1": mpar,
        "dens_placa": dens_placa if dens_placa is not None else MATERIAIS.get(placa, MATERIAIS["PLA"]).densidade,
        "dens_parafina": dens_parafina if dens_parafina is not None else MATERIAIS[fundo].densidade,
        "contraste_tecido_cm-1": mg - mf, "contraste_placa_cm-1": mp - mpar,
        "fator_material": (mg - mf) / (mp - mpar),
    }


def relatorio(r: dict) -> str:
    return (f"Energia {r['energia_keV']:g} keV · placa {r['placa']}\n"
            f"  mu glandular = {r['mu_glandular_cm-1']:.4f} cm^-1   mu gordura  = {r['mu_gordura_cm-1']:.4f} cm^-1\n"
            f"  mu placa     = {r['mu_placa_cm-1']:.4f} cm^-1   mu {r['fundo']} = {r['mu_parafina_cm-1']:.4f} cm^-1\n"
            f"  contraste do tecido (mu_g - mu_f) = {r['contraste_tecido_cm-1']:.4f} cm^-1\n"
            f"  contraste da placa  (mu_p - mu_par) = {r['contraste_placa_cm-1']:.4f} cm^-1\n"
            f"  => fator_material = {r['fator_material']:.3f} mm de placa por mm de tecido glandular")


def varredura(e0: float, e1: float, passo: float = 1.0, **kw) -> list[dict]:
    n = int(round((e1 - e0) / passo)) + 1
    return [fator_material(e0 + i * passo, **kw) for i in range(n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--energia", type=float, default=20.0, help="energia efetiva em keV (padrão 20)")
    ap.add_argument("--placa", choices=["PLA", "PETG", "ABS"], default="PLA")
    ap.add_argument("--dens-placa", type=float, help="densidade real da placa impressa em g/cm³ (infill/vazios!)")
    ap.add_argument("--fundo", choices=["parafina", "cera_abelha"], default="parafina",
                    help="material que imita a gordura ao redor das placas (padrão parafina)")
    ap.add_argument("--dens-parafina", type=float, help="densidade do material de fundo em g/cm³ (parafina 0,82; cera 0,96)")
    ap.add_argument("--mu-placa", type=float, help="mu linear medido da placa (cm^-1); ignora o cálculo")
    ap.add_argument("--varredura", type=float, nargs=2, metavar=("E0", "E1"), help="tabela de E0 a E1 keV")
    ap.add_argument("--png", help="salva o gráfico da varredura")
    a = ap.parse_args()
    kw = dict(placa=a.placa, dens_placa=a.dens_placa, dens_parafina=a.dens_parafina, mu_placa=a.mu_placa, fundo=a.fundo)
    if a.varredura:
        linhas = varredura(a.varredura[0], a.varredura[1], **kw)
        print("keV   mu_gland  mu_gord  mu_placa  mu_paraf  fator")
        for r in linhas:
            print(f"{r['energia_keV']:4g}  {r['mu_glandular_cm-1']:8.4f} {r['mu_gordura_cm-1']:8.4f} "
                  f"{r['mu_placa_cm-1']:8.4f} {r['mu_parafina_cm-1']:8.4f}  {r['fator_material']:.3f}")
        if a.png:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            e = [r["energia_keV"] for r in linhas]
            for k, rot in (("mu_glandular_cm-1", "glandular"), ("mu_gordura_cm-1", "gordura"),
                           ("mu_placa_cm-1", a.placa), ("mu_parafina_cm-1", "parafina")):
                ax[0].plot(e, [r[k] for r in linhas], label=rot)
            ax[0].set_yscale("log"); ax[0].set_xlabel("energia (keV)"); ax[0].set_ylabel("mu (cm$^{-1}$)")
            ax[0].legend(); ax[0].set_title("atenuação linear")
            ax[1].plot(e, [r["fator_material"] for r in linhas], color="C3")
            ax[1].set_xlabel("energia (keV)"); ax[1].set_ylabel("mm de placa / mm de tecido")
            ax[1].set_title("fator de material"); ax[1].grid(alpha=.3)
            fig.tight_layout(); fig.savefig(a.png, dpi=130)
            print("gráfico:", a.png)
    else:
        print(relatorio(fator_material(a.energia, **kw)))


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
