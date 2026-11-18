"""Gera a planilha de viabilidade (xlsx) a partir de saida/05_comparacoes/quantizacao_etapas/quantizacao.csv.

Abas: Recomendação (ranking com pesos editáveis), Configurações (21 = 7 etapas × 3 limiares, média das 8 imagens),
Por imagem (dados brutos) e Premissas. O escore é calculado por FÓRMULAS do Excel: mude os pesos e o ranking recalcula.

Uso: python src/planilha_viabilidade.py [saida.xlsx]
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from collections import defaultdict

from openpyxl import Workbook
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAIDA = os.path.join(RAIZ, "saida")
CSV = os.path.join(SAIDA, "05_comparacoes", "quantizacao_todas_etapas", "quantizacao.csv")
ROTULO = {"01": "Antes · fator 1,0 · 3,4 mm", "02": "Depois · atenuação · 3,4 mm", "03": "Depois · atenuação · 1,0 mm",
          "04": "Depois · atenuação · 1,4 mm", "06": "Depois · atenuação · 2,0 mm",
          "07": "Cera de abelha · 1,0 mm", "08": "Cera de abelha · 2,0 mm"}
ALVO_FG = 0.137  # limite inferior da faixa 13,7–25,6 % (Yaffe 2009 via Bliznakova 2024)

AZUL = PatternFill("solid", fgColor="1F3864")
CINZA = PatternFill("solid", fgColor="F2F2F2")
VERDE = PatternFill("solid", fgColor="E2F0D9")
AMAR = PatternFill("solid", fgColor="FFF2CC")
FINA = Side(style="thin", color="BFBFBF")
BORDA = Border(left=FINA, right=FINA, top=FINA, bottom=FINA)


def num(s):
    return float(str(s).replace(",", "."))


def cabecalho(ws, linha, titulos, larguras=None):
    for j, t in enumerate(titulos, 1):
        c = ws.cell(linha, j, t)
        c.font = Font(bold=True, color="FFFFFF")
        c.fill = AZUL
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = BORDA
    if larguras:
        for j, w in enumerate(larguras, 1):
            ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[linha].height = 42


def carregar():
    linhas = list(csv.DictReader(open(CSV, encoding="utf-8")))
    for r in linhas:
        for k in ("fator", "espessura_mm", "n_placas", "pecas", "filamento_cm3", "fracao_glandular", "degrau_tecido_mm",
                  "quant_medio_mm", "quant_rmse_mm", "quant_p95_mm", "cadeia_medio_mm"):
            r[k] = num(r[k])
        r["cod"] = r["etapa"][:2]
        r["limiar_n"] = num(r["limiar"])
    # cabe na mesa (info.json)
    for r in linhas:
        p = glob.glob(os.path.join(SAIDA, r["etapa"], r["imagem"], f"{r['imagem']}_pf{int(round(r['limiar_n'] * 10)):02d}",
                                   r["imagem"], "info.json"))
        r["cabe"] = None
        if p:
            r["cabe"] = bool(json.load(open(p[0], encoding="utf-8"))["impressao"]["cabe_na_mesa"])
    return linhas


def main():
    destino = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SAIDA, "05_comparacoes", "planilha_viabilidade_v2.xlsx")
    dados = carregar()
    wb = Workbook()
    ws_rec, ws_cfg, ws_img, ws_pre = wb.active, wb.create_sheet("Configurações"), wb.create_sheet("Por imagem"), wb.create_sheet("Premissas")
    ws_rec.title = "Recomendação"

    # ------------------------------------------------------------- Por imagem
    cols = ["imagem", "etapa", "limiar", "fator", "espessura_mm", "n_placas", "pecas", "filamento_cm3", "fracao_glandular",
            "degrau_tecido_mm", "quant_medio_mm", "quant_rmse_mm", "quant_p95_mm", "cadeia_medio_mm"]
    tit = ["Imagem", "Etapa", "Limiar fibroso", "Fator de material", "Espessura da placa (mm)", "Placas", "Peças soltas",
           "Filamento (cm³)", "Fração glandular impressa", "Degrau (mm de tecido)", "Erro médio de quantização (mm tecido)",
           "RMSE (mm tecido)", "p95 (mm tecido)", "Erro da cadeia (mm tecido)", "Cabe na mesa 220×220"]
    cabecalho(ws_img, 1, tit, [9, 32, 9, 10, 11, 8, 9, 11, 12, 12, 16, 11, 11, 14, 11])
    for i, r in enumerate(sorted(dados, key=lambda x: (x["cod"], x["limiar_n"], x["imagem"])), 2):
        vals = [r[c] if c not in ("limiar",) else r["limiar_n"] for c in cols] + ["sim" if r["cabe"] else "não"]
        for j, v in enumerate(vals, 1):
            c = ws_img.cell(i, j, v)
            c.border = BORDA
    ws_img.freeze_panes = "A2"
    ws_img.auto_filter.ref = f"A1:{get_column_letter(len(tit))}{len(dados) + 1}"
    n_img = len(dados) + 1

    # ------------------------------------------------------------- Configurações (médias por fórmula)
    configs = sorted({(r["cod"], r["limiar_n"]) for r in dados})
    t2 = ["Etapa", "Limiar", "Espessura (mm)", "Fator", "Placas (média)", "Peças soltas (média)", "Filamento (cm³, média)",
          "Fração glandular (média)", "Degrau (mm tecido)", "Erro médio quant. (mm tecido)", "p95 (mm tecido)",
          "Imagens que cabem na mesa (de 8)", "Escore: erro", "Escore: p95", "Escore: peças", "Escore: filamento",
          "Escore: fração (perto de 13,7 %)", "ESCORE FINAL (0–100)", "Posição", "|fração − alvo| (auxiliar)", "Posição só pelo erro"]
    cabecalho(ws_cfg, 1, t2, [34, 8, 11, 8, 9, 10, 12, 12, 11, 14, 10, 13, 9, 9, 9, 10, 13, 12, 9, 12, 11])
    ref = lambda col: f"'Por imagem'!${col}$2:${col}${n_img}"  # noqa: E731
    n = len(configs)
    r0, r1 = 2, n + 1
    for i, (cod, lim) in enumerate(configs, 2):
        crit = f'{ref("B")},"{[r["etapa"] for r in dados if r["cod"] == cod][0]}",{ref("C")},{lim}'
        ws_cfg.cell(i, 1, ROTULO[cod])
        ws_cfg.cell(i, 2, lim)
        for j, col in ((3, "E"), (4, "D"), (5, "F"), (6, "G"), (7, "H"), (8, "I"), (9, "J"), (10, "K"), (11, "M")):
            ws_cfg.cell(i, j, f"=AVERAGEIFS({ref(col)},{crit})")
        ws_cfg.cell(i, 12, f'=COUNTIFS({crit},{ref("O")},"sim")')
        # escores 0–1 (min-max entre as configurações; 1 = melhor)
        ws_cfg.cell(i, 13, f"=(MAX($J${r0}:$J${r1})-J{i})/(MAX($J${r0}:$J${r1})-MIN($J${r0}:$J${r1}))")
        ws_cfg.cell(i, 14, f"=(MAX($K${r0}:$K${r1})-K{i})/(MAX($K${r0}:$K${r1})-MIN($K${r0}:$K${r1}))")
        ws_cfg.cell(i, 15, f"=(MAX($F${r0}:$F${r1})-F{i})/(MAX($F${r0}:$F${r1})-MIN($F${r0}:$F${r1}))")
        ws_cfg.cell(i, 16, f"=(MAX($G${r0}:$G${r1})-G{i})/(MAX($G${r0}:$G${r1})-MIN($G${r0}:$G${r1}))")
        ws_cfg.cell(i, 20, f"=ABS(H{i}-Recomendação!$C$12)")
        ws_cfg.cell(i, 20).number_format = "0.000"
        ws_cfg.cell(i, 21, f"=RANK(J{i},$J${r0}:$J${r1},1)")
        ws_cfg.cell(i, 17, f"=(MAX($T${r0}:$T${r1})-T{i})/(MAX($T${r0}:$T${r1})-MIN($T${r0}:$T${r1}))")
        ws_cfg.cell(i, 18, f"=100*(Recomendação!$C$5*M{i}+Recomendação!$C$6*N{i}+Recomendação!$C$7*O{i}+"
                           f"Recomendação!$C$8*P{i}+Recomendação!$C$9*Q{i})/SUM(Recomendação!$C$5:$C$9)")
        ws_cfg.cell(i, 19, f"=RANK(R{i},$R${r0}:$R${r1})")
        for j in range(1, 20):
            c = ws_cfg.cell(i, j)
            c.border = BORDA
            if j in (3, 4):
                c.number_format = "0.00" if j == 3 else "0.000"
            elif j in (5, 6, 7, 9, 10, 11):
                c.number_format = "0.0" if j in (5, 6, 7) else "0.00"
            elif j == 8:
                c.number_format = "0.0%"
            elif 13 <= j <= 17:
                c.number_format = "0.00"
            elif j == 18:
                c.number_format = "0.0"
                c.font = Font(bold=True)
    ws_cfg.freeze_panes = "C2"
    ws_cfg.conditional_formatting.add(f"R{r0}:R{r1}", ColorScaleRule(start_type="min", start_color="F8CBAD", end_type="max", end_color="A9D08E"))
    ws_cfg.conditional_formatting.add(f"J{r0}:J{r1}", ColorScaleRule(start_type="min", start_color="A9D08E", end_type="max", end_color="F8CBAD"))

    # ------------------------------------------------------------- Recomendação
    ws = ws_rec
    ws.sheet_view.showGridLines = False
    ws.column_dimensions["A"].width = 3
    ws.column_dimensions["B"].width = 44
    for col in "CDEFGH":
        ws.column_dimensions[col].width = 16
    ws["B1"] = "Phantom de mama em placas fibrosas — qual configuração é mais viável?"
    ws["B1"].font = Font(bold=True, size=15)
    ws["B2"] = ("Comparação de 7 etapas × 3 limiares fibrosos (média de 8 mamografias). Erro em mm de tecido equivalente. "
                "Altere os pesos (verde) e o ranking recalcula.")
    ws["B2"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("B2:H2")
    ws.row_dimensions[2].height = 32
    ws["B4"], ws["C4"] = "Pesos do escore (editáveis)", "Peso"
    for c in ("B4", "C4"):
        ws[c].font = Font(bold=True, color="FFFFFF")
        ws[c].fill = AZUL
    for i, (rot, p) in enumerate([("Erro médio de quantização", 40), ("Erro p95 (pior caso local)", 10),
                                  ("Peças soltas (manuseio/montagem)", 20), ("Filamento (custo/tempo)", 10),
                                  ("Fração glandular próxima da faixa da literatura", 20)], 5):
        ws.cell(i, 2, rot).border = BORDA
        c = ws.cell(i, 3, p)
        c.fill, c.border = VERDE, BORDA
    ws["B11"], ws["C11"] = "Alvo (referência)", "Valor"
    for c in ("B11", "C11"):
        ws[c].font = Font(bold=True, color="FFFFFF")
        ws[c].fill = AZUL
    ws["B12"], ws["C12"] = "Fração glandular alvo (limite inferior da literatura, 13,7–25,6 %)", ALVO_FG
    ws["C12"].number_format, ws["C12"].fill, ws["C12"].border = "0.0%", VERDE, BORDA
    ws["B12"].alignment = Alignment(wrap_text=True)
    ws.row_dimensions[12].height = 30

    ws["B14"] = "Ranking (as 5 melhores configurações)"
    ws["B14"].font = Font(bold=True, size=12)
    for j, t in enumerate(["Configuração", "Limiar", "Posição", "Escore (0–100)", "Erro médio (mm tecido)", "Placas",
                           "Peças soltas", "Filamento (cm³)"], 2):
        c = ws.cell(15, j, t)
        c.font, c.fill, c.border = Font(bold=True, color="FFFFFF"), AZUL, BORDA
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = 15 if j > 2 else 44
    ws.row_dimensions[15].height = 34
    for k in range(1, 6):
        i = 15 + k
        m = f"MATCH({k},'Configurações'!$S${r0}:$S${r1},0)"
        ws.cell(i, 2, f"=INDEX('Configurações'!$A${r0}:$A${r1},{m})")
        ws.cell(i, 3, f"=INDEX('Configurações'!$B${r0}:$B${r1},{m})")
        ws.cell(i, 4, k)
        ws.cell(i, 5, f"=INDEX('Configurações'!$R${r0}:$R${r1},{m})")
        ws.cell(i, 6, f"=INDEX('Configurações'!$J${r0}:$J${r1},{m})")
        ws.cell(i, 7, f"=INDEX('Configurações'!$E${r0}:$E${r1},{m})")
        ws.cell(i, 8, f"=INDEX('Configurações'!$F${r0}:$F${r1},{m})")
        ws.cell(i, 9, f"=INDEX('Configurações'!$G${r0}:$G${r1},{m})")
        for j in range(2, 10):
            ws.cell(i, j).border = BORDA
        ws.cell(i, 3).number_format = "0.0"
        ws.cell(i, 5).number_format = "0.0"
        ws.cell(i, 6).number_format = "0.00"
        ws.cell(i, 7).number_format = "0.0"
        ws.cell(i, 9).number_format = "0.0"
        if k == 1:
            for j in range(2, 10):
                ws.cell(i, j).fill, ws.cell(i, j).font = VERDE, Font(bold=True)
    # ---- conclusão (texto fixo, escrito a partir dos números medidos)
    ws["B21"] = "Melhor apenas pelo erro de quantização:"
    ws["B21"].font = Font(bold=True)
    mm = f"MATCH(1,'Configurações'!$U${r0}:$U${r1},0)"
    ws["C21"] = (f"=INDEX('Configurações'!$A${r0}:$A${r1},{mm})&\" · limiar \"&TEXT(INDEX('Configurações'!$B${r0}:$B${r1},{mm}),\"0.0\")"
                 f"&\" · erro \"&TEXT(INDEX('Configurações'!$J${r0}:$J${r1},{mm}),\"0.00\")&\" mm\"")
    ws.merge_cells("C21:I21")
    ws["B23"] = "Leitura dos resultados (valores medidos)"
    ws["B23"].font = Font(bold=True, size=12)
    notas = [
        "1) Pelo erro de quantização apenas: cera de abelha + limiar 1,0 + placas de 1,0 mm (erro médio 0,24 mm de tecido; 21 placas; 74 peças soltas). Com a parafina, o melhor erro é 0,29 mm (17 placas; 59 peças). Com os pesos padrão (que premiam também fração glandular próxima da literatura e poucas peças), lidera parafina + limiar 0,6 + 1,0 mm; veja o ranking acima.",
        "1b) A cera de abelha (fator 0,478) exige ~20 % mais PLA que a parafina (fator 0,387) e gera mais placas e peças; o μ da cera é calculado a partir de C46H92O2 e ainda não foi medido.",
        "2) Alternativa robusta: placas de 2,0 mm (10 camadas): erro 0,52 mm e só 30 peças soltas, contra 59 em 1,0 mm; porém pior que o 'antes' (0,36 mm).",
        "3) Placas de 3,4 mm com o fator de atenuação NÃO servem: sobram 5 placas e o erro é 0,81 mm (2,2× o 'antes').",
        "4) O fator de atenuação reduz o filamento em cerca de 60 % em qualquer espessura (23 contra 61 cm³ no limiar 1,0).",
        "5) Limiar 1,0 dá menor erro, mas fração glandular de 6,3–9,9 %, abaixo da faixa de 13,7–25,6 %; limiar 0,6 chega a 11–15 % com mais erro e mais peças.",
        "6) Vistas MLO (lmlo1, lmlo2, rmlo3, rmlo4) não cabem na mesa de 220 × 220 mm em nenhuma configuração.",
    ]
    for i, t in enumerate(notas, 24):
        ws.cell(i, 2, t).alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=i, start_column=2, end_row=i, end_column=9)
        ws.row_dimensions[i].height = 32
    ws["B31"] = ("Atenção: 'mais viável' depende dos pesos. O erro de quantização é só uma parte do erro total; "
                 "o fator de material (0,387) é calculado e ainda não validado por radiografia de um degrau impresso.")
    ws["B31"].alignment = Alignment(wrap_text=True, vertical="top")
    ws["B31"].fill = AMAR
    ws.merge_cells("B31:I31")
    ws.row_dimensions[31].height = 34

    # ------------------------------------------------------------- Premissas
    ws_pre.column_dimensions["A"].width = 120
    ws_pre["A1"] = "Premissas e limites"
    ws_pre["A1"].font = Font(bold=True, size=13)
    for i, t in enumerate([
        "Erro de quantização = |dq_sem_limpeza − d| / fator de material, na região da mama, em mm de TECIDO equivalente (mm de placa ÷ fator).",
        "Fator de material = (μ_glandular − μ_gordura) / (μ_PLA − μ_parafina) a 20 keV = 0,387 (src/atenuacao.py; validado contra NIST/ICRU-44, erro ≤ 0,03 %).",
        "Etapa 'antes' = fator 1,0 (sem atenuação). Etapas 'depois' = fator 0,387 (PLA 1,24 g/cm³, parafina 0,82 g/cm³, 20 keV monoenergético).",
        "Escore = soma ponderada de 5 critérios normalizados 0–1 entre as 21 configurações (1 = melhor); pesos editáveis na aba Recomendação.",
        "Fração glandular: referência 13,7–25,6 % (Yaffe 2009 via Bliznakova 2024); escore premia proximidade de 13,7 %.",
        "Pixel (0,07 / 0,1167 mm) e espessura da mama (50 mm) são suposições; as imagens provavelmente são 'for presentation'.",
        "Não incluídos: suavização, limpeza de manufatura, vazios de impressão, tempo de impressão real, fragilidade de placas finas (apenas indireta, via nº de peças).",
        "O G-code gerado pelo projeto nunca foi executado em impressora.",
        "Dados: saida/05_comparacoes/quantizacao_etapas/quantizacao.csv; geração: python src/etapas.py comparar; planilha: python src/planilha_viabilidade.py.",
    ], 3):
        ws_pre.cell(i, 1, "• " + t).alignment = Alignment(wrap_text=True, vertical="top")
        ws_pre.row_dimensions[i].height = 30

    wb.save(destino)
    print("planilha:", destino)

    # ---- espelho em Python do mesmo escore, para conferir o resultado das fórmulas
    por = defaultdict(list)
    for r in dados:
        por[(r["cod"], r["limiar_n"])].append(r)
    med = {k: {m: sum(x[m] for x in v) / len(v) for m in ("quant_medio_mm", "quant_p95_mm", "pecas", "filamento_cm3", "fracao_glandular")}
           for k, v in por.items()}

    def norm(f, inv=True):
        vs = [f(m) for m in med.values()]
        lo, hi = min(vs), max(vs)
        return {k: ((hi - f(m)) if inv else (f(m) - lo)) / (hi - lo) for k, m in med.items()}
    e, p95, pc, fl = (norm(lambda m: m["quant_medio_mm"]), norm(lambda m: m["quant_p95_mm"]), norm(lambda m: m["pecas"]),
                      norm(lambda m: m["filamento_cm3"]))
    fr = norm(lambda m: abs(m["fracao_glandular"] - ALVO_FG))
    esc = {k: 100 * (40 * e[k] + 10 * p95[k] + 20 * pc[k] + 10 * fl[k] + 20 * fr[k]) / 100 for k in med}
    print("\nescore (espelho em Python, pesos padrão):")
    for k, v in sorted(esc.items(), key=lambda x: -x[1])[:6]:
        print(f"  {ROTULO[k[0]]:<32} limiar {k[1]}: {v:5.1f}   erro {med[k]['quant_medio_mm']:.3f}  peças {med[k]['pecas']:.0f}")


if __name__ == "__main__":
    main()
