# Roteiro para a apresentação do grupo (imagens médicas)

Sugestão de 10 a 12 slides. Todos os números abaixo foram medidos neste repositório; o que
é suposição está marcado.

| # | Slide | Conteúdo | De onde tirar |
|---|---|---|---|
| 1 | Problema | Phantom de mama a partir de mamografias reais; imprimir só o tecido fibroglandular, a gordura vira parafina | README seção 1 |
| 2 | Física | Mamografia = projeção (Beer–Lambert): só a espessura atenuante total importa; 2D→3D é mal posto | README 1.1 |
| 3 | Estado da arte | mammoreplicator (Clark/Badal), Bliznakova 2024, AAPM TG-234 | README 1.2–1.5 e tabela da seção 2 |
| 4 | Nossa ideia | Fatiar o mapa de espessura em **placas** de espessura `e`; placa `k` cobre `d ≥ (k−½)e`; pilha no recipiente com parafina | `saida/mostrador_pf/` (visualizador) |
| 5 | **Fator de material (novo)** | `fator = (μ_gland − μ_gordura)/(μ_placa − μ_parafina)` por atenuação linear; 20 keV, PLA → **0,387** | `python src/atenuacao.py --varredura 16 34` |
| 6 | Validação do cálculo | Reproduz a tabela NIST/ICRU-44 (adiposo e mama, 15–40 keV) com erro ≤ 0,03 % | `python tests/testar_atenuacao.py` |
| 7 | Efeito no modelo (rcc2) | Fator 1,0 → 0,387: pilha de 44,2 → 18,0 mm, filamento de 60,4 → 23,1 cm³ (−62 %), fração glandular igual (9,7–9,8 %) | comparar `info.json` |
| 8 | Sensibilidade ao limiar | Limiar 0,6 / 0,8 / 1,0: fração glandular 12–15 % / 9–12 % / 7–10 %; só o 0,6 se aproxima da faixa de 13,7–25,6 % da literatura | `saida/mostrador_pf/index.html` |
| 9 | Espessura e nº de placas | 24→15 placas altera pouco o contraste (SSIM igual); o erro dominante é suavização/limpeza | medido na rcc2 (comparação 24×2,0 / 15×3,2 / 14×3,4 mm); ainda não está no README |
| 10 | Reprodutibilidade | Cada modelo guarda `info.json` com parâmetros, fator de material e `atenuacao` (μ de cada material); os comandos estão em `docs/USO.md` | `info.json` das pastas de saída |
| 11 | Limitações (honestidade) | Cinza→espessura é heurístico (imagens provavelmente "for presentation"); pixel e 50 mm de mama são suposições; cálculo monoenergético; PLA tem vazios; MLO não cabe na mesa; G-code nunca rodou numa impressora | README seção 4 e 6 |
| 12 | Próximos passos | Imprimir **degrau** em PLA e radiografar em parafina para validar o 0,387; DICOM RAW; placas focalizadas; polienergético | `docs/CALIBRACAO_RADIOGRAFICA.md` |

## Perguntas que a banca pode fazer

- **"De onde vem o 0,387?"** Δμ do tecido = 0,163 cm⁻¹ (mama ICRU-44 − adiposo ICRU-44) dividido por
  Δμ da placa = 0,423 cm⁻¹ (PLA 1,24 g/cm³ − parafina 0,82 g/cm³), a 20 keV.
- **"E o espectro polienergético?"** O cálculo é monoenergético na energia efetiva; a varredura
  15–40 keV mostra a sensibilidade (0,40 a 0,31 para PLA entre 16 e 34 keV). Um espectro real exigiria
  integrar a atenuação ponderada pelo espectro (limitação declarada).
- **"A densidade do PLA impresso é 1,24?"** Não necessariamente: vazios a reduzem. Há `--dens-placa`
  para informar o valor medido.
- **"Por que não 1:1?"** Porque o PLA atenua muito mais que a parafina (Δμ maior) do que a glândula
  atenua a mais que a gordura; logo menos espessura de plástico reproduz o mesmo contraste.
