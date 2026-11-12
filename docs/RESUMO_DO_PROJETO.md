# Resumo do projeto: o que fizemos e por quê

Documento para estudo e anotações. Linguagem simples primeiro; detalhes técnicos depois.
Quando um número vem de simulação ou de suposição, isso está dito.

---

## 1. O objetivo em uma frase

Transformar uma **mamografia 2D** em peças que se **imprimem em 3D** (Ender 3), de modo que,
colocadas num recipiente com **parafina líquida** (que faz o papel da gordura), a peça se
pareça numa radiografia com a mamografia original. Só o **tecido fibroglandular** é impresso,
em **placas planas** empilháveis.

## 2. A ideia física (o que sustenta tudo)

- Um raio X que atravessa a mama perde intensidade conforme o material por onde passa
  (lei de Beer–Lambert): `I = I0 · exp(−∫μ dz)`. O detector registra só a **soma** da atenuação
  ao longo do raio. A profundidade `z` se perde.
- Por isso a mamografia **não permite reconstruir a anatomia 3D**. O que ela diz, em cada pixel,
  é **quanta espessura de material denso** o raio atravessou. A posição desse material em `z`
  não muda a imagem.
- Consequência usada no projeto: o relevo de cinza vira um **mapa de espessura** `d(x,y)`, e esse
  mapa é **fatiado em placas**. A placa `k` cobre onde `d ≥ (k − ½)·e` (`e` = espessura da placa).
  Empilhadas, as placas somam `d` em degraus de `e`.
- **Gordura = parafina líquida.** A placa substitui parafina, não ar. Por isso o que importa é o
  **contraste** (μ_placa − μ_parafina), e não a atenuação do plástico sozinha.

## 3. O pipeline passo a passo (`src/mama3d.py`)

1. **Leitura** (`src/dicom_io.py`): DICOM ou TIFF; tamanho do pixel em mm; inverte imagens em que
   o branco é o valor baixo; converte imagem bruta (RAW) em atenuação `ln(I0/I)`.
2. **Máscara da mama:** separa mama do fundo (limiar do método do triângulo), fica com o maior
   componente e preenche buracos.
3. **Orientação:** parede torácica sempre à esquerda (espelha se preciso).
4. **MLO** (vista oblíqua): remove o **músculo peitoral** e o **tronco** abaixo do sulco inframamário.
5. **Mapa fibroso** (`mapa_fibroso`): onde está a linha entre gordura e tecido fibroso? Um limiar
   (Otsu dentro da mama, multiplicado por `--limiar-fibroso`) vira o "zero"; o percentil 99,5 vira "100 %".
   Abaixo do zero é gordura e **não é impressa**.
6. **Faixas excluídas:** 5 mm junto à pele (a pele é a parede do recipiente) e 3 mm na borda da imagem
   (artefato de borda, que antes virava lascas fora da mama).
7. **Fatiamento em placas** (`--espessura-placa`, `--camada`): quantiza `d`, limpa traços finos demais
   e ilhas pequenas (`--traco-min`, `--area-min`).
8. **Recipiente:** bandeja com o contorno da mama (parede 3 mm, fundo 2 mm), para a parafina.
9. **Exportação:** um STL por placa, mais recipiente, guia de montagem em escala 1:1 (PDF), visualizador
   HTML e `info.json` com todas as medidas.
10. **G-code** (`src/gcode.py`): fatia cada placa para a Ender 3 no formato dos arquivos do pendrive.

## 4. O que fizemos, em ordem

| # | Etapa | Resultado |
|---|---|---|
| 1 | Discutimos as duas ideias (emaranhado de pixels × ImageJ 3D) | Com uma imagem 2D as duas dão um mapa de altura; escolhemos Python |
| 2 | Leitor de DICOM e TIFF com escala física | Funciona com 8/16 bits, RGB, RAW, pilhas |
| 3 | Modelos 3D (sólido, relevo, composto, volume DBT) e STL fechado, imprimível | Malhas válidas por *marching cubes*; redução de faces opcional |
| 4 | Clonamos e validamos o **mammoreplicator** (FDA) | A física dele está correta (R² 0,9998); a malha dele não é fechada |
| 5 | Você leu os artigos e mudou o foco: **só o fibroso, em placas** | Modo `placas` virou o padrão |
| 6 | Medimos a fidelidade: macro do **ImageJ** e script Python | Acerto de ~90 % (mede a conversão, não a física) |
| 7 | **Teste de perdas** por etapa (`src/perdas.py`) | Maiores perdas: suavização e limpeza de manufatura |
| 8 | Configuração para **Ender 3** (camada 0,2 mm, traço 0,8 mm) | Ganho de ~3 pontos de fidelidade |
| 9 | **Fusion** via MCP (porta 27183): abrir placas, medir, comparar | As medidas do Fusion bateram com as nossas |
| 10 | Corrigimos o **transbordamento** (faixa de borda da imagem) | Margem de 3 mm e verificação de que nada fica fora da mama |
| 11 | Lemos o **pendrive**: arquivos `.gcode` do Cura 5.13 (Marlin) | Criamos um gerador de G-code no mesmo formato |
| 12 | Estimativas de filamento, tempo e custo | Ver seção 6 |
| 13 | Imagens novas: `rmlo4` (MLO), `rcc`; comparações de espessura e de limiar | Ver seção 6 |
| 14 | Repositório **git** local e README com a teoria dos artigos | Um commit; mudanças posteriores ainda sem commit |

## 5. Parâmetros que você vai ouvir falar

| Parâmetro | O que faz | Efeito de aumentar |
|---|---|---|
| `--limiar-fibroso` | Multiplica o corte gordura × fibroso (padrão 1,0) | **Menos** tecido (corte mais alto); **baixar** inclui mais |
| `--espessura-placa` | Espessura de cada placa (mm) | Menos placas, mais degraus visíveis |
| `--max-placas` | Limite de placas (a espessura sobe se passar) | Mais placas permitidas |
| `--camada` | Altura de camada do fatiador (0,2 mm) | Força a espessura a ser múltiplo dela |
| `--traco-min` / `--area-min` | Menor traço e menor ilha impressáveis | Menos detalhe fino, menos peças soltas |
| `--suavizacao` | Suavização do relevo (mm) | Menos ruído, menos detalhe |
| `--margem-pele` / `--margem-borda` | Faixas sem placa junto à pele e à borda da imagem | Menos tecido na periferia |
| `--fator-material` | mm de placa por mm de tecido (hoje 1,0, **não calibrado**) | Ver `CALIBRACAO_RADIOGRAFICA.md` |
| `--fracao-alvo` | Reescala para uma fração glandular desejada | Fixa a fração, em vez do corte |
| `--pixel-mm`, `--espessura` | Tamanho do pixel e espessura da mama (suposições) | Mudam a escala de tudo |

**A variável que sua professora questionou** é o limiar de seleção (`zero`, em `src/mama3d.py`, função
`mapa_fibroso`). Com o padrão (1,0) só ~27 % da mama vira placa (fração glandular 8,8 % na `lcc1`).
Com `0,6` a fração vai a ~13,5 %, dentro da faixa da literatura (13,7–25,6 %).

## 6. Resultados principais

Todos são **simulações**. Nada foi impresso nem radiografado.

**`lcc1` (CC, 3328 × 4096, pixel 0,07 mm assumido)**

| Versão | Placas | Peças | Relevo preservado | Observação |
|---|---|---|---|---|
| Config. A (antes das melhorias) | 14 × 3,2 mm | 44 | 86,9 % | |
| Config. Ender 3 (traço 0,8, camada 0,2) | 14 × 3,4 mm | 91 | 90,6 % | 58,2 cm³; fração 8,8 % |
| Placas de 2,0 mm | 24 × 2,0 mm | 151 | 90,8 % | quase o mesmo plástico |
| Placas de 1,0 mm | 48 × 1,0 mm | 300 | não medido | |
| **Limiar 0,6** | 14 × 3,4 mm | 93 | 92,1 % | 89,5 cm³; fração 13,5 %; ~10,6 h; ~119 g |

**`rcc` (CC)**

| Versão | Placas | Relevo preservado | Tempo | Filamento |
|---|---|---|---|---|
| 13 placas de 3,4 mm | 13 | 88,2 % | ≈ 7,0 h | 74 g |
| 21 placas de 2,2 mm | 21 | 88,5 % | ≈ 7,9 h | 78 g |
| Limiar 0,6 | 14 × 3,4 mm | 91,1 % | ≈ 10,7 h | 120 g |

**`rmlo4` (MLO, 2457 × 1996, pixel 0,1167 mm inferido)**: 14 × 3,4 mm, 114 peças, 48,5 cm³,
fração 6,9 %, relevo preservado 82,3 %, ≈ 6,5 h, 67,5 g. O recipiente (119 × 234 mm) só cabe na
mesa girado na diagonal.

**Conclusões que se repetiram**
- **A espessura da placa quase não muda a fidelidade** (diferenças de décimos de ponto). Mais placas
  significam mais peças soltas, mais tempo (+12 %) e pouco mais plástico (+5 %).
- **O limiar de seleção muda muito mais** o resultado (quantidade de tecido, plástico e tempo).
- O "acerto pixel a pixel" **depende de como "denso" é definido**. Com o limiar 0,6 ele cai (≈ 80 %)
  por causa do critério da medida, e não por erro do modelo. Para esse caso vale o "relevo preservado".

## 7. Imprimir as placas juntas ou uma a uma?

**Dá para imprimir várias juntas, mas hoje o projeto gera um G-code por placa.**

- A mesa da Ender 3 tem 220 × 220 mm (484 cm²). A placa 1 mede cerca de 93 × 159 mm; as seguintes
  são menores (cada uma cabe dentro da anterior).
- **Encaixando só as caixas das placas** (com brim de 4 mm e folga de 2 mm), calculei:

| Conjunto | Pegada total | Camas de 220 × 220 mm |
|---|---|---|
| 14 placas de 3,4 mm | 171 cm² | **3** |
| 24 placas de 2,0 mm | 291 cm² | **4** |

  Ou seja, **24 placas não cabem numa mesa só**; caberiam em cerca de 4 impressões. Um encaixe mais fino
  (acompanhando o formato real das peças) poderia reduzir isso.
- **Tempo:** quase o mesmo no total, porque o tempo vem do volume extrudado. Juntar economiza o
  aquecimento e a troca de arquivos de cada placa.
- **Vantagens de juntar:** menos impressões para acompanhar e menos aquecimentos.
- **Desvantagens:** se uma impressão falhar, perdem-se várias placas; mais deslocamentos entre peças
  (fios de filamento); é preciso identificar qual placa é qual (elas se parecem).
- **Como fazer hoje:** no Cura ou Creality Print, importe vários STL juntos e use *Organizar tudo*
  (as placas de um mesmo conjunto têm a mesma altura, então as camadas coincidem). Juntar vários STL no
  nosso gerador de G-code ainda **não** está implementado.
- **Recomendação para o primeiro teste:** imprima **uma a uma**, começando pela placa 1, para validar
  o fluxo e a qualidade antes de juntar.

## 8. O que está validado e o que não está

**Validado (no computador)**
- STLs fechados, sem vértice fora da mama, dimensões batendo com o Fusion.
- G-code conferido lendo o arquivo (limites da mesa, Z, filamento, cobertura de 97–98 % da seção).
- Equivalência entre a macro do ImageJ e o Python (diferença < 0,5 ponto).

**Não validado**
- Nenhuma peça foi impressa; o G-code é de gerador próprio e **nunca rodou numa impressora**.
- Nenhuma radiografia da peça física.
- **Suposições:** pixel de 0,07 mm (e 0,1167 mm nas imagens menores), espessura da mama de 50 mm.
- **μ do plástico não medido** (`--fator-material` = 1,0); a fração glandular vem de limiar, não de física.
- As imagens são 8 bits processadas (algumas "C-View"); a literatura recomenda imagem bruta.

## 9. Próximos passos sugeridos

1. Confirmar com a professora o limiar de seleção (0,6 ou outro) e a fração alvo.
2. Imprimir a placa 1 (Ender 3) e conferir medidas, bolhas e aderência.
3. Imprimir um degrau do material e radiografar, para medir μ e calibrar `--fator-material`.
4. Obter imagem bruta (DICOM "for processing") e as informações reais de pixel e espessura.
5. Focalizar as placas para o foco do tubo (hoje empilham na vertical).
6. Decidir se vale implementar o G-code de várias placas numa mesma mesa.

## 10. Glossário rápido

- **Fibroglandular:** tecido denso da mama (glândulas e estroma); aparece claro na mamografia.
- **CC / MLO:** vistas da mamografia (craniocaudal e oblíqua mediolateral).
- **RAW / "for processing":** imagem bruta, proporcional à intensidade. "For presentation": já processada.
- **Beer–Lambert:** lei que relaciona a intensidade transmitida à atenuação do material.
- **Contraste da placa:** quanto a placa atenua a mais que a parafina que ela substitui.
- **Brim:** faixa extra de filamento na primeira camada para a peça aderir à mesa.
- **Otsu / triângulo:** métodos automáticos para escolher um limiar numa imagem.
- ***Marching cubes*:** algoritmo que gera a superfície 3D a partir de um volume.

## 11. Onde está cada coisa

```
README.md                       visão geral, teoria dos artigos, viabilidade
docs/USO.md                     comandos e parâmetros
docs/ENDER3.md                  fatiamento e impressão na Ender 3
docs/CALIBRACAO_RADIOGRAFICA.md como radiografar o material e calcular o fator
docs/RESUMO_DO_PROJETO.md       este arquivo
src/mama3d.py                   pipeline principal
src/gcode.py                    gerador de G-code
src/comparar.py, perdas.py      medidas de fidelidade e de perdas
saida/                          resultados (não versionados)
dados/                          imagens (não versionadas)
```
