# Radiografar o material: para que serve e como fazer

## 1. A pergunta que a radiografia responde

O gerador decide **quantos milímetros de plástico** colocar em cada ponto para que a
peça "escureça" a radiografia como o tecido fibroglandular faria. Essa conta depende de
quanto o **plástico, do jeito que a sua impressora o deposita**, atenua o feixe de raios X
em comparação com o tecido.

**Atualização (orientação da pesquisa):** o fator passa a ser **calculado pela atenuação linear**
(seção 1b). Antes disso, `--fator-material` valia 1,0, ou seja, o código assumia que 1 mm de
plástico equivale a 1 mm de tecido, o que era errado. A radiografia continua necessária, mas
como **validação** do cálculo (e para capturar o que a tabela não vê: vazios, endurecimento do
feixe, densidade real da peça impressa).

## 1b. Cálculo pela atenuação linear

```
μ (cm⁻¹) = ρ (g/cm³) · (μ/ρ) (cm²/g)            (μ/ρ de mistura = Σ w_i · (μ/ρ)_i)

fator_material = (μ_glandular − μ_gordura) / (μ_placa − μ_parafina)      [mm de placa por mm de tecido]
```

Implementado em `src/atenuacao.py`; no gerador: `--energia-kev 20 --material-placa PLA`.
Os μ/ρ vêm da base Elam (pacote `xraydb`) combinada por frações de massa, e o código foi
conferido contra a tabela NIST XCOM/ICRU-44 (tecido adiposo e mama, 15–40 keV): erro ≤ 0,03 %
(`tests/testar_atenuacao.py`). Composições: ICRU-44 (NIST) para gordura e tecido mamário;
PLA = (C₃H₄O₂)ₙ, 1,24 g/cm³; parafina = C₂₅H₅₂, 0,82 g/cm³.

| Energia | Δμ tecido (cm⁻¹) | Δμ placa PLA (cm⁻¹) | fator |
|---|---|---|---|
| 20 keV | 0,163 | 0,423 | **0,387** |

O fator varre de ~0,40 (16 keV) a ~0,31 (34 keV) para PLA (`python src/atenuacao.py --varredura 16 34`).
Consequência prática: a espessura máxima de placa para 50 mm de mama cai de 50 para ~19 mm.

**Limites do cálculo:** é monoenergético (use a energia efetiva do seu espectro); a composição do
"fibroglandular" é a do ICRU-44 (outros autores usam 100 % glandular com composição diferente); PETG e
ABS são aproximados pela unidade de repetição; e a densidade usada é a do material maciço, então
vazios de impressão reduzem μ e exigem informar `--dens-placa`.

## 2. Por que tabelas não bastam

- A atenuação depende da **energia** do feixe (kVp, alvo e filtro). Plástico e tecido têm
  composição diferente, então variam de modo diferente com a energia: um fator medido a
  70 kVp não vale para 28 kVp.
- Uma impressora FDM não entrega a densidade do filamento: há **microfrestas** entre
  trilhas de extrusão e entre camadas, bolhas de ar (Bliznakova et al. 2024) e um padrão
  de camadas que aparece na radiografia (Clark et al. 2016). A atenuação real é a da peça
  impressa, não a do material na ficha técnica.
- Clark et al. (2016) caracterizaram o material exatamente assim: imprimiram **degraus**
  de espessura conhecida, radiografaram no mamógrafo e estimaram a atenuação média.

## 3. A equivalência correta (o plástico desloca parafina, não ar)

Dentro do recipiente, cada coluna de altura `t` é parafina com uma placa de espessura `d`
no lugar de parte dela. Na mama real, a mesma coluna é gordura com uma fração `f` de
tecido glandular. Igualando só a parte que **varia** de uma coluna para outra:

```
(μ_placa − μ_parafina) · d  =  (μ_glandular − μ_gordura) · f · t

d = f · t · (μ_glandular − μ_gordura) / (μ_placa − μ_parafina)
                \_____________ fator_material _____________/
```

Duas consequências:

1. O fator é uma razão de **contrastes**, não de atenuações. O que importa para a placa é o
   quanto ela atenua **a mais** que a parafina que ela substitui.
2. Como a diferença entre plástico e parafina é, em geral, bem maior que a diferença entre
   glândula e gordura, o fator tende a ser **menor que 1** (placas mais finas que as atuais).
   Calculado (seção 1b, 20 keV, PLA): Δμ do tecido = 0,163 cm⁻¹ e Δμ da placa = 0,423 cm⁻¹,
   fator 0,387. Só a medição confirma o valor para o filamento e a energia reais.

## 4. Procedimento

**O que imprimir**
- Uma pilha de **placas-teste** de 20 × 20 mm, com a mesma espessura, o mesmo filamento, o
  mesmo bico e os mesmos ajustes das placas reais (100 % de preenchimento, mesma altura de
  camada). Use de 1 a 14 placas empilhadas, como será a montagem real; assim as interfaces
  entre placas entram na medida.
- Se possível, uma versão **monolítica** (um degrau único de mesma altura) para comparar e
  ver o efeito das interfaces.

**Como radiografar**
1. No recipiente, com a **parafina até a altura final** (50 mm), coloque os degraus
   afastados entre si e das paredes. Mantenha também uma região só de parafina (referência).
2. Use **o mesmo equipamento e a mesma técnica** (kVp, alvo/filtro, compressão) que serão
   usados com o phantom e, de preferência, o modo de imagem **bruta** ("for processing"),
   em que o valor do pixel é proporcional à intensidade.
3. Repita 3 a 5 exposições e use a média, como fazem Bliznakova et al. (2024).

**Como analisar** (ImageJ ou Python)
1. Meça o valor médio do pixel em uma região de cada degrau (`I_k`) e na região só de
   parafina (`I_0`).
2. Calcule `Δμ_k = −ln(I_k / I_0) / t_k`, com `t_k` em cm. Se os valores forem
   praticamente constantes, `Δμ_placa` é a média. Se caírem com a espessura, é
   endurecimento do feixe: use a curva `espessura → −ln(I/I_0)` inteira (tabela de busca)
   no lugar de um número único.
3. Obtenha `μ_glandular − μ_gordura` no **mesmo espectro** (tabelas publicadas, por
   exemplo Johns & Yaffe, *Phys. Med. Biol.* 32:675, 1987, ou integrando dados do NIST sobre
   o espectro usado; esta fonte não foi lida aqui). Divida para obter `fator_material`.
4. Rode o gerador com `--fator-material <valor>`.

**Exemplo numérico (ILUSTRATIVO, não é medida)**

| Placas empilhadas | t (mm) | `I_k` | `−ln(I_k/I_0)/t` |
|---|---|---|---|
| 0 (só parafina) | 0 | 1000 | n/a |
| 1 | 3,2 | 790 | 0,74 cm⁻¹ |
| 2 | 6,4 | 622 | 0,74 cm⁻¹ |
| 3 | 9,6 | 490 | 0,74 cm⁻¹ |

Exemplo ilustrativo (números de tabela de medição, não nossos): com Δμ do tecido de 0,34 cm⁻¹ (suposto), o fator seria 0,34 / 0,74 ≈ **0,46**: 1 mm de
tecido glandular a mais que a gordura equivale a 0,46 mm de placa.

## 5. O que a medição também revela

- **Bolhas e microfrestas:** a radiografia da própria placa mostra vazios (pontos escuros)
  e padrão de camadas. Se a peça não for homogênea, mude o filamento, seque-o, ajuste o
  fluxo ou troque por resina.
- **Escala em mm por pixel:** incluir uma régua impressa ao lado dos degraus calibra o
  tamanho do pixel da radiografia.
- **Se a peça final reproduz a imagem:** compare a radiografia do phantom completo com a
  mamografia original (`src/comparar.py` e a macro do ImageJ servem para isso).

## 6. Sem acesso a um mamógrafo

A medida só vale no espectro em que foi feita. Alternativas, em ordem de confiança:
1. Mamógrafo de um hospital ou laboratório de física médica (a melhor opção).
2. Calcular o fator por tabelas de atenuação (NIST, composição do PLA e dos tecidos)
   integradas no espectro do mamografia desejado, assumindo plástico **sem vazios**. Dá um
   limite superior da atenuação real.
3. Usar o valor da resina de Clark et al. (0,93 cm⁻¹ a 28 kVp, Mo) apenas como ponto de
   partida. A resina deles não é o seu filamento.
