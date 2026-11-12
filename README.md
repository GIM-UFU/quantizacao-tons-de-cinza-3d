# Phantom de mama a partir de mamografia 2D: placas fibrosas

Pipeline em Python que transforma uma mamografia (DICOM ou TIFF) em peças imprimíveis
em 3D: **só o tecido fibroglandular**, fatiado em **placas planas** impressas uma a uma, e
um **recipiente** com o contorno da mama. O tecido adiposo é substituído por um líquido
(parafina líquida) colocado no recipiente.

> **Status:** protótipo de pesquisa. Testado com uma imagem real (`lcc1.tif`) e com
> imagens sintéticas. **Nada foi impresso nem radiografado ainda.** Veja a
> [seção 4](#4-a-primeira-impressão-é-viável) antes de imprimir.

Guia de comandos e parâmetros: [docs/USO.md](docs/USO.md).

---

## 1. Fundamentação teórica

Esta seção resume o que cada artigo da pasta `Literatura/` afirma e o que o código usa.
Onde o código faz algo que **não** vem dos artigos, isso é dito na [seção 3](#3-o-que-é-nosso-e-não-vem-dos-artigos).

### 1.1 A mamografia é uma projeção (Beer–Lambert)

O detector registra a intensidade transmitida, que depende da atenuação **somada ao longo
do raio** (Bliznakova et al. 2024, Eq. 1):

$$I = I_0 \, e^{-\int_0^{l} \mu(E,x)\,dx}$$

Consequências para este projeto:

- A profundidade `z` é integrada e some. Reconstruir a anatomia 3D a partir de **uma**
  imagem é um problema inverso mal posto; Clark et al. (2016) ressalvam que o método não
  se destina a criar mamas 3D anatomicamente corretas.
- Para a imagem 2D importa só a **espessura total de material atenuante** que cada raio
  atravessa. A posição em `z` dessa espessura não altera a projeção (ponto de partida do
  fatiamento em placas, seção 3).
- Isso vale para **um** espectro e **uma** geometria de aquisição. O objeto só reproduz a
  imagem se for radiografado com o espectro e a geometria usados no projeto.

### 1.2 Do tom de cinza à espessura de material (mammoreplicator)

Clark, Ghammraoui & Badal (SPIE 9783, 2016) e Badal, Clark & Ghammraoui (*J. Med.
Imaging* 5:033501, 2018) descrevem o método do `mammoreplicator` (código da FDA, em
`externo/`, arquivado):

- Cada pixel é transformado em uma **coluna de material** cuja altura produz a mesma
  atenuação que a observada. Como qualquer material de atenuação conhecida serve, não é
  preciso material equivalente a tecido para reproduzir o feixe primário.
- As colunas são **inclinadas para o foco do tubo** (SOD de 66 cm nos testes), para que cada
  coluna se projete em um único pixel e não vire um traço.
- Conversão usada no código (`mammoreplicator.py`), com `T` a espessura comprimida,
  `I_med` a intensidade média da mama, `I_ar` a do ar e `mfp` o livre caminho médio
  (1,284 cm para a mama a 16,5 keV):

  `h(pixel) = T / ln(I_med / I_ar) · (mfp_plástico / mfp_mama) · ln(pixel / I_ar)`

- **Exige imagem RAW ("for processing")**. A imagem "for presentation" não serve porque
  inclui escala logarítmica, achatamento do contorno e outros processamentos.
- **Limites declarados:** o espalhamento da peça não representa o do órgão real; ruído
  quântico, borramento do detector e grade antiespalhamento ficam embutidos na imagem de
  origem; a peça só vale para o espectro usado no projeto.
- A ideia vem de Caldwell & Yaffe (1990), o "Rachel phantom": bloco esculpido de material
  equivalente a tecido (baixas frequências) mais filme dopado com mercúrio (altas).

### 1.3 O que as impressoras conseguem reproduzir (Clark et al. 2016)

| Tecnologia | Impressora | Resultado relatado |
|---|---|---|
| Extrusão (FDM), ABS | Stratasys uPrint SE | Aparece um **padrão quadriculado** das camadas de extrusão na radiografia; os autores a consideram inadequada |
| Estereolitografia (SLA) | Formlabs Form1+ | Preserva detalhe fino; a peça testada saiu ruim por problemas de impressão |
| Jato de tinta (PolyJet) | Stratasys Objet260 | Objetos de tamanho real e detalhados; cerca de 1 kg e ~10 h por peça |

- **Resolução:** nenhum detalhe abaixo de ~300 µm foi reproduzido com fidelidade; acima
  disso, microcalcificações e textura glandular saíram bem. O melhor agrupamento de
  pixels foi **3×3** (colunas de ~0,3 mm para pixels de 0,1 mm).
- **Material (SLA):** atenuação média de **0,93 cm⁻¹** para um feixe de 28 kVp (Mo), semelhante
  ao PMMA, porém ~19 % menor entre 10 e 20 keV. Cores diferentes têm perfis muito parecidos.
- **Caracterização do material:** feita com **degraus impressos** (espessuras conhecidas)
  radiografados no sistema, e por espectroscopia com detector de CdTe a 25 kVp.

### 1.4 Materiais e montagem física (Bliznakova et al. 2024)

Estudo de 7 materiais impressos (PLA, ABS, PET-G; resinas White v4, Flexible 80A, Model
v2, Wax40) como lesões, radiografados a 28 kVp (W/Rh) dentro de um recipiente com **água**
(mama densa) ou **parafina líquida** (mama gordurosa).

| Item | Valor no artigo |
|---|---|
| Forma da mama | semicilindro de **50 mm** de espessura comprimida |
| Pele | **5 mm** (recipiente de PLA de 3 mm + envelope de 2 mm) |
| Densidades (g/cm³) | adiposo 0,93 · glandular 1,04 · 50/50 0,98 · pele 1,09 |
| Materiais impressos (g/cm³) | PLA 1,24 · ABS 1,04 · PET-G 1,23 · resinas 1,15–1,20 após cura |
| Preenchimentos (g/cm³) | água 1,00 · parafina líquida 0,82 |

Achados que importam aqui:

- PLA e White v4 foram os mais densos; ABS e Flexible 80A, os menos densos. Em fundo de
  água, a lesão de PLA ficou quase invisível (atenuação próxima da água); em parafina, a de
  PLA teve o maior contraste.
- **Bolhas de ar** foram o problema recorrente: PLA (várias tentativas, mesmo após 1 h de
  imersão), PET-G (porosidade e umidade) e Model v2. White v4 e Wax40 tiveram menos
  vazios e menos ruído. ABS flutuou na água e deslocou a peça entre projeções.
- Conclusão dos autores: as **resinas** (Flexible 80A e White v4) foram o melhor material
  para lesões de baixo e alto contraste, considerando as dificuldades práticas.
- A posição das lesões no recipiente foi difícil de repetir entre exposições.

### 1.5 Composição da mama (Bliznakova 2024 e AAPM TG-234)

- A fração fibroglandular média (com pele) varia de **13,7 % a 25,6 %** entre grupos de
  mama (Yaffe et al., *Med. Phys.* 36:5437, 2009, citado por Bliznakova). Mama 100 %
  glandular é irreal; modelos experimentais usam até 75 %.
- **AAPM TG-234** (Bakic et al., doi:10.1002/mp.70215) trata de phantoms
  **virtuais**; o que aproveitamos dele:
  - os phantoms UMass convertem valores reconstruídos em **percentual de tecido glandular
    por voxel**, análogo 3D da nossa "fração fibroglandular × espessura";
  - os phantoms Duke foram fabricados em versões com 1 e 2 materiais; a de 1 material cobriu
    mais faixa dinâmica (64 % de densidade glandular) mas teve **mais bolhas de ar**;
  - o phantom VICTRE já foi impresso por jato de tinta sobre papel-manteiga (Ikejimba et al.).

---

## 2. Do artigo para o código

| O que o artigo diz | Onde está no código | Observação |
|---|---|---|
| Beer–Lambert; imagem RAW exige `ln(I0/I)` | `linearizar()` em `src/mama3d.py` | detecta RAW pela borda da imagem e converte |
| "For presentation" não serve (Clark 2016) | aviso e limitações (seção 6) | a imagem `lcc1.tif` parece ser desse tipo |
| Espessura comprimida de 50 mm (Bliznakova) | `--espessura` (padrão 50 mm no TIFF) | DICOM traz `BodyPartThickness` |
| Fração fibroglandular 13,7–25,6 % | `--fracao-alvo` | opcional; o padrão não força nenhum valor |
| Parafina líquida como gordura | `recipiente.stl` e `parafina_estimada_mL` em `info.json` | |
| Pele de 5 mm no recipiente | `--parede` (padrão 3 mm) e `--margem-pele` (padrão 5 mm) | a parede padrão é menor que a do artigo |
| Menor detalhe confiável ~0,3 mm (SLA/PolyJet) | `--traco-min` (padrão 1,2 mm) | piso mais conservador, pensando em FDM |
| Material do degrau e μ do material (Clark) | `--fator-material` | **ainda não calibrado** (seção 4) |
| Resinas como melhor material para lesões (Bliznakova) | sem parâmetro | escolha do material é feita na impressão |
| Quantidade de tecido glandular por voxel (TG-234, UMass) | `mapa_fibroso()` | análogo 2D |

## 3. O que é nosso e não vem dos artigos

Estas decisões são **nossas**; nenhum dos artigos as valida:

1. **Fatiar em placas planas.** Os artigos imprimem colunas de altura variável
   (Clark) ou lesões inteiras (Bliznakova). Aqui o mapa de espessura densa `d(x,y)` é
   quantizado em `N` placas de espessura `e`: a placa `k` cobre a região onde
   `d ≥ (k − ½)·e`. Empilhadas, somam `d` em degraus de `e`. Parte da premissa da seção 1.1.
2. **Separar gordura × fibroso** pelo limiar de Otsu dentro da mama. Gordura vira zero
   (não é impressa); o percentil 99,5 vira 100 %. É heurístico.
3. **Perfil de espessura** semicircular na borda (`mapa_espessura`): a mama comprimida entre
   duas placas tem espessura `H` no centro e cai na borda.
4. **Segmentação** da mama pelo método do triângulo (Zack et al., 1977), porque o Otsu
   cortava a gordura em imagens reais.
5. **Malhas** por *marching cubes* (Lorensen & Cline, 1987) e redução opcional de faces por
   quádricas (Garland & Heckbert, 1997).

## 4. A primeira impressão é viável?

**Resposta curta:** **sim, como impressão de teste de engenharia** (geometria, manuseio,
montagem no recipiente, comportamento do material); **ainda não, como phantom
quantitativo** que reproduza a mamografia.

### O que a literatura sustenta

- Imprimir só a estrutura densa e preencher o resto com líquido é a montagem de
  Bliznakova et al. (2024); parafina líquida (0,82 g/cm³) e PLA/resina são materiais
  testados por eles.
- Em 2D, só a espessura total de material importa (Clark 2016), então empilhar placas é
  coerente com a física, desde que o material seja calibrado.
- O tamanho das nossas estruturas (traço mínimo 1,2 mm; placas de 3,2 mm) está acima do
  piso de ~0,3 mm das impressoras testadas.

### O que impede de chamar de fiel (e o que fazer)

| Risco | Evidência | Para a 1ª impressão |
|---|---|---|
| **A imagem provavelmente não é RAW.** `lcc1.tif` é um TIFF de 8 bits, sem metadados; se for "for presentation", o cinza não é proporcional à atenuação | Clark 2016: "for presentation" não é adequada | tratar como protótipo **morfológico**; para fidelidade, usar DICOM RAW |
| **Escala desconhecida.** Pixel de 0,07 mm e espessura de 50 mm são **suposições** | TIFF não traz calibração | confirmar os dois valores antes de imprimir |
| **μ do material não medido.** `--fator-material` está em 1,0 | Clark: degraus impressos + radiografia | imprimir um **degrau** do material escolhido e medir |
| **Bolhas de ar** em PLA e PET-G; padrão quadriculado em FDM | Bliznakova 2024; Clark 2016; TG-234 (Duke) | preferir **resina** (White v4); em FDM, 100 % de preenchimento, secar o filamento e degaseificar sob a parafina |
| **Placas afundam e deslizam** na parafina (PLA 1,24 e ABS 1,04 contra 0,82 g/cm³) | densidades de Bliznakova; dificuldade de posicionar peças | usar o `guia_montagem_1para1.pdf`; considerar guias de registro no recipiente |
| **Colunas não focalizadas.** As placas empilham na vertical, sem inclinação para o foco | Clark 2016: colunas focalizadas (SOD 66 cm) | estimativa geométrica: deslocamento de até ~8 mm na borda das placas e ~2 mm na região densa central (estimativa nossa, SOD 660 mm). Correção possível: encolher cada placa no plano em `(1 − h/SOD)` |
| **Degraus de contraste.** Cada placa de 3,2 mm de resina (μ ≈ 0,93 cm⁻¹) reduz a intensidade em ~26 % | cálculo com μ de Clark (a resina de Clark não é a nossa) | reduzir `--espessura-placa` (mais placas) se o banding atrapalhar |
| **Fração glandular baixa**: 8,8 % no `lcc1`, abaixo da média da literatura | Yaffe 2009 | pode ser mama gordurosa; conferir com `--fracao-alvo` |
| **Espalhamento e espectro** não equivalentes ao tecido | Clark 2016 | aceitável para o teste; imagem só no espectro previsto |

### Plano sugerido (do mais barato ao mais caro)

1. **Antes de imprimir (horas):** confirmar tamanho do pixel e espessura da mama; decidir
   material (resina White v4 recomendada pela literatura) e a impressora (a maior placa mede
   93 × 159 mm; em uma Form 3, de 145 mm, só cabe na diagonal, a conferir).
2. **Degrau de calibração** no material escolhido, radiografado no mamógrafo (se houver acesso).
   Dá o μ real e ajusta `--fator-material`.
3. **Teste parcial:** imprimir 2 ou 3 placas grandes (01 a 03) e o recipiente. Verificar encaixe,
   bolhas (radiografia ou transiluminação) e o comportamento na parafina.
4. **Montagem completa** com as 14 placas, só depois de 2 e 3.

**Critério para seguir:** sem bolhas visíveis nas placas, encaixe sem folga excessiva no
recipiente e μ medido dentro do esperado.

---

## 5. Resultados medidos até agora (`lcc1.tif`)

Medidas do pipeline, **não** da peça física:

| Item | Valor |
|---|---|
| Placas | 14 × 3,2 mm (pilha de 44,7 mm), 44 peças soltas, 58,6 cm³ |
| Recipiente | 109 × 204 × 53 mm; parafina ≈ 706 mL |
| Fração glandular impressa | 8,8 % (heurística) |
| Acerto pixel a pixel, modelo × imagem | 90,5 % na mama (ImageJ) e 90,6 % (Python) |
| Variação do relevo preservada no STL | 86,9 % (volume retido 90,6 %) |
| Maiores perdas na cadeia | suavização 5,5 % e limpeza de manufatura 4,4 % |
| Configuração para Ender 3 ([docs/ENDER3.md](docs/ENDER3.md)) | variação preservada 90,6 %, acerto 91,8 % (Python), 14 placas × 3,4 mm, 91 peças; com margem de 3 mm na borda da imagem (artefato de borda removido) |

O acerto usa o **mesmo critério** (limiar de Otsu) que gera as placas; mede a perda na
conversão, não a fidelidade física.

## 6. Limitações e pendências

- Calibrar o cinza para fração glandular (hoje heurística por Otsu) e medir o μ do material.
- Usar imagem **RAW**; para RAW, a conversão física da fração glandular ainda não foi
  implementada (o teste sintético superestima: 50,6 % contra 28,7 % reais).
- Focalizar as placas para o foco do tubo.
- Nenhuma validação por radiografia da peça física.
- Escalonamento do protótipo existe (`--escala-impressao`), mas foi **deixado de lado** por decisão do projeto.

## 7. Estrutura do repositório

```
src/          mama3d.py (pipeline) · dicom_io.py (leitura) · comparar.py + comparar_imagej.ijm · perdas.py · inspecionar.py
tests/        geradores de imagens sintéticas e scripts de validação do mammoreplicator
docs/USO.md   comandos, parâmetros e modos
externo/      código de terceiros (não versionado): veja externo/LEIAME.md
Literatura/   artigos (PDFs não versionados): veja Literatura/LEIAME.md
dados/        imagens (não versionadas)
```

## 8. Referências

1. Clark M, Ghammraoui B, Badal A. *Reproducing 2D breast mammography images with 3D printed phantoms.* Proc. SPIE 9783, 97830B (2016). doi:10.1117/12.2217215
2. Badal A, Clark M, Ghammraoui B. *Reproducing two-dimensional mammograms with three-dimensional printed phantoms.* J Med Imaging 5(3):033501 (2018). doi:10.1117/1.JMI.5.3.033501. Código: <https://github.com/DIDSR/mammoreplicator> (arquivado).
3. Bliznakova K, Dukov N, Toshkova-Velikova O, Bliznakov Z. *Assessment of a method for manufacturing realistic breast lesions for experimental investigations.* Front Phys 12:1441740 (2024). doi:10.3389/fphy.2024.1441740 (CC BY)
4. Bakic PR et al. *AAPM task group 234 report: Virtual tools for the evaluation of new 3D/4D breast imaging systems.* Med Phys, doi:10.1002/mp.70215 (aceito em 04/11/2025).
5. Caldwell CB, Yaffe MJ. *Development of an anthropomorphic breast phantom.* Med Phys 17:273–280 (1990). doi:10.1118/1.596506
6. Yaffe MJ et al. *The myth of the 50-50 breast.* Med Phys 36:5437 (2009). doi:10.1118/1.3250863 (citado em 3)
7. Ivanov D et al. *Suitability of low density materials for 3D printing of physical breast phantoms.* Phys Med Biol 63:175020 (2018). doi:10.1088/1361-6560/aad315 (citado em 3)
8. Lorensen WE, Cline HE. Marching cubes. SIGGRAPH (1987). · Garland M, Heckbert P. Surface simplification using quadric error metrics. SIGGRAPH (1997). · Otsu N. IEEE Trans SMC 9:62 (1979). · Zack GW et al. J Histochem Cytochem 25:741 (1977).

As referências 5 a 7 foram lidas apenas por citação nos artigos 1 a 4, não no texto original.
