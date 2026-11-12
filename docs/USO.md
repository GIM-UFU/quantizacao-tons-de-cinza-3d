# Guia de uso

> Fundamentação teórica, viabilidade e visão geral: veja o [README](../README.md).

## Instalação (uma vez)
```bash
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

## Uso
Coloque os DICOMs e/ou TIFFs em `dados/`. O formato é detectado pelo conteúdo,
não pela extensão.

### TIFF
TIFF não carrega os metadados clínicos do DICOM. Duas informações exigem atenção:
- **Tamanho do pixel:** é lido dos metadados do ImageJ ou das tags de resolução
  (dpi ou px/cm). Os valores 72 e 96 dpi são ignorados, porque são padrão de tela e
  não medida física. Sem calibração, o script assume 0,1 mm e avisa; nesse caso,
  passe `--pixel-mm`.
- **Espessura da mama comprimida:** não existe em TIFF, então o padrão é 50 mm.
  Passe `--espessura` se souber o valor real.

Também são tratados: 8 ou 16 bits, RGB (convertido para cinza), MinIsWhite (invertido),
compressão LZW/Deflate e pilhas multi-página (convertidas em volume, como a DBT).

```bash
.venv\Scripts\python src\mama3d.py dados --pixel-mm 0.07 --espessura 48
```

```bash
.venv\Scripts\python src\inspecionar.py dados
.venv\Scripts\python src\mama3d.py dados --modo solido
.venv\Scripts\python src\mama3d.py dados --modo relevo
```

Saída em `saida/<nome>/`: `impressao.stl` (pronto para o fatiador), `modelo.glb` (com cor),
`visualizar.html` (abrir no navegador, precisa de internet para o three.js),
`painel.png` (etapas 2D) e `info.json`.

## Modo placas (padrão): só o tecido fibroglandular
Segue Bliznakova et al., *Front. Phys.* 12:1441740 (2024): imprime-se apenas o tecido
denso (PLA ou resina, 100 % de infill). O tecido adiposo é parafina líquida dentro de
um recipiente impresso com o contorno da mama, cuja parede faz o papel da pele.

O mapa de espessura densa `d(x,y)` é fatiado em placas planas: a placa `k` cobre a
região onde `d ≥ (k−½)·e`. Em mamografia 2D só a espessura total atravessada pelo
raio importa (Beer–Lambert), então a pilha de placas reproduz `d` em passos de `e`.

```bash
.venv\Scripts\python src\mama3d.py dados\lcc1.tif --voxel 0.5 --pixel-mm 0.07
```

Saídas em `saida/<nome>/`:
- `placas/placa_01.stl … placa_NN.stl`: imprimir deitadas, uma a uma (01 fica embaixo).
- `recipiente.stl`: bandeja com o contorno da mama (parede 3 mm, fundo 2 mm).
- `guia_montagem_1para1.pdf`: uma página por placa, em escala real, para posicionar as peças.
- `visualizar.html`, `painel.png`, `info.json` (peças por placa, volumes, parafina estimada).

Parâmetros principais:
- `--espessura-placa 2` e `--max-placas 15`: se passar do limite, a espessura aumenta.
- `--fracao-alvo 0.15`: reescala para uma fração glandular volumétrica (literatura:
  13,7–25,6 % em média, Yaffe et al., 2009).
- `--limiar-fibroso 0.8`: inclui feixes fibrosos mais tênues (`>1` deixa só o mais denso).
- `--base-placa 0.6`: lâmina que une as ilhas de cada placa numa peça só.
- `--fator-material`: mm de material por mm de tecido glandular (μ_glandular/μ_material).
- `--margem-pele`, `--traco-min`, `--area-min`, `--parede`, `--fundo`, `--folga`.

**Limitações:** em imagem "FOR PRESENTATION" a escala de cinza não é física; o corte
gordura × fibroso é o Otsu dentro da mama. Em imagem RAW a fração sai superestimada
(falta a conversão física `f = (A/t − μ_a)/(μ_g − μ_a)`).

## Protótipo escalonado (opcional; fora do fluxo principal)
`--escala-impressao 0.5` gera as placas e o recipiente em metade do tamanho (1/8 do
volume). As restrições de manufatura (`--traco-min`, `--area-min`, `--parede`, `--fundo`,
`--folga`) continuam valendo no tamanho IMPRESSO; por isso, ao reduzir, o gerador
descarta mais detalhe fino no modelo. O guia de montagem também sai na escala impressa.

```bash
.venv\Scripts\python src\mama3d.py dados\lcc1.tif --voxel 0.5 --pixel-mm 0.07 --escala-impressao 0.5 --saida saida\placas_esc50
```

## Quanto se perde do relevo de cinza (por etapa)
```bash
.venv\Scripts\python src\perdas.py saida\placas\lcc1
```
Compara o relevo "ideal" (o que o cinza pede, pixel a pixel) com cada etapa até o STL:
faixa da pele, suavização, quantização em placas, limpeza de manufatura e malha final.
Mede volume retido e variação preservada (R²). Grava `resultado_perdas.json` e `perdas.png`.

## Comparar o modelo com a imagem (ImageJ)
A mamografia é uma projeção, então o modelo é projetado de volta em 2D (soma das
espessuras das placas) e comparado com a imagem original.

```bash
.venv\Scripts\python src\comparar.py saida\placas\lcc1
```

Gera em `saida/placas/lcc1/comparacao/`:
- `projecao_modelo.tif`: mm de material denso por pixel, extraído do `visualizar.html`
  (`--stl` usa os STL das placas).
- `resultado_python.json` e `mapa_acerto_python.png`: métricas de referência.
- `comando_imagej.txt`: linha de comando para rodar `src/comparar_imagej.ijm` em lote.

O ImageJ 1.54 (com Java próprio) está em `externo/imagej/ImageJ`; a linha do
`comando_imagej.txt` roda a macro em lote com ele.

No ImageJ/Fiji, de forma interativa: *Plugins → Macros → Run…*, escolher `src/comparar_imagej.ijm`, depois a
imagem original e o `projecao_modelo.tif`. A macro grava `resultado_imagej.csv` e
`mapa_acerto_imagej.png` (verde = acerto, vermelho = só no modelo, azul = só na imagem).

## Outros modos
| modo | eixo z | uso |
|---|---|---|
| `composto` (padrão) | sólido comprimido + relevo externo pelos tons + interior dividido em adiposo/denso | protótipo com variação interna |
| `solido` | espessura de compressão (`BodyPartThickness`), perfil semicircular na borda | sólido fechado em mm, cor = intensidade |
| `relevo` | intensidade normalizada × `--escala-relevo` | visualização tipo "terreno" |
| `volume` | real (DBT multi-frame) | marching cubes; automático se o DICOM tiver vários quadros |

**Limitação:** de uma mamografia 2D, a forma 3D é uma hipótese geométrica, não uma
reconstrução anatômica. A profundidade é perdida na projeção.

## Modo composto (interior pelos tons de cinza)
Modelo de duas componentes: cada pixel é uma coluna de espessura `t` com uma fração
`p` (tom de cinza normalizado) de tecido denso. O tecido denso ocupa uma camada
centrada de espessura `p·t`.
- A quantidade de tecido denso por coluna reproduz o tom de cinza.
- A posição em z é hipótese: a projeção não informa a profundidade.

Saídas:
- `impressao.stl`: peça inteira, para um único material.
- `partes/adiposo.stl` e `partes/denso*.stl`: mesmo referencial, para importar juntas.
- `impressao_multimaterial.3mf`: partes já agrupadas. No fatiador, atribua um
  material ou infill diferente a cada parte.

Parâmetros:
- `--niveis 4`: subdivide o denso em 3 faixas de cinza.
- `--relevo-superficie 3`: relevo externo em mm (0 desliga).

**Atenção:** `p` é normalizado pelos percentis 1–99 da própria imagem, não é densidade
calibrada. A fração densa resultante fica perto de 50 % e deve ser recalibrada com
a literatura de referência.

## Imagens RAW ("FOR PROCESSING")
Em imagem RAW o ar é o pixel mais claro e o tecido é escuro. O script detecta isso
pela borda da imagem e converte para atenuação, `A = ln(I0/I)` (Beer–Lambert, como
no mammoreplicator). `info.json → convertida_para_atenuacao` indica quando a
conversão foi aplicada.

## Referência: mammoreplicator (FDA)
Badal, Clark & Ghammraoui, *J. Med. Imaging* 5:033501 (2018). O código original
está em `externo/mammoreplicator` (Python 2 + PyPy, arquivado).

Para gerar e rodar o porte Python 3 (só mudanças mecânicas, registradas em
`externo/mammoreplicator_py3.diff`):

```bash
.venv\Scripts\python tests\portar_mammoreplicator.py
.venv\Scripts\python tests\sintetico_raw.py tmp_raw\raw_CC.dcm
.venv\Scripts\python externo\mammoreplicator_py3.py params.json
.venv\Scripts\python tests\validar_mammoreplicator.py <saida.ply> tmp_raw\raw_CC_gabarito.npz params.json
```

## Impressão 3D
`impressao.stl` sai em **mm**, já orientado e apoiado em z = 0. `info.json → impressao`
traz a validação (estanque, normais, corpo único, faces degeneradas, cabe na mesa).

- `solido`: em pé, com a face da parede torácica na mesa. Quase sem balanços: dispensa suportes.
- `relevo`: deitado sobre a base plana (`--base`, padrão 3 mm).
- STL não guarda cor. Em FDM a peça sai em uma cor só.
- O fatiador preenche o interior com infill; não é preciso ocar o modelo.
  Em resina (SLA), oque no fatiador e abra furos de drenagem.

```bash
.venv\Scripts\python src\mama3d.py dados --mesa 256 256 256          # volume da sua impressora
.venv\Scripts\python src\mama3d.py dados --escala-impressao 0.5      # meia escala
```

## Fusion 360
*Insert → Insert Mesh*, unidade **mm**. Para editar como sólido (*Convert Mesh*),
reduza antes o número de faces:

```bash
.venv\Scripts\python src\mama3d.py dados --max-faces 20000
```

A redução usa decimação por quádricas (Garland & Heckbert, 1997): elimina faces
primeiro nas regiões planas. Nos testes sintéticos, de 153 mil para 5 a 50 mil faces,
o desvio ficou abaixo da resolução da medição (cerca de 0,1 mm) e o volume variou
menos de 0,03 %.

## Parâmetros úteis
- `--voxel 0.5`: malha mais fina (≈4× mais faces)
- `--limiar 0.7`: se a segmentação cortar tecido adiposo perto da pele
- `--espessura 55`: força H quando o DICOM não tem `BodyPartThickness`

## Teste sem dados reais
```bash
.venv\Scripts\python tests\gerar_sintetico.py tmp_sint
.venv\Scripts\python src\mama3d.py tmp_sint
```
