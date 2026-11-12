# Imprimindo as placas na Ender 3

Valores de referência da Creality Ender-3: área de impressão 220 × 220 × 250 mm, bico de
0,4 mm, cartão microSD para os arquivos. As recomendações de fatiamento abaixo são nossas
(a literatura só indica 100 % de preenchimento, Bliznakova et al. 2024); ajuste com a
radiografia do material (veja [CALIBRACAO_RADIOGRAFICA.md](CALIBRACAO_RADIOGRAFICA.md)).

## 1. Como gerar as placas para a Ender 3

```bash
.venv\Scripts\python src\mama3d.py dados\lcc1.tif --pixel-mm 0.07 --voxel 0.25 --suavizacao 0.4 --traco-min 0.8 --area-min 10 --camada 0.2 --mesa 220 220 250 --saida saida\placas_ender3
```

- `--camada 0.2`: a espessura da placa vira múltiplo da altura de camada (3,4 mm = 17 camadas).
- `--traco-min 0.8`: dois filetes de 0,4 mm; `--area-min 10`: ilhas de pelo menos 10 mm².
- `--voxel 0.25` e `--suavizacao 0.4`: mais detalhe (arquivos maiores: a placa 01 tem 21 MB).

Comparação das configurações testadas no `lcc1.tif` (variação do relevo preservada no STL
final, volume retido e nº de peças soltas; "A" é a configuração anterior):

| Config. | Placas × espessura | Peças | Variação preservada | Volume retido |
|---|---|---|---|---|
| A: voxel 0,5, suavização 0,75, traço 1,2, área 25 | 14 × 3,2 mm | 44 | 86,9 % | 90,6 % |
| B: A com suavização 0,4 | 15 × 3,2 mm | 45 | 88,5 % | 87,3 % |
| **C: voxel 0,25, suavização 0,4, traço 0,8, área 10, camada 0,2** | **14 × 3,4 mm** | **98** | **90,0 %** | **92,3 %** |
| D: C com placas de ~2,6 mm | 18 × 2,6 mm | 127 | 90,1 % | 92,4 % |
| E: C com placas de 2,0 mm | 24 × 2,0 mm | 161 | 90,2 % | 92,5 % |

- C ganha 3,1 pontos sobre A e também no acerto contra a imagem (92,4 % contra 90,6 %).
- **Placas mais finas (D, E) quase não ganham fidelidade** e multiplicam placas e peças.
  Por isso a configuração recomendada é a C.
- O custo de C são mais peças soltas (98 contra 44). Se as menores (cerca de 3 × 3 mm)
  soltarem da mesa, suba `--area-min` para 25.

## 2. Qual arquivo usar

Use o **STL gerado pelo script**, direto no fatiador:

```
saida/placas_ender3/lcc1/placas/placa_01.stl   (e placa_02.stl, placa_03.stl, ...)
```

- Já está em **milímetros**, **deitado** (espessura na vertical) e com a base em z = 0.
- A placa 01 mede **90,0 × 121,5 × 3,4 mm**, tem 14 peças e 12,52 cm³. (A medida feita no
  Fusion, 93,0 × 158,8 × 3,19 mm, 5 peças, é de uma versão anterior.)
- **Correção da modelagem:** a borda da imagem tem uma faixa clara de ~2–3 mm que não é tecido e virava
  lascas de placa fora do contorno da mama. Agora há uma margem de 3 mm na borda (--margem-borda) e o
  gerador confere que nenhum vértice das placas fica fora da mama (ertices_das_placas_fora_da_mama = 0).
- O Fusion não é necessário para imprimir. Serve para inspecionar ou editar. O eixo
  "para cima" do Fusion é o Y; os fatiadores usam Z. Por isso a cópia aberta no Fusion
  (`fusion_yup/placa_01_yup.stl`) **não** deve ir para o fatiador.
- Se editar a placa no Fusion: clique com o botão direito no corpo, *Save As Mesh*, formato
  **STL (binário)**, unidade **milímetro** (o nome do menu pode variar com a versão).
  Depois confira no fatiador que a placa aparece deitada e mede 90 × 121,5 × 3,4 mm.

A placa 01 cabe na mesa de 220 × 220 mm; duas placas lado a lado também cabem, mas imprima
**uma por vez** para controlar a qualidade.

## 3. Fatiador (Cura ou Creality Print)

1. Selecione a máquina **Creality Ender-3** (já existe no Cura).
2. Abra o STL. **Não gire nem escale** o modelo; mantenha 100 %.
3. Ajustes recomendados:

| Ajuste | Valor | Motivo |
|---|---|---|
| Altura de camada | **0,2 mm** | a placa de 3,4 mm = 17 camadas exatas |
| Preenchimento | **100 %** | a radiografia deve ver material sólido |
| Paredes | 2 a 3 | |
| Suportes | **nenhum** | a placa é plana |
| Aderência | *brim* de 3 a 5 mm (ou *skirt*) | ilhas pequenas soltam da mesa |
| Primeira camada | mais lenta (cerca de 20 mm/s) | |
| Fluxo | 100 % para começar | só aumente com base na radiografia do degrau |
| Ventilador | padrão do perfil | |

4. Use o **mesmo filamento, a mesma temperatura e os mesmos ajustes** em todas as placas e
   nas placas-teste do degrau. Seque o PLA antes de imprimir: a umidade gera bolhas de ar
   (Bliznakova et al. 2024).
5. Confira a pré-visualização camada a camada: a placa 01 tem 14 peças soltas, e todas
   precisam estar presentes.

## 4. Salvar para a impressora

- No Cura, *Salvar em disco* (ou *Salvar no dispositivo removível*) gera um arquivo
  **`.gcode`**; é ele que a Ender 3 lê. O STL não vai para a impressora.
- Nome curto, sem acentos nem espaços: `placa01.gcode`.
- Cartão microSD formatado em **FAT32**.

## 5. Depois de imprimir

- Meça com paquímetro: o esperado é **90,0 × 121,5 × 3,4 mm**, com tolerância de cerca de
  ±0,3 mm. Diferenças maiores indicam calibração de passos, fluxo ou encolhimento.
- Olhe a placa contra a luz: bolhas, falhas de trilha e microfrestas aparecem.
- Guarde as placas numeradas; a ordem de empilhamento e o `guia_montagem_1para1.pdf`
  definem a posição de cada peça.

## 6. Estimativa de filamento

```bash
.venv\Scripts\python src\filamento.py saida\placas_ender3\lcc1\placas\placa_01.stl --preco-kg 100
```

Placa 01 (PLA 1,75 mm, ρ = 1,24 g/cm³, 100 % de preenchimento, brim de 4 mm):

| Cenário | Massa | Comprimento | % de 1 kg |
|---|---|---|---|
| Só o modelo (12,52 cm³) | 15,5 g | 5,20 m | 1,6 % |
| Modelo + brim + skirt + purga | 17,0 g | 5,69 m | 1,7 % |
| Com reserva de 10 % | 18,7 g | 6,26 m | 1,9 % |

O custo do exemplo usa R$ 100/kg, valor ilustrativo; troque pelo do seu filamento.
O tempo (≈ 1,2 h) é só ordem de grandeza. O fatiador dá o número final. As 14 placas
somam 58,2 cm³, ou cerca de 72 g só de modelo.

## 7. G-code pronto (mesmo formato do pendrive)

Os arquivos do pendrive são `.gcode` (Marlin, extrusão absoluta), gerados pelo Cura 5.13 com o
perfil Ender-3 e nomeados `CE3E3V2_<nome>.gcode` (a convenção do Creality Print). Para a placa 01:

```bash
.venv\Scripts\python src\gcode.py saida\placas_ender3\lcc1\placas\placa_01.stl saida\gcode\CE3E3V2_placa01.gcode
.venv\Scripts\python tests\verificar_gcode.py saida\gcode\CE3E3V2_placa01.gcode saida\placas_ender3\lcc1\placas\placa_01.stl
```

- **Início e fim do arquivo:** copiados dos arquivos do pendrive (PLA a 200 °C e mesa a 60 °C, linhas de
  purga, retração, ventilador em rampa nas camadas 1 a 3, final com `M84`).
- **Fatiamento:** 17 camadas de 0,2 mm, brim de 4 mm, duas paredes e preenchimento 100 % a 45°/135°,
  peça centralizada na mesa (110, 110). Retração de 5 mm só em deslocamentos que saem da peça.
- **Conferência feita lendo o arquivo:** limites dentro da mesa de 220 × 220 mm, Z de 0,2 a 3,4 mm,
  filamento 5,54 m (16,5 g), linhas cobrem 98,2 % da seção e só 0,06 % do material cai fora dela.
- **Atenção:** este G-code vem de um gerador próprio e **nunca rodou numa impressora**. Abra o arquivo no
  Creality Print (ou outro visualizador de G-code) e confira a prévia camada a camada antes de imprimir,
  e acompanhe a primeira camada. O Cura continua sendo a referência; use-o se preferir.