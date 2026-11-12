// Comparacao modelo 3D x mamografia original, no ImageJ/Fiji.
//
// Entradas:
//   original  = mamografia (TIFF ou DICOM)
//   projecao  = comparacao/projecao_modelo.tif, gerado por src/comparar.py
//               (32 bits; cada pixel = mm de material denso do modelo, somados em z)
// Saidas (na pasta da projecao):
//   resultado_imagej.csv   metricas
//   mapa_acerto_imagej.png verde = acerto (denso nos dois), vermelho = so no modelo,
//                          azul = so na imagem
//
// Uso interativo: Plugins > Macros > Run... e escolher os dois arquivos.
// Uso em lote:    java -jar ij.jar -batch comparar_imagej.ijm "original|projecao|pasta_saida|margem_px"
//                 (linha completa em comparacao/comando_imagej.txt)
//
// Definicoes (identicas as de src/comparar.py, para os numeros serem comparaveis):
//   mama            = limiar automatico "Triangle" + preenchimento de buracos
//   denso na imagem = cinza > Otsu calculado so com os pixels da mama
//   denso no modelo = espessura projetada > 0
//   interior        = mama erodida de `margem` pixels (o modelo nao imprime a faixa da pele)

// globais: na linguagem de macro, funcoes so enxergam variaveis declaradas com `var`
var w, h, N, limiar;
var o, e;

arg = getArgument();
if (arg == "") {
    orig = File.openDialog("Mamografia original (TIFF/DICOM)");
    proj = File.openDialog("projecao_modelo.tif");
    out = File.getParent(proj) + File.separator;
    margem = getNumber("Faixa da pele a ignorar (pixels da projecao)", 10);
} else {
    p = split(arg, "|");
    orig = p[0]; proj = p[1]; out = p[2]; margem = parseInt(p[3]);
}

setBatchMode(true);
setOption("BlackBackground", true);

// --- 1. modelo e original na mesma grade ---
open(proj); rename("modelo");
w = getWidth(); h = getHeight();
open(orig); run("32-bit"); rename("orig_cheia");
run("Scale...", "x=- y=- width=" + w + " height=" + h + " interpolation=Bilinear average create title=original");
close("orig_cheia");

// imagem RAW (ar claro, tecido escuro)? -> inverte, para tecido = claro
selectImage("original");
makeRectangle(0, 0, w, h); getStatistics(areaT, mediaT);
soma = 0; n = 0; k = maxOf(2, floor(minOf(w, h) / 50));
for (x = 0; x < w; x++) { for (j = 0; j < k; j++) { soma += getPixel(x, j) + getPixel(x, h - 1 - j); n += 2; } }
for (y = 0; y < h; y++) { for (j = 0; j < k; j++) { soma += getPixel(w - 1 - j, y); n++; } }
run("Select None");
if (soma / n > mediaT) { run("Invert"); print("Imagem com ar claro (RAW): invertida."); }

// --- 2. mascara da mama e interior ---
selectImage("original"); run("Duplicate...", "title=mascara");
setAutoThreshold("Triangle dark");
run("Convert to Mask");
run("Fill Holes");
run("Duplicate...", "title=interior");
if (margem > 0) {
    run("Options...", "iterations=" + margem + " count=1 black pad");  // pad: nao erode a partir da borda da imagem
    run("Erode");
    run("Options...", "iterations=1 count=1 black");
}

// --- 3. le os pixels ---
N = w * h;
o = ler("original"); e = ler("modelo"); m = ler("mascara"); mi = ler("interior");

// --- 4. Otsu dentro da mama ---
mn = 1e30; mx = -1e30;
for (i = 0; i < N; i++) { if (m[i] > 0) { if (o[i] < mn) mn = o[i]; if (o[i] > mx) mx = o[i]; } }
nb = 256; hist = newArray(nb); tot = 0;
for (i = 0; i < N; i++) { if (m[i] > 0) { b = floor((o[i] - mn) / (mx - mn) * (nb - 1) + 0.5); hist[b]++; tot++; } }
somaT = 0; for (b = 0; b < nb; b++) somaT += b * hist[b];
w0 = 0; s0 = 0; melhor = -1; bOtsu = 0;
for (b = 0; b < nb - 1; b++) {
    w0 += hist[b]; s0 += b * hist[b]; w1 = tot - w0;
    if (w0 > 0 && w1 > 0) {
        d = s0 / w0 - (somaT - s0) / w1; entre = w0 * w1 * d * d;  // variancia entre classes
        if (entre > melhor) { melhor = entre; bOtsu = b; }
    }
}
limiar = mn + (bOtsu + 0.5) * (mx - mn) / (nb - 1);

// --- 5. metricas ---
csv = "regiao,pixels,limiar_otsu,pearson_r,r2,acerto_pixel_a_pixel,dice,iou,sensibilidade,precisao,TP,FP,FN,TN\n";
csv = csv + medir("mama_inteira", m);
csv = csv + medir("interior_sem_faixa_da_pele", mi);
File.saveString(csv, out + "resultado_imagej.csv");
print(csv);

// --- 6. mapa de acerto ---
newImage("mapa_acerto", "RGB black", w, h, 1);
for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
        i = y * w + x;
        g = floor(140 * (o[i] - mn) / (mx - mn)); if (g < 0) g = 0; if (g > 140) g = 140;
        r = g; gg = g; bb = g;
        if (m[i] > 0) {
            di = o[i] > limiar; dm = e[i] > 0;
            if (di && dm) { r = 40; gg = 200; bb = 60; }
            else if (dm) { r = 230; gg = 50; bb = 50; }
            else if (di) { r = 60; gg = 110; bb = 255; }
        }
        setPixel(x, y, (r << 16) + (gg << 8) + bb);
    }
}
updateDisplay();
saveAs("PNG", out + "mapa_acerto_imagej.png");
print("Saidas em: " + out);
if (arg != "") eval("script", "System.exit(0);");

function ler(titulo) {
    selectImage(titulo);
    a = newArray(w * h);
    for (y = 0; y < h; y++) for (x = 0; x < w; x++) a[y * w + x] = getPixel(x, y);
    return a;
}

function medir(nome, reg) {
    n = 0; so = 0; se = 0; soo = 0; see = 0; soe = 0; tp = 0; fp = 0; fn = 0; tn = 0;
    for (i = 0; i < N; i++) {
        if (reg[i] > 0) {
            n++; so += o[i]; se += e[i]; soo += o[i] * o[i]; see += e[i] * e[i]; soe += o[i] * e[i];
            di = o[i] > limiar; dm = e[i] > 0;
            if (di && dm) tp++; else if (dm) fp++; else if (di) fn++; else tn++;
        }
    }
    cov = soe / n - (so / n) * (se / n);
    r = cov / sqrt((soo / n - (so / n) * (so / n)) * (see / n - (se / n) * (se / n)));
    return nome + "," + n + "," + d2s(limiar, 2) + "," + d2s(r, 4) + "," + d2s(r * r, 4) + ","
        + d2s((tp + tn) / n, 4) + "," + d2s(2 * tp / (2 * tp + fp + fn), 4) + "," + d2s(tp / (tp + fp + fn), 4) + ","
        + d2s(tp / (tp + fn), 4) + "," + d2s(tp / (tp + fp), 4) + "," + tp + "," + fp + "," + fn + "," + tn + "\n";
}
