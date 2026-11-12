"""Gera uma versão Python 3 do mammoreplicator SEM editar o original.

O original (externo/mammoreplicator/mammoreplicator.py) é Python 2 + PyPy (`numpypy`)
e usa o pacote antigo `dicom`. Aqui aplicamos só substituições mecânicas, cada uma
listada abaixo, para que a lógica numérica permaneça idêntica e o diff seja auditável:

  1. print 'x', y            → print('x', y)                 (sintaxe Py3)
  2. import numpypy ...      → import numpy as numpypy ...   (alternativa do próprio autor, linha 16)
  3. import dicom            → import pydicom as dicom; dicom.read_file → dicom.dcmread
  4. N/output_binning        → N//output_binning              (divisão inteira do Py2)
  5. mammo.pixel_array       → PIX (mesmo array, lido uma vez; evita reprocessar a cada acesso)
  6. parâmetros fixos no código (arquivo, espessura, intensidades, linhas/colunas,
     binning, pixel) → lidos de um JSON passado na linha de comando.

Uso: python tests/portar_mammoreplicator.py  → escreve externo/mammoreplicator_py3.py
"""
import difflib
import os
import re

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIG = os.path.join(RAIZ, "externo", "mammoreplicator", "mammoreplicator.py")
DEST = os.path.join(RAIZ, "externo", "mammoreplicator_py3.py")

PARAMS = ["image_file_name", "compression_thickness", "average_breast_intensity",
          "air_intensity", "output_binning", "rows", "columns", "pixel_size",
          "mfp_breast", "mfp_plastic"]


def portar(src: str) -> str:
    out = []
    for ln in src.splitlines():
        s = ln
        # 2. numpy no lugar de numpypy
        if s.startswith("import numpypy;"):
            s = "#" + s
        elif s.startswith("#import numpy as numpypy"):
            s = s[1:]
        # 3. pydicom
        if s.strip() == "import dicom":
            s = "import pydicom as dicom; import json; P = json.load(open(sys.argv[1]))"
        s = s.replace("dicom.read_file(", "dicom.dcmread(")
        # 6. parâmetros
        m = re.match(r"^(\w+)\s*=", s)
        if m and m.group(1) in PARAMS:
            s = f"{m.group(1)} = P['{m.group(1)}'] if '{m.group(1)}' in P else {s.split('=', 1)[1].split('#')[0].strip()}"
        # 4. divisão inteira
        s = re.sub(r"(num_columns|num_rows)/output_binning", r"\1//output_binning", s)
        # 5. array de pixels
        s = s.replace("mammo.pixel_array", "PIX")
        # 1. print statement → função (só linhas cujo código começa com print)
        pm = re.match(r"^(\s*)print\s+(.*)$", s)
        if pm and not pm.group(2).startswith("("):
            s = f"{pm.group(1)}print({pm.group(2)})"
        out.append(s)
        if s.startswith("mammo = dicom.dcmread("):
            out.append("PIX = mammo.pixel_array  # 5. mesmo array do DICOM")
    # `import sys` vem depois de `import dicom` no original: move o json/argv para depois
    txt = "\n".join(out) + "\n"
    txt = txt.replace("import pydicom as dicom; import json; P = json.load(open(sys.argv[1]))",
                      "import pydicom as dicom")
    txt = txt.replace("import sys\n", "import sys\nimport json\nP = json.load(open(sys.argv[1]))\n", 1)
    return txt


if __name__ == "__main__":
    with open(ORIG, encoding="utf-8") as f:
        src = f.read()
    novo = portar(src)
    with open(DEST, "w", encoding="utf-8") as f:
        f.write(novo)
    diff = list(difflib.unified_diff(src.splitlines(), novo.splitlines(), "original", "py3", n=0, lineterm=""))
    with open(DEST.replace(".py", ".diff"), "w", encoding="utf-8") as f:
        f.write("\n".join(diff) + "\n")
    print(f"porte escrito em {DEST}  ({sum(1 for d in diff if d.startswith('+') and not d.startswith('+++'))} linhas alteradas)")
