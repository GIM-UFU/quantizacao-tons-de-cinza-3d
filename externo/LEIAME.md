# Código de terceiros (não versionado)

Pasta ignorada pelo git. Para reproduzir o ambiente:

```bash
# mammoreplicator (FDA, domínio público; repositório arquivado)
git clone --depth 1 https://github.com/DIDSR/mammoreplicator externo/mammoreplicator

# porte para Python 3 (só mudanças mecânicas, com diff em externo/mammoreplicator_py3.diff)
.venv\Scripts\python tests\portar_mammoreplicator.py
```

**ImageJ 1.54 (Windows, com Java 8)**: baixar `ij154-win-java8.zip` de
<https://wsr.imagej.net/distros/win/ij154-win-java8.zip> e extrair em `externo/imagej/`
(resulta em `externo/imagej/ImageJ/ij.jar`). SHA-256 do zip usado:
`32B05356D3A7906E10CA91491C2BE06BFE0397C9E7A3D3B52FDAD9E3EF48276C`.
