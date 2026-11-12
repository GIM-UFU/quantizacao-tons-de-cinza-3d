"""Leitura de mamografias em DICOM ou TIFF, com as correções necessárias para análise.

Por que cada passo existe:
- Modality LUT (DICOM): converte o valor armazenado em valor "real" (RescaleSlope/Intercept).
- MONOCHROME1 (DICOM) / MinIsWhite (TIFF): branco = valor baixo. Invertemos para que,
  em toda a pipeline, tecido denso (mais atenuante) = valor ALTO.
- Espaçamento de pixel: sem ele a malha não tem escala em mm. O DICOM traz no
  cabeçalho; o TIFF nem sempre (ver `_pixel_mm_tiff`).
"""
from __future__ import annotations

import glob
import os
import warnings
from dataclasses import dataclass, field

import numpy as np
import pydicom
import tifffile
from pydicom.pixels import apply_modality_lut

PIXEL_PADRAO_MM = 0.1  # típico de mamografia digital (50–100 µm)


@dataclass
class Mamografia:
    caminho: str
    imagem: np.ndarray             # 2D (linhas, colunas) ou 3D (fatias, linhas, colunas)
    pixel_mm: tuple[float, float]  # (linha, coluna)
    espessura_mm: float | None     # espessura da mama comprimida (só DICOM)
    fatia_mm: float | None = None  # só para volumes (DBT)
    meta: dict = field(default_factory=dict)

    @property
    def eh_volume(self) -> bool:
        return self.imagem.ndim == 3


# ─────────────────────────────── DICOM ───────────────────────────────

def _f(ds, tag, default=None):
    v = ds.get(tag, None)
    if v is None or v == "":
        return default
    return v


def ler_metadados(ds) -> dict:
    return {
        "Formato": "DICOM",
        "Fabricante": str(_f(ds, "Manufacturer", "")),
        "Modalidade": str(_f(ds, "Modality", "")),
        "Intencao": str(_f(ds, "PresentationIntentType", "")),  # FOR PROCESSING / FOR PRESENTATION
        "Vista": str(_f(ds, "ViewPosition", "")),               # CC / MLO
        "Lateralidade": str(_f(ds, "ImageLaterality", _f(ds, "Laterality", ""))),
        "Fotometria": str(_f(ds, "PhotometricInterpretation", "")),
        "Linhas": int(_f(ds, "Rows", 0)),
        "Colunas": int(_f(ds, "Columns", 0)),
        "Quadros": int(_f(ds, "NumberOfFrames", 1)),
        "BitsArmazenados": int(_f(ds, "BitsStored", 0)),
        "PixelSpacing": [float(x) for x in _f(ds, "PixelSpacing", [])],
        "ImagerPixelSpacing": [float(x) for x in _f(ds, "ImagerPixelSpacing", [])],
        "EspessuraMama_mm": float(_f(ds, "BodyPartThickness", 0)) or None,
        "ForcaCompressao_N": float(_f(ds, "CompressionForce", 0)) or None,
        "kVp": float(_f(ds, "KVP", 0)) or None,
        "Exposicao_mAs": float(_f(ds, "Exposure", 0)) or None,
        "TransferSyntax": str(getattr(ds.file_meta, "TransferSyntaxUID", "")),
    }


def _carregar_dicom(caminho: str) -> Mamografia:
    ds = pydicom.dcmread(caminho)
    meta = ler_metadados(ds)

    img = apply_modality_lut(ds.pixel_array, ds).astype(np.float32)
    if meta["Fotometria"] == "MONOCHROME1":
        img = img.max() - img

    esp = meta["PixelSpacing"] or meta["ImagerPixelSpacing"]
    pixel_mm = (esp[0], esp[1]) if len(esp) == 2 else None
    meta["FontePixel"] = "cabeçalho DICOM" if pixel_mm else "ausente"

    fatia_mm = None
    if img.ndim == 3:
        fatia_mm = float(_f(ds, "SpacingBetweenSlices", 0)) or float(_f(ds, "SliceThickness", 0)) or 1.0

    return Mamografia(caminho, img, pixel_mm, meta["EspessuraMama_mm"], fatia_mm, meta)


# ─────────────────────────────── TIFF ───────────────────────────────

_PARA_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "inch": 25.4, "in": 25.4,
            "micron": 1e-3, "um": 1e-3, "µm": 1e-3, "µm": 1e-3, "nm": 1e-6}


def _racional(v) -> float:
    """Tags de resolução vêm como (numerador, denominador)."""
    if isinstance(v, tuple):
        return v[0] / v[1] if v[1] else 0.0
    return float(v)


def _pixel_mm_tiff(tf: tifffile.TiffFile) -> tuple[tuple[float, float] | None, float | None, str]:
    """Descobre o tamanho do pixel em mm, em ordem de confiabilidade:
    1. metadados do ImageJ ('unit' + XResolution em pixels/unidade) — inclui o
       espaçamento entre fatias ('spacing') em pilhas;
    2. ResolutionUnit do TIFF em cm/polegada — MAS 72 e 96 dpi são padrões de tela
       que muitos programas gravam sem significado físico, então são ignorados;
    3. nada → None (quem chama decide o padrão e avisa)."""
    p = tf.pages[0]
    tags = p.tags
    xr = _racional(tags["XResolution"].value) if "XResolution" in tags else 0.0
    yr = _racional(tags["YResolution"].value) if "YResolution" in tags else xr

    ij = tf.imagej_metadata or {}
    unidade = str(ij.get("unit", "")).strip().lower().replace("\\u00b5", "µ")
    if unidade in _PARA_MM and xr > 0:
        k = _PARA_MM[unidade]
        z = float(ij["spacing"]) * k if "spacing" in ij else None
        return (k / yr, k / xr), z, f"metadados ImageJ ({unidade})"

    ru = int(tags["ResolutionUnit"].value) if "ResolutionUnit" in tags else 2
    if xr > 0 and ru in (2, 3):
        if ru == 2 and round(xr) in (1, 72, 96):
            return None, None, f"ignorado ({xr:.0f} dpi é padrão de tela, não medida física)"
        k = 25.4 if ru == 2 else 10.0
        return (k / yr, k / xr), None, f"tag ResolutionUnit ({'dpi' if ru == 2 else 'px/cm'} = {xr:g})"
    return None, None, "ausente"


def _carregar_tiff(caminho: str) -> Mamografia:
    with tifffile.TiffFile(caminho) as tf:
        arr = tf.asarray()
        p = tf.pages[0]
        fotometria = p.photometric.name if p.photometric is not None else ""
        pixel_mm, fatia_mm, fonte = _pixel_mm_tiff(tf)
        bits = p.bitspersample

    img = np.asarray(arr)
    # RGB/RGBA (ex.: exportado de visualizador) → luminância
    if img.ndim >= 3 and img.shape[-1] in (3, 4) and fotometria in ("RGB", "YCBCR"):
        img = img[..., :3].astype(np.float32) @ np.array([0.299, 0.587, 0.114], np.float32)
    img = np.squeeze(img).astype(np.float32)
    if img.ndim > 3:
        raise RuntimeError(f"TIFF com {img.ndim} dimensões {img.shape}; esperado 2D ou pilha 3D.")
    if fotometria == "MINISWHITE":
        img = img.max() - img

    if img.ndim == 3 and fatia_mm is None:
        fatia_mm = 1.0
        warnings.warn(f"{os.path.basename(caminho)}: pilha sem espaçamento entre fatias; assumindo 1 mm.")

    meta = {
        "Formato": "TIFF", "Fabricante": "", "Modalidade": "", "Intencao": "desconhecida",
        "Vista": "", "Lateralidade": "", "Fotometria": fotometria,
        "Linhas": int(img.shape[-2]), "Colunas": int(img.shape[-1]),
        "Quadros": int(img.shape[0]) if img.ndim == 3 else 1,
        "BitsArmazenados": int(bits), "EspessuraMama_mm": None, "kVp": None,
        "FontePixel": fonte,
    }
    return Mamografia(caminho, img, pixel_mm, None, fatia_mm, meta)


# ─────────────────────────────── geral ───────────────────────────────

def carregar(caminho: str, pixel_mm: float | None = None) -> Mamografia:
    """Lê DICOM ou TIFF (detectado pelo conteúdo, não pela extensão).
    `pixel_mm` força o tamanho do pixel — necessário para TIFF sem calibração."""
    m = _carregar_tiff(caminho) if _eh_tiff(caminho) else _carregar_dicom(caminho)
    if pixel_mm:
        m.pixel_mm = (pixel_mm, pixel_mm)
        m.meta["FontePixel"] = "argumento --pixel-mm"
    elif m.pixel_mm is None:
        m.pixel_mm = (PIXEL_PADRAO_MM, PIXEL_PADRAO_MM)
        warnings.warn(f"{os.path.basename(caminho)}: tamanho do pixel desconhecido "
                      f"({m.meta['FontePixel']}); assumindo {PIXEL_PADRAO_MM} mm. "
                      "Use --pixel-mm para a escala correta.")
        m.meta["FontePixel"] += f" → padrão {PIXEL_PADRAO_MM} mm"
    return m


def listar_imagens(entradas: list[str]) -> list[str]:
    """Aceita arquivos, pastas e curingas (o PowerShell não expande *.dcm)."""
    arquivos: list[str] = []
    for e in entradas:
        if os.path.isdir(e):
            candidatos = sorted(glob.glob(os.path.join(e, "**", "*"), recursive=True))
        else:
            candidatos = sorted(glob.glob(e)) or [e]
        for c in candidatos:
            if os.path.isfile(c) and (_eh_tiff(c) or _eh_dicom(c)):
                arquivos.append(c)
    return arquivos


listar_dicoms = listar_imagens  # compatibilidade


def _eh_tiff(caminho: str) -> bool:
    try:
        with open(caminho, "rb") as f:
            return f.read(4) in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")  # TIFF e BigTIFF
    except OSError:
        return False


def _eh_dicom(caminho: str) -> bool:
    try:
        with open(caminho, "rb") as f:
            f.seek(128)
            if f.read(4) == b"DICM":
                return True
        pydicom.dcmread(caminho, stop_before_pixels=True)
        return True
    except Exception:
        return False
