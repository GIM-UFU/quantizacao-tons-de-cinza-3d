"""Mamografia RAW ("FOR PROCESSING") sintética com gabarito físico.

Física (Beer–Lambert, feixe monoenergético — simplificação):
    I(x,y) = I0 · exp(−t(x,y) · μ(x,y)),   μ = (1−f)·μ_adiposo + f·μ_glandular
- t: espessura local (mm→cm), mesma geometria de compressão da pipeline.
- f: fração glandular (gabarito, salvo em .npz).
- μ ≈ 0.46 cm⁻¹ (adiposo) e 0.80 cm⁻¹ (glandular) — ordem de grandeza perto de
  20 keV (cf. Johns & Yaffe, Phys. Med. Biol. 1987); valores aproximados, só para teste.
Em RAW o ar é o pixel MAIS CLARO (I0) e o tecido denso é ESCURO — o oposto de uma
imagem "FOR PRESENTATION".

Uso: python tests/sintetico_raw.py <saida.dcm>
"""
import os
import sys

import numpy as np
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid
from scipy import ndimage as ndi

I0, T_MM, PX_MM = 16383, 50.0, 0.1
MU_A, MU_G = 0.46, 0.80  # cm⁻¹


def gerar(h=2304, w=1912, seed=7):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:h, 0:w]
    mask = (xx / (w * 0.75)) ** 2 + ((yy - h / 2) / (h * 0.42)) ** 2 <= 1
    d = ndi.distance_transform_edt(mask) * PX_MM
    R = T_MM / 2
    t = 2 * np.sqrt(np.maximum(R**2 - (R - np.minimum(d, R)) ** 2, 0))  # mm
    g = ndi.gaussian_filter(rng.random((h, w)), 30)
    g = (g - g.min()) / (g.max() - g.min())
    f = np.clip(1.6 * g - 0.5, 0, 0.6) * mask  # fração glandular 0–60 %
    mu = (1 - f) * MU_A + f * MU_G
    I = I0 * np.exp(-(t / 10) * mu)
    I = rng.poisson(I / 4) * 4.0  # ruído quântico aproximado
    return np.clip(I, 0, I0).astype(np.uint16), t, f, mask


def salvar(caminho, pix):
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.1.2.1"  # DX Mammo For Processing
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = FileDataset(caminho, {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID, ds.SOPInstanceUID = meta.MediaStorageSOPClassUID, meta.MediaStorageSOPInstanceUID
    ds.Modality, ds.Manufacturer = "MG", "SINTETICO_RAW"
    ds.PatientName, ds.PatientID = "Teste^Raw", "000"
    ds.PresentationIntentType = "FOR PROCESSING"
    ds.ViewPosition, ds.ImageLaterality = "CC", "L"
    ds.BodyPartThickness = T_MM
    ds.ImagerPixelSpacing = [PX_MM, PX_MM]
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.SamplesPerPixel, ds.BitsAllocated, ds.BitsStored, ds.HighBit = 1, 16, 14, 13
    ds.PixelRepresentation = 0
    ds.Rows, ds.Columns = pix.shape
    ds.PixelData = pix.tobytes()
    ds.save_as(caminho, enforce_file_format=True)


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    pix, t, f, mask = gerar()
    salvar(out, pix)
    np.savez_compressed(out.replace(".dcm", "_gabarito.npz"), t=t, f=f, mask=mask)
    centro = mask & (t >= T_MM - 1e-6)
    print(f"ok: {pix.shape}, ar={pix[~mask].mean():.0f}, mama média(centro)={pix[centro].mean():.0f}, "
          f"f média={f[mask].mean():.3f}")
