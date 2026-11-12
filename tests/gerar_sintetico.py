"""Gera DICOMs sintéticos de mamografia só para testar a pipeline (não são dados reais).

Uso: python tests/gerar_sintetico.py <pasta_saida>
"""
import os
import sys

import numpy as np
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid
from scipy import ndimage as ndi

DX_MAMMO = "1.2.840.10008.5.1.4.1.1.1.2"


def mama_2d(h=2400, w=1800, seed=0):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:h, 0:w]
    # semi-elipse apoiada na parede torácica (x = 0)
    mask = (xx / (w * 0.75)) ** 2 + ((yy - h / 2) / (h * 0.42)) ** 2 <= 1
    tecido = ndi.gaussian_filter(rng.random((h, w)), 25)
    tecido = (tecido - tecido.min()) / (tecido.max() - tecido.min())
    img = np.where(mask, 600 + 2200 * tecido + 300 * (1 - xx / w), 80)
    img += rng.normal(0, 15, img.shape)
    img[50:120, w - 300:w - 100] = 3800  # marcador/etiqueta no fundo
    return np.clip(img, 0, 4095).astype(np.uint16)


def salvar(caminho, pix, vista, lado, mono1=False, frames=1):
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = DX_MAMMO
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = FileDataset(caminho, {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID, ds.SOPInstanceUID = DX_MAMMO, meta.MediaStorageSOPInstanceUID
    ds.Modality, ds.Manufacturer = "MG", "SINTETICO"
    ds.PatientName, ds.PatientID = "Teste^Sintetico", "000"
    ds.PresentationIntentType = "FOR PROCESSING"
    ds.ViewPosition, ds.ImageLaterality = vista, lado
    ds.BodyPartThickness, ds.KVP = 52, 29
    ds.ImagerPixelSpacing = [0.1, 0.1]
    if mono1:
        pix = (4095 - pix).astype(np.uint16)
    ds.PhotometricInterpretation = "MONOCHROME1" if mono1 else "MONOCHROME2"
    ds.SamplesPerPixel, ds.BitsAllocated, ds.BitsStored, ds.HighBit = 1, 16, 12, 11
    ds.PixelRepresentation = 0
    if frames > 1:
        ds.NumberOfFrames, ds.SpacingBetweenSlices = frames, 1.0
        ds.Rows, ds.Columns = pix.shape[1:]
    else:
        ds.Rows, ds.Columns = pix.shape
    ds.PixelData = pix.tobytes()
    ds.save_as(caminho, enforce_file_format=True)


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    salvar(os.path.join(out, "L_CC.dcm"), mama_2d(seed=1), "CC", "L")
    salvar(os.path.join(out, "R_CC.dcm"), mama_2d(seed=2)[:, ::-1], "CC", "R", mono1=True)
    salvar(os.path.join(out, "L_MLO.dcm"), mama_2d(seed=3), "MLO", "L")
    # DBT pequeno: 40 fatias de 1 mm, pixel 0.5 mm
    base = mama_2d(480, 360, seed=4).astype(np.float32)
    z = np.arange(40)[:, None, None]
    vol = np.where(np.abs(z - 20) < 18, base[None], 80).astype(np.uint16)
    salvar(os.path.join(out, "L_DBT.dcm"), vol, "CC", "L", frames=40)
    gerar_tiffs(os.path.join(out, "tiff"))
    print("ok:", os.listdir(out))


def gerar_tiffs(out):
    """Cobre as variações de TIFF que aparecem na prática."""
    import tifffile
    os.makedirs(out, exist_ok=True)
    img16 = mama_2d(seed=5)
    img8 = (img16 >> 4).astype(np.uint8)
    # 1. ImageJ 16 bits, calibrado em mm (10 px/mm = 0.1 mm)
    tifffile.imwrite(os.path.join(out, "1_imagej_mm.tif"), img16, imagej=True,
                     resolution=(10, 10), metadata={"unit": "mm"})
    # 2. 8 bits, LZW, 72 dpi (padrão de tela → deve ser ignorado)
    tifffile.imwrite(os.path.join(out, "2_tela_72dpi.tif"), img8, compression="lzw",
                     resolution=(72, 72), resolutionunit="INCH")
    # 3. RGB, 254 dpi (= 0.1 mm), mama espelhada
    tifffile.imwrite(os.path.join(out, "3_rgb_254dpi.tiff"), np.repeat(img8[:, ::-1, None], 3, 2),
                     photometric="rgb", resolution=(254, 254), resolutionunit="INCH")
    # 4. MinIsWhite (equivalente ao MONOCHROME1), 100 px/cm
    tifffile.imwrite(os.path.join(out, "4_miniswhite_cm.tif"), (4095 - img16).astype(np.uint16),
                     photometric="miniswhite", resolution=(100, 100), resolutionunit="CENTIMETER")
    # 5. Pilha ImageJ (DBT exportada): 40 fatias, pixel 0.5 mm em µm, fatia 1 mm
    base = mama_2d(480, 360, seed=6)
    z = np.arange(40)[:, None, None]
    vol = np.where(np.abs(z - 20) < 18, base[None], 80).astype(np.uint16)
    tifffile.imwrite(os.path.join(out, "5_pilha_dbt.tif"), vol, imagej=True,
                     resolution=(1 / 500, 1 / 500), metadata={"unit": "micron", "spacing": 1000.0, "axes": "ZYX"})


if __name__ == "__main__":
    main()
