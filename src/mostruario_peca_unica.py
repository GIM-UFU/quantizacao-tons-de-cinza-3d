"""Mostruário 3D do phantom final de peça única: peça vazada (PLA) + cera de abelha + recipiente e tampa (resina).

Gera um HTML autônomo (GLB embutido, three.js por CDN) com:
  ligar/desligar cada parte · corte transversal (plano Y) · "até a camada N" da peça (camadas de 0,2 mm)
  · tampa aberta/fechada · números do info.json.

Uso:
    python src/mostruario_peca_unica.py saida/09_peca_unica/lcc1_pf06
Requer, na pasta: peca_unica.stl, recipiente.stl, tampa.stl, info.json, info_recipiente.json
(gerados por peca_unica.py e recipiente_peca_unica.py). A cera é reconstruída da máscara em
<origem>/intermediarios.npz. Grava <pasta>/mostruario.html.
"""
from __future__ import annotations

import base64
import io
import json
import os
import sys

import numpy as np
import trimesh
from scipy import ndimage as ndi
from skimage import morphology

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mama3d import reduzir, solido_entre  # noqa: E402
from peca_unica import reamostrar_media  # noqa: E402

FACES = {"peca": 160000, "recipiente": 50000, "tampa": 40000, "cera": 30000}

PAGINA = """<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Phantom peça única</title>
<style>
html,body{margin:0;height:100%;background:#101215;color:#ddd;font:13px system-ui;overflow:hidden}
#barra{position:absolute;left:0;right:0;top:0;padding:8px 12px;background:#000b;display:flex;gap:16px;
 align-items:center;flex-wrap:wrap;z-index:2}
#barra b{color:#fff}
#info{position:absolute;left:10px;bottom:10px;z-index:2;background:#000b;padding:8px 12px;border-radius:6px;
 line-height:1.55;font-size:12px;max-width:min(420px,calc(100% - 44px))}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:-1px}
input[type=range]{width:150px;vertical-align:middle}
</style>
<script type="importmap">{"imports":{
 "three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
 "three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}</script></head>
<body>
<div id="barra"><b>__NOME__ · peça única</b>
 <label><input id="v_peca" type="checkbox" checked><span class="sw" style="background:#d9665f"></span>peça (PLA)</label>
 <label><input id="v_cera" type="checkbox" checked><span class="sw" style="background:#e8b84a"></span>cera</label>
 <label><input id="v_recipiente" type="checkbox" checked><span class="sw" style="background:#9fc6e0"></span>recipiente (resina)</label>
 <label><input id="v_tampa" type="checkbox" checked><span class="sw" style="background:#b9d6ea"></span>tampa (resina)</label>
 <label><input id="aberta" type="checkbox" checked> tampa aberta</label>
 <span>corte <input id="corte" type="range" min="0" max="1000" value="1000"></span>
 <span>até a camada <input id="cam" type="range" min="1" max="__NC__" value="__NC__"> <span id="camv">__NC__</span>/__NC__</span>
 <span style="opacity:.7">arraste = girar · roda = zoom · botão direito = mover</span></div>
<div id="info">__INFO__</div>
<script type="module">
import * as THREE from "three";
import {OrbitControls} from "three/addons/controls/OrbitControls.js";
import {GLTFLoader} from "three/addons/loaders/GLTFLoader.js";
const B64="__GLB__",FUNDO=__FUNDO__,CAMADA=__CAMADA__,ABRE=__ABRE__;
const r=new THREE.WebGLRenderer({antialias:true});r.setPixelRatio(devicePixelRatio);
r.setSize(innerWidth,innerHeight);r.localClippingEnabled=true;document.body.appendChild(r.domElement);
const s=new THREE.Scene();s.background=new THREE.Color(0x101215);
s.add(new THREE.HemisphereLight(0xffffff,0x334455,1.5));
const d=new THREE.DirectionalLight(0xffffff,1.5);d.position.set(1,-1,2);s.add(d);
const d2=new THREE.DirectionalLight(0xffffff,.5);d2.position.set(-1,1,-1);s.add(d2);
const c=new THREE.PerspectiveCamera(40,innerWidth/innerHeight,0.1,5000);c.up.set(0,0,1);
const ctl=new OrbitControls(c,r.domElement);ctl.enableDamping=true;
const pCorte=new THREE.Plane(new THREE.Vector3(0,-1,0),0),pCam=new THREE.Plane(new THREE.Vector3(0,0,-1),0);
const MAT={
 peca:{color:0xd9665f,roughness:.55,opacity:1},
 cera:{color:0xe8b84a,roughness:.3,opacity:.28},
 recipiente:{color:0x9fc6e0,roughness:.15,opacity:.32},
 tampa:{color:0xb9d6ea,roughness:.15,opacity:.32}};
const P={};let ext=null,cen=null;
const bin=Uint8Array.from(atob(B64),ch=>ch.charCodeAt(0)).buffer;
new GLTFLoader().parse(bin,"",g=>{g.scene.traverse(n=>{if(!n.isMesh)return;
  const nm=(n.name+" "+(n.parent?n.parent.name:"")).toLowerCase();
  const k=Object.keys(MAT).find(k=>nm.includes(k));if(!k)return;const m=MAT[k],tr=m.opacity<1;
  n.geometry.computeVertexNormals();
  n.material=new THREE.MeshStandardMaterial({color:m.color,roughness:m.roughness,metalness:0,side:THREE.DoubleSide,
   transparent:tr,opacity:m.opacity,depthWrite:!tr,
   clippingPlanes:k==="peca"?[pCorte,pCam]:[pCorte]});
  n.renderOrder={peca:0,cera:1,recipiente:2,tampa:3}[k];P[k]=n;});
 const bx=new THREE.Box3().setFromObject(g.scene);cen=bx.getCenter(new THREE.Vector3());ext=bx;
 g.scene.position.sub(cen);s.add(g.scene);
 const sz=bx.getSize(new THREE.Vector3()).length();c.position.set(sz*.7,-sz*.95,sz*.75);ctl.update();aplica();});
const $=id=>document.getElementById(id);
function aplica(){if(!ext)return;
 for(const k of Object.keys(MAT))if(P[k])P[k].visible=$("v_"+k).checked;
 if(P.tampa)P.tampa.position.z=$("aberta").checked?ABRE:0;
 const f=+$("corte").value/1000,y=ext.min.y+(ext.max.y-ext.min.y)*f;pCorte.constant=(y-cen.y)+1e-3;
 const n=+$("cam").value;$("camv").textContent=n;pCam.constant=(FUNDO+n*CAMADA+1e-3)-cen.z;}
for(const id of ["v_peca","v_cera","v_recipiente","v_tampa","aberta","corte","cam"])$(id).oninput=$(id).onchange=aplica;
addEventListener("resize",()=>{r.setSize(innerWidth,innerHeight);c.aspect=innerWidth/innerHeight;c.updateProjectionMatrix();});
(function loop(){requestAnimationFrame(loop);ctl.update();r.render(s,c);})();
</script></body></html>"""


def cera(info: dict, info_r: dict, raiz: str) -> trimesh.Trimesh:
    """Prisma da cavidade (máscara dilatada pela folga), de z = fundo até fundo + H."""
    z = np.load(os.path.join(raiz, info["origem"], "intermediarios.npz"))
    v0, v = float(z["voxel"]), float(info_r["voxel_malha_mm"])
    f = max(1, int(round(v / v0)))
    mask = reamostrar_media(z["mask"].astype(np.float32), f) > 0.5
    folga = float(info_r["folga_mm"])
    if folga > 0:
        mask = ndi.binary_dilation(mask, morphology.disk(max(1, int(round(folga / v)))))
    fundo, H = float(info_r["fundo_mm"]), float(info_r["altura_cavidade_mm"])
    return solido_entre(np.full(mask.shape, fundo), np.full(mask.shape, fundo + H), mask, v * f / f)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    pasta = sys.argv[1]
    raiz = os.getcwd()
    info = json.load(open(os.path.join(pasta, "info.json"), encoding="utf-8"))
    info_r = json.load(open(os.path.join(pasta, "info_recipiente.json"), encoding="utf-8"))
    fundo, et, H = float(info_r["fundo_mm"]), float(info_r["tampa_mm"]), float(info_r["altura_cavidade_mm"])

    peca = trimesh.load(os.path.join(pasta, "peca_unica.stl"))
    peca.apply_translation([0, 0, fundo])                      # a peça apoia no fundo do recipiente
    rec = trimesh.load(os.path.join(pasta, "recipiente.stl"))
    tam = trimesh.load(os.path.join(pasta, "tampa.stl"))
    # tampa em uso: espelhada em z (aba para dentro da cavidade), placa de fundo+H até fundo+H+et
    espelho = np.diag([1.0, 1.0, -1.0, 1.0])
    espelho[2, 3] = fundo + H + et
    tam.apply_transform(espelho)
    cer = cera(info, info_r, raiz)

    cena = trimesh.Scene()
    for nome, m in (("peca", peca), ("cera", cer), ("recipiente", rec), ("tampa", tam)):
        m = reduzir(m, FACES[nome])
        cena.add_geometry(m, node_name=nome, geom_name=nome)
        print(nome, len(m.faces), "faces")
    buf = io.BytesIO(cena.export(file_type="glb"))
    b64 = base64.b64encode(buf.getvalue()).decode()

    nome = os.path.basename(os.path.normpath(pasta))
    linhas = [
        f"<b>Peça única vazada (PLA)</b>: {info['n_camadas']} camadas × {info['camada_mm']} mm = {info['altura_peca_mm']} mm · "
        f"{info['PLA_g']} g · {info['malha']['corpos']} corpo · vazios presos: {info['vazios_presos']}",
        f"erro de quantização {info['erro_medio_mm_PLA']} mm de PLA (placas de 1 mm: {info['erro_medio_placas_origem_mm_PLA']} mm) · "
        f"fração glandular {info['fracao_glandular_impressa'] * 100:.1f} % · fator de material {info['fator_material']}",
        f"<b>Recipiente e tampa (resina)</b>: parede {info_r['parede_mm']} mm · fundo {info_r['fundo_mm']} mm · "
        f"tampa {info_r['tampa_mm']} mm · cavidade {info_r['altura_cavidade_mm']:g} mm · "
        f"folga peça→parede ≥ {info_r['folga_minima_peca_parede_mm']} mm",
        f"<b>Cera de abelha</b>: ≈ {info_r['cera_mL']} mL · {info_r['cera_acima_da_peca_mm']} mm acima da peça",
        "<span style='opacity:.7'>visualização com malhas simplificadas; medidas nos STL originais</span>",
    ]
    html = (PAGINA.replace("__NOME__", nome).replace("__NC__", str(info["n_camadas"]))
            .replace("__INFO__", "<br>".join(linhas)).replace("__FUNDO__", f"{fundo}")
            .replace("__CAMADA__", f"{info['camada_mm']}").replace("__ABRE__", f"{0.8 * H:.1f}")
            .replace("__GLB__", b64))
    out = os.path.join(pasta, "mostruario.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print("ok:", os.path.abspath(out), f"{len(html) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
