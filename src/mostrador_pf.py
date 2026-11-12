"""Mostrador lado a lado: limiar fibroso 0,6 × 0,8 × 1,0 (pastas <img>_pf06, _pf08 e _pf10).

Gera, para cada imagem, um HTML autônomo (GLB de cada modelo embutido) com câmera
sincronizada, controle "até a placa N", ligar/desligar a gordura (recipiente) e os números
de cada variante; e um index.html com todos os casos. Variantes ausentes são ignoradas
(o mostrador precisa de ao menos duas).

Uso: python src/mostrador_pf.py [saida]   (procura <grupo>/<nome>_pfXX; grava em saida/mostrador_pf)
"""
from __future__ import annotations

import base64
import glob
import json
import os
import sys

VARIANTES = [("pf06", "0,6", "mais tênues incluídos"),
             ("pf08", "0,8", "só o mais denso"),
             ("pf10", "1,0", "padrão do gerador")]

PAGINA = """<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>__NOME__ · limiar __TIT__</title>
<style>
html,body{margin:0;height:100%;background:#111;color:#ddd;font:13px system-ui;overflow:hidden}
#barra{position:absolute;left:0;right:0;top:0;padding:8px 12px;background:#000a;display:flex;gap:18px;
 align-items:center;flex-wrap:wrap;z-index:2}
#barra b{color:#fff}
.rot{position:absolute;top:44px;z-index:2;background:#000a;padding:6px 10px;border-radius:6px;pointer-events:none;
 line-height:1.5;font-size:12px}
.sep{position:absolute;top:0;bottom:0;width:1px;background:#555;z-index:1}
input[type=range]{width:180px;vertical-align:middle}
a{color:#8cf}
</style>
<script type="importmap">{"imports":{
 "three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
 "three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}</script></head>
<body>
<div id="barra"><b>__NOME__</b>
 <span>até a placa <input id="n" type="range" min="1" max="__NP__" value="__NP__"> <span id="nv">__NP__</span></span>
 <label><input id="gord" type="checkbox" checked> gordura (parafina)</label>
 <span><a href="index.html">← todas as imagens</a></span>
 <span style="opacity:.7">arraste = girar · roda = zoom · botão direito = mover (câmera sincronizada)</span></div>
__ROTULOS__
<script type="module">
import * as THREE from "three";
import {OrbitControls} from "three/addons/controls/OrbitControls.js";
import {GLTFLoader} from "three/addons/loaders/GLTFLoader.js";
const B=[__BS__],K=B.length,cenas=[],plc=B.map(()=>[]),gord=B.map(()=>[]);let ref=null;
const r=new THREE.WebGLRenderer({antialias:true});r.setPixelRatio(devicePixelRatio);
r.setSize(innerWidth,innerHeight);r.setScissorTest(true);document.body.appendChild(r.domElement);
const c=new THREE.PerspectiveCamera(40,innerWidth/K/innerHeight,0.1,5000);
const ctl=new OrbitControls(c,r.domElement);ctl.enableDamping=true;
function nova(){const s=new THREE.Scene();s.background=new THREE.Color(0x111111);
 s.add(new THREE.HemisphereLight(0xffffff,0x333344,1.6));
 const d=new THREE.DirectionalLight(0xffffff,1.4);d.position.set(1,1,2);s.add(d);return s;}
B.forEach((b64,i)=>{const s=nova();cenas.push(s);
 const bin=Uint8Array.from(atob(b64),ch=>ch.charCodeAt(0)).buffer;
 new GLTFLoader().parse(bin,"",g=>{const o=g.scene;o.traverse(n=>{if(!n.isMesh)return;
  n.geometry.computeVertexNormals();const nm=(n.name+" "+(n.parent?n.parent.name:"")).toLowerCase();
  const adip=nm.includes("adiposo");
  n.material=new THREE.MeshStandardMaterial({vertexColors:!!n.geometry.attributes.color,
   color:n.geometry.attributes.color?0xffffff:0xd9b8a8,roughness:.6,side:THREE.DoubleSide,
   transparent:adip,opacity:adip?0.18:1,depthWrite:!adip});
  if(adip){n.renderOrder=1;gord[i].push(n);}
  else{const m=nm.match(/placa_(\\d+)/);if(m)plc[i][parseInt(m[1])]=n;}});
  if(!ref){const bx=new THREE.Box3().setFromObject(o);ref=bx.getCenter(new THREE.Vector3());
   const sz=bx.getSize(new THREE.Vector3()).length();
   c.position.set(sz*.35,-sz*.9,sz*.7);c.up.set(0,0,1);ctl.update();}
  o.position.sub(ref);s.add(o);});});
const n=document.getElementById("n"),nv=document.getElementById("nv"),gc=document.getElementById("gord");
function aplica(){nv.textContent=n.value;for(let i=0;i<K;i++){
 plc[i].forEach((m,k)=>{if(m)m.visible=k<=+n.value;});gord[i].forEach(m=>m.visible=gc.checked);}}
n.oninput=gc.onchange=aplica;
addEventListener("resize",()=>{r.setSize(innerWidth,innerHeight);c.aspect=innerWidth/K/innerHeight;c.updateProjectionMatrix();});
(function loop(){requestAnimationFrame(loop);ctl.update();aplica();
 const w=innerWidth/K,h=innerHeight;
 for(let i=0;i<K;i++){r.setViewport(i*w,0,w,h);r.setScissor(i*w,0,w,h);r.render(cenas[i],c);}})();
</script></body></html>"""


def b64(caminho: str) -> str:
    with open(caminho, "rb") as f:
        return base64.b64encode(f.read()).decode()


def rotulo(titulo: str, info: dict, esq_pct: float, larg_pct: float) -> str:
    return (f"<div class='rot' style='left:calc({esq_pct:.3f}% + 8px);max-width:calc({larg_pct:.3f}% - 24px)'>"
            f"<b>{titulo}</b><br>{info['n_placas']} placas × {info['espessura_placa_mm']} mm · "
            f"fração glandular {info['fracao_glandular_impressa'] * 100:.1f} %<br>"
            f"limiar gordura×fibroso {info['limiar_gordura_glandula']:.1f} · "
            f"pilha {info['altura_pilha_mm']} mm · parafina ≈ {info['parafina_estimada_mL']} mL</div>")


def main():
    raiz = sys.argv[1] if len(sys.argv) > 1 else "saida"
    dest = os.path.join(raiz, "mostrador_pf")
    os.makedirs(dest, exist_ok=True)
    # layout: <raiz>/<grupo>/<nome>_pfXX/<nome>/modelo.glb
    casos: dict[str, dict[str, str]] = {}
    for sufixo, _, _ in VARIANTES:
        for g in glob.glob(os.path.join(raiz, "*", f"*_{sufixo}", "*", "modelo.glb")):
            p = os.path.dirname(g)
            casos.setdefault(os.path.basename(p), {})[sufixo] = p
    linhas = []
    for nome, pastas in sorted(casos.items()):
        vs = [(s, v, d) for s, v, d in VARIANTES if s in pastas]
        if len(vs) < 2:
            print("pulado (menos de 2 variantes):", nome)
            continue
        infos = {s: json.load(open(os.path.join(pastas[s], "info.json"), encoding="utf-8")) for s, _, _ in vs}
        k = len(vs)
        rotulos = "".join(rotulo(f"limiar {v} ({d})", infos[s], 100 * i / k, 100 / k)
                          for i, (s, v, d) in enumerate(vs))
        seps = "".join(f"<div class='sep' style='left:{100 * i / k:.3f}%'></div>" for i in range(1, k))
        html = (PAGINA.replace("__NOME__", nome).replace("__TIT__", " × ".join(v for _, v, _ in vs))
                .replace("__NP__", str(max(i["n_placas"] for i in infos.values())))
                .replace("__ROTULOS__", rotulos + seps)
                .replace("__BS__", ",".join('"' + b64(os.path.join(pastas[s], "modelo.glb")) + '"' for s, _, _ in vs)))
        with open(os.path.join(dest, f"{nome}.html"), "w", encoding="utf-8") as f:
            f.write(html)
        cel = "".join(f"<td>{infos[s]['fracao_glandular_impressa'] * 100:.1f} % · {infos[s]['n_placas']} pl.</td>"
                      if s in infos else "<td>—</td>" for s, _, _ in VARIANTES)
        linhas.append(f"<tr><td><a href='{nome}.html'>{nome}</a></td>{cel}</tr>")
        print("ok", nome, [v for _, v, _ in vs])
    cab = "".join(f"<th>limiar {v}</th>" for _, v, _ in VARIANTES)
    idx = ("<!doctype html><meta charset='utf-8'><title>Mostrador do limiar fibroso</title>"
           "<style>body{font:15px system-ui;background:#111;color:#ddd;margin:24px}a{color:#8cf}"
           "td,th{padding:6px 16px;text-align:left;border-bottom:1px solid #333}</style>"
           "<h2>Limiar fibroso 0,6 × 0,8 × 1,0</h2><p>fração glandular impressa · nº de placas</p>"
           f"<table><tr><th>imagem</th>{cab}</tr>" + "".join(linhas) + "</table>")
    with open(os.path.join(dest, "index.html"), "w", encoding="utf-8") as f:
        f.write(idx)
    print("em:", os.path.abspath(dest))


if __name__ == "__main__":
    main()
