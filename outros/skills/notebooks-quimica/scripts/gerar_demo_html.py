"""Gera uma página HTML autocontida com moléculas 3D interativas.

Serve para (a) o aluno ver o resultado sem abrir o notebook e (b) você conferir
o visual antes de escrever as células do notebook.

Uso:
    python gerar_demo_html.py ../exemplos/demo_3d.html
"""

from __future__ import annotations

import argparse
import html as html_lib
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

MOLECULAS = [
    ("CCO", "Etanol", "C₂H₅OH — ligação O–H polar"),
    ("CC(=O)O", "Ácido acético", "CH₃COOH — grupo carboxila"),
    ("c1ccccc1", "Benzeno", "C₆H₆ — anel aromático plano"),
    ("O", "Água", "H₂O — angular, 104,5°"),
]


def molblock(smiles: str) -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parametros = AllChem.ETKDGv3()
    parametros.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, parametros) != 0:
        raise RuntimeError(f"não foi possível gerar a geometria 3D de {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)


def visor_html(smiles: str, largura: int = 420, altura: int = 360) -> str:
    """HTML do py3Dmol para uma molécula (embute o molblock, não baixa nada)."""
    import py3Dmol

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(molblock(smiles), "mol")
    visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.28, "colorscheme": "Jmol"}})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    visor.spin("y", 0.4)   # giro lento: mostra de cara que é 3D e arrastável
    return visor._make_html()


def montar_pagina() -> str:
    cartoes = []
    for smiles, nome, descricao in MOLECULAS:
        formula = rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))
        cartoes.append(f"""
      <figure class="cartao">
        {visor_html(smiles)}
        <figcaption>
          <strong>{html_lib.escape(nome)}</strong>
          <span class="formula">{html_lib.escape(formula)}</span>
          <br><span class="dica">{html_lib.escape(descricao)}</span>
        </figcaption>
      </figure>""")
    return f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<title>Moléculas em 3D — arraste para girar</title>
<style>
  body {{ font-family: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
         margin: 0; padding: 32px 24px 56px; background: #f6f7f9; color: #1c1f23; }}
  header {{ max-width: 960px; margin: 0 auto 28px; }}
  h1 {{ font-size: 1.6rem; margin: 0 0 8px; }}
  p.instrucao {{ margin: 0; color: #47505a; }}
  .galeria {{ max-width: 960px; margin: 0 auto; display: grid;
             grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 20px; }}
  .cartao {{ margin: 0; background: #fff; border: 1px solid #e3e6ea; border-radius: 12px;
            padding: 12px; box-shadow: 0 1px 2px rgba(16, 24, 40, .04); }}
  figcaption {{ padding: 8px 6px 2px; font-size: .92rem; line-height: 1.5; }}
  .formula {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
             color: #6b7480; margin-left: 6px; }}
  .dica {{ color: #47505a; }}
  .aviso {{ max-width: 960px; margin: 0 auto 20px; padding: 12px 16px;
            background: #fff8e6; border: 1px solid #f2dfae; border-radius: 10px;
            font-size: .9rem; color: #5c4a1c; }}
  footer {{ max-width: 960px; margin: 32px auto 0; color: #6b7480; font-size: .85rem; }}
</style>
</head>
<body>
  <header>
    <h1>Moléculas em 3D</h1>
    <p class="instrucao">
      <strong>Arraste com o mouse</strong> para girar em 360° em qualquer direção.
      Use a roda do mouse (ou dois dedos) para dar zoom, e o botão direito para mover.
      As cores seguem a convenção dos químicos: carbono cinza, hidrogênio branco,
      oxigênio vermelho.
    </p>
  </header>
  <div class="aviso">
    <strong>Se as moléculas não aparecerem:</strong> a biblioteca 3Dmol.js vem de
    um CDN e esta página precisa de internet. No Jupyter/Colab acontece o mesmo.
  </div>
  <main class="galeria">{''.join(cartoes)}
  </main>
  <footer>
    Gerado por <code>outros/skills/notebooks-quimica/scripts/gerar_demo_html.py</code>
    (RDKit + 3Dmol.js). Para funcionar, a página precisa de internet: o 3Dmol.js vem de CDN.
  </footer>
</body>
</html>
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("saida", type=Path, nargs="?", default=Path("demo_3d.html"))
    args = parser.parse_args(argv)
    args.saida.parent.mkdir(parents=True, exist_ok=True)
    args.saida.write_text(montar_pagina(), encoding="utf-8")
    print(f"página escrita em {args.saida}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
