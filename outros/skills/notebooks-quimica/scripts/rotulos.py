"""Rótulos com acento/subscrito para desenhos 2D do RDKit.

Por que este arquivo existe
---------------------------
O desenho do RDKit é montado em SVG/Cairo com a API de fonte do sistema. Nesta
combinação (RDKit 2026.09 + Cairo 1.18), caracteres fora do ASCII **saem como
quadradinhos** em legendas e rótulos: "ácido acético" vira "cido ac tico",
"H₂O" vira "H" + caixinhas. Isso é feio em material didático e acontece no
Colab também.

Solução adotada neste repositório (duas camadas):

1. `legenda_segura(texto)` — deixa a legenda em ASCII para o desenho do RDKit:
   remove acentos, troca subscrito por dígito normal. Assim a imagem sempre fica
   legível.
2. `figura_com_legenda(...)` — quando o rótulo precisa sair bonito (com acento e
   subscrito de verdade), desenha a molécula **sem** legenda e coloca o título
   com **matplotlib**, que renderiza Unicode e LaTeX mathtext corretamente.

Use a camada 1 dentro de loops/galerias rápidas e a camada 2 em toda figura que
vai para o material final.
"""

from __future__ import annotations

import io
import unicodedata

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

_SUBSCRITOS = str.maketrans("₀₁₂₃₄₅₆₇₈₉₊₋", "0123456789+-")
_SOBRESCRITOS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")


def legenda_segura(texto: str) -> str:
    """Deixa o texto legível no motor de desenho do RDKit (ASCII, sem acento)."""
    texto = texto.translate(_SUBSCRITOS).translate(_SOBRESCRITOS)
    texto = texto.replace("°", " graus").replace("º", "o").replace("ª", "a")
    # remove acentos: "ácido" -> "acido", "ligação" -> "ligacao"
    texto = unicodedata.normalize("NFKD", texto)
    return texto.encode("ascii", "ignore").decode("ascii")


def desenho_svg(mol: Chem.Mol, largura: int = 520, altura: int = 360,
                legenda: str = "", espessura: float = 2.2, seguro: bool = True):
    """Mol -> SVG. `seguro=True` aplica `legenda_segura` (recomendado)."""
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    opcoes = desenho.drawOptions()
    opcoes.bondLineWidth = espessura
    opcoes.addStereoAnnotation = True
    opcoes.padding = 0.08
    texto = legenda_segura(legenda) if seguro else legenda
    desenho.DrawMolecule(mol, legend=texto)
    desenho.FinishDrawing()
    from IPython.display import SVG
    return SVG(desenho.GetDrawingText())


def figura_com_legenda(smiles: str, titulo: str, largura: int = 520,
                       altura: int = 360, mostrar_h: bool = False,
                       tamanho_figura: tuple[float, float] | None = None):
    """Figura com a molécula + título em matplotlib (acentos e subscritos ok).

    O título aceita mathtext: use `"Ácido acético (CH$_3$COOH)"` para subscrito.
    Devolve a `matplotlib.figure.Figure` (o notebook exibe automaticamente).
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if mostrar_h:
        mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol)

    desenho = rdMolDraw2D.MolDraw2DCairo(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.drawOptions().addStereoAnnotation = True
    desenho.DrawMolecule(mol)
    desenho.FinishDrawing()
    imagem = np.asarray(Image.open(io.BytesIO(desenho.GetDrawingText())))

    # 100 px por polegada em dpi padrão: a figura não fica borrada
    largura_in, altura_in = tamanho_figura or (largura / 100, (altura + 46) / 100)
    figura, eixo = plt.subplots(figsize=(largura_in, altura_in))
    eixo.imshow(imagem)
    eixo.set_axis_off()
    eixo.set_title(titulo, fontsize=11, pad=8)
    figura.tight_layout()
    return figura


def galeria_com_legendas(smiles, titulos, por_linha: int = 3,
                         tamanho_celula: tuple[int, int] = (320, 260)):
    """Grade de estruturas 2D com títulos em matplotlib (acentos corretos)."""
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image

    linhas = (len(smiles) + por_linha - 1) // por_linha
    figura, eixos = plt.subplots(linhas, por_linha,
                                 figsize=(3.4 * por_linha, 2.9 * linhas))
    eixos = np.atleast_1d(eixos).ravel()

    for eixo in eixos[len(smiles):]:
        eixo.set_axis_off()

    for eixo, smi, titulo in zip(eixos, smiles, titulos):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"SMILES inválido: {smi!r}")
        rdDepictor.Compute2DCoords(mol)
        desenho = rdMolDraw2D.MolDraw2DCairo(*tamanho_celula)
        desenho.drawOptions().bondLineWidth = 2.0
        desenho.drawOptions().padding = 0.12   # nada de "-O" colado na borda
        desenho.DrawMolecule(mol)
        desenho.FinishDrawing()
        eixo.imshow(np.asarray(Image.open(io.BytesIO(desenho.GetDrawingText()))))
        eixo.set_axis_off()
        eixo.set_title(titulo, fontsize=10, pad=6)

    figura.tight_layout()
    return figura
