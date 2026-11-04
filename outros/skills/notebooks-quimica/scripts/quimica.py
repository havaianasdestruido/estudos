"""Funções prontas para os notebooks de química deste repositório.

Copie a função que precisar para dentro do notebook (o material didático deve
ser autocontido). Este arquivo é a fonte única: teste aqui antes de colar.

Dependências: rdkit, py3Dmol, matplotlib (opcional: plotly).
"""

from __future__ import annotations

from typing import Iterable, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

try:  # legenda sem acento: o desenho do RDKit corrompe "á", "ç", "₂", ver rotulos.py
    from rotulos import galeria_com_legendas, figura_com_legenda, legenda_segura
except ImportError:  # notebook autocontido: cole a função lá do SKILL.md
    def legenda_segura(texto: str) -> str:  # type: ignore[misc]
        import unicodedata
        tabela = str.maketrans("₀₁₂₃₄₅₆₇₈₉°", "0123456789o")
        texto = unicodedata.normalize("NFKD", texto.translate(tabela))
        return texto.encode("ascii", "ignore").decode("ascii")

rdDepictor.SetPreferCoordGen(True)

_SEED = 0xC0FFEE


# ---------------------------------------------------------------------------
# 3D interativo
# ---------------------------------------------------------------------------

def mol3d(smiles: str, seed: int = _SEED) -> Chem.Mol:
    """SMILES -> mol com hidrogênios, geometria 3D reproduzível e otimizada.

    Levanta ValueError se o SMILES for inválido e RuntimeError se a geração da
    geometria falhar (acontece em moléculas muito rígidas ou com metais).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    mol = Chem.AddHs(mol)  # sem H a geometria 3D fica errada e a figura fica feia

    params = AllChem.ETKDGv3()
    params.randomSeed = seed  # mesmo desenho em toda execução
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"não foi possível gerar a geometria 3D de {smiles!r}")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)  # relaxa a estrutura
    else:  # alguns íons não têm parâmetros MMFF; UFF cobre o resto
        AllChem.UFFOptimizeMolecule(mol)
    return mol


ESTILOS_3D = {
    # nome: spec do 3Dmol.js
    "ball-and-stick": {"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                       "sphere": {"scale": 0.28, "colorscheme": "Jmol"}},
    "stick": {"stick": {"radius": 0.15, "colorscheme": "Jmol"}},
    "space-filling": {"sphere": {"scale": 0.95, "colorscheme": "Jmol"}},
    "wireframe": {"line": {"linewidth": 2.0, "colorscheme": "Jmol"}},
}


def molecula_3d(smiles: str, estilo: str = "ball-and-stick", largura: int = 640,
                altura: int = 420, rodar: bool = False, fundo: str = "#ffffff",
                rotulo: str | None = None):
    """Molécula 3D colorida e INTERATIVA (arraste para girar 360°, rolo = zoom).

    estilo: "ball-and-stick" (padrão), "stick", "space-filling" ou "wireframe".
    rodar: True liga o giro automático (use só em uma figura de destaque).
    """
    import py3Dmol  # importado aqui para o erro aparecer só quando for usar

    if estilo not in ESTILOS_3D:
        raise ValueError(f"estilo deve ser um de {sorted(ESTILOS_3D)}")

    mol = mol3d(smiles)
    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle(ESTILOS_3D[estilo])
    visor.setBackgroundColor(fundo)
    if rodar:
        visor.spin("y", 0.5)
    if rotulo:
        visor.addLabel(rotulo, {"fontColor": "black", "backgroundColor": "white",
                                "backgroundOpacity": 0.75, "fontSize": 13})
    visor.zoomTo()
    return visor


def mostrar_molecula(smiles: str, estilo: str = "ball-and-stick", largura: int = 640,
                     altura: int = 420, rotulo: str | None = None):
    """Mostra a molécula: 3D interativo se der; se não, cai para o 2D.

    Útil quando o notebook vai rodar em máquina sem internet (o 3Dmol.js vem de
    CDN) ou quando alguém abriu sem instalar o py3Dmol.
    """
    try:
        return molecula_3d(smiles, estilo=estilo, largura=largura, altura=altura,
                           rotulo=rotulo)
    except ImportError:
        print("py3Dmol não instalado — mostrando a estrutura 2D no lugar do 3D.")
        return estrutura_2d(smiles, legenda=rotulo, largura=largura, altura=altura)


def geometria_3d(smiles: str, atomo: int = 0, largura: int = 640, altura: int = 460,
                 cor_vetor: str = "#1f6feb"):
    """3D interativo com setas do átomo central aos vizinhos e o ângulo medido.

    Devolve (visor, angulo_em_graus) usando os dois primeiros vizinhos.
    """
    import numpy as np
    import py3Dmol

    mol = mol3d(smiles)
    conf = mol.GetConformer()
    centro = np.array(conf.GetAtomPosition(atomo))
    vizinhos = [n.GetIdx() for n in mol.GetAtomWithIdx(atomo).GetNeighbors()]
    if len(vizinhos) < 2:
        raise ValueError("o átomo escolhido precisa ter pelo menos 2 vizinhos")

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle(ESTILOS_3D["ball-and-stick"])
    visor.setBackgroundColor("#ffffff")

    for i in vizinhos:
        p = np.array(conf.GetAtomPosition(i))
        visor.addArrow({"start": {"x": float(centro[0]), "y": float(centro[1]),
                                  "z": float(centro[2])},
                        "end": {"x": float(p[0]), "y": float(p[1]), "z": float(p[2])},
                        "radius": 0.045, "radiusRatio": 1.8, "mid": 0.85,
                        "color": cor_vetor})

    a, b = vizinhos[0], vizinhos[1]
    va = np.array(conf.GetAtomPosition(a)) - centro
    vb = np.array(conf.GetAtomPosition(b)) - centro
    angulo = float(np.degrees(np.arccos(
        np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))))
    meio = centro + 0.6 * (va + vb) / 2
    visor.addLabel(f"{angulo:.1f}°", {"position": {"x": float(meio[0]),
                                                   "y": float(meio[1]),
                                                   "z": float(meio[2])},
                                      "fontColor": cor_vetor, "fontSize": 14,
                                      "backgroundOpacity": 0.0})
    visor.zoomTo()
    return visor, angulo


# ---------------------------------------------------------------------------
# 2D "com risquinhos"
# ---------------------------------------------------------------------------

def desenho_2d(mol: Chem.Mol, largura: int = 520, altura: int = 360,
               legenda: str = "", espessura: float = 2.2, legenda_livre: bool = False):
    """Mol -> SVG do RDKit (traço preto, duplas com 2 riscos, etc.).

    A legenda é convertida para ASCII por padrão (o desenho do RDKit corrompe
    acentos/subscritos). Para título bonito, use `figura_com_legenda`.
    """
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    opts = desenho.drawOptions()
    opts.bondLineWidth = espessura
    opts.addStereoAnnotation = True   # mostra R/S quando houver carbono quiral
    opts.padding = 0.08
    desenho.DrawMolecule(mol, legend=legenda if legenda_livre else legenda_segura(legenda))
    desenho.FinishDrawing()
    from IPython.display import SVG
    return SVG(desenho.GetDrawingText())


def estrutura_2d(smiles: str, legenda: str | None = None, mostrar_h: bool = False,
                 largura: int = 520, altura: int = 360):
    """Estrutura 2D de uma molécula a partir do SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if mostrar_h:
        mol = Chem.AddHs(mol)
    return desenho_2d(mol, largura=largura, altura=altura,
                      legenda=legenda if legenda is not None else formula(smiles))


def galeria_2d(smiles: Sequence[str], legendas: Iterable[str] | None = None,
               por_linha: int = 3, tamanho: tuple[int, int] = (300, 240)):
    """Várias estruturas 2D numa grade (para comparar funções orgânicas)."""
    mols = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise ValueError(f"SMILES inválido: {s!r}")
        mols.append(mol)
    from IPython.display import SVG
    textos = [legenda_segura(t) for t in legendas] if legendas else None
    imagem = Draw.MolsToGridImage(mols, molsPerRow=por_linha, subImgSize=tamanho,
                                  legends=textos, useSVG=True)
    # com IPython ativo o RDKit devolve um objeto SVG; em script puro, uma string
    return imagem if hasattr(imagem, "data") else SVG(imagem)


def estrutura_realcada(smiles: str, padrao_smarts: str, largura: int = 520,
                       altura: int = 340, cor: tuple[float, float, float] = (0.99, 0.75, 0.30)):
    """Estrutura 2D com um grupo funcional destacado (SMARTS)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    padrao = Chem.MolFromSmarts(padrao_smarts)
    if padrao is None:
        raise ValueError(f"SMARTS inválido: {padrao_smarts!r}")
    correspondencia = mol.GetSubstructMatch(padrao)
    if not correspondencia:
        raise ValueError(f"{smiles!r} não contém o padrão {padrao_smarts!r}")

    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.DrawMolecule(mol, highlightAtoms=list(correspondencia),
                         highlightAtomColors={i: cor for i in correspondencia})
    desenho.FinishDrawing()
    from IPython.display import SVG
    return SVG(desenho.GetDrawingText())


def ligacoes_hidrogenio(pares: Sequence[tuple[int, int]]):
    """Decorador simples: recebe um mol e devolve função que liga pares com
    ligação de ordem zero (o RDKit desenha pontilhado = ligação de hidrogênio).

    Exemplo:
        mol = Chem.AddHs(Chem.MolFromSmiles("O.O"))
        mol = ligacoes_hidrogenio([(1, 3)])(mol)
    """
    def aplicar(mol: Chem.Mol) -> Chem.Mol:
        rw = Chem.RWMol(mol)
        for i, j in pares:
            if i >= rw.GetNumAtoms() or j >= rw.GetNumAtoms():
                raise IndexError(f"átomo fora da molécula: par ({i}, {j})")
            rw.AddBond(int(i), int(j), Chem.BondType.ZERO)
        novo = rw.GetMol()
        novo.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(novo, Chem.SanitizeFlags.SANITIZE_ALL
                         ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return novo
    return aplicar


# ---------------------------------------------------------------------------
# Dados / utilidades
# ---------------------------------------------------------------------------

def formula(smiles: str) -> str:
    """Fórmula molecular (Hill) — ex.: 'C2H4O2'."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    return rdMolDescriptors.CalcMolFormula(mol)


def massa_molar(smiles: str) -> float:
    """Massa molar em g/mol."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    return Descriptors.MolWt(mol)


def descricao(smiles: str) -> dict:
    """Resumo rápido para a célula de 'conferindo o cálculo'."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    return {
        "smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "massa_molar_g_mol": round(Descriptors.MolWt(mol), 2),
        "atomos_pesados": mol.GetNumHeavyAtoms(),
        "ligacoes": mol.GetNumBonds(),
        "anéis": Chem.GetSSSR(mol) and len(Chem.GetSSSR(mol)) or 0,
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "logp": round(Descriptors.MolLogP(mol), 2),
    }


def eh_valido(smiles: str) -> bool:
    """Teste rápido de SMILES (use no exercício: 'o aluno digitou certo?')."""
    return Chem.MolFromSmiles(smiles) is not None
