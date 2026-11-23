# -*- coding: utf-8 -*-
"""Rascunho-fonte de quimica/bases/01_bases_hidroxilas_amonia_poh.ipynb.

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/fontes/bases.py quimica/bases/01_bases_hidroxilas_amonia_poh.ipynb
"""

# %% md
# Bases: hidroxilas, amônia e pOH
### Funções inorgânicas: as bases

**Objetivo:** Reconhecer a função base (NaOH, KOH, Ca(OH)$_2$, Mg(OH)$_2$, NH$_3$, aminas), escrever a dissociação em água, separar bases fortes de fracas pelo $K_b$ e calcular o pOH e o pH de uma solução.

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. As cores seguem a convenção dos químicos (Na = roxo, K = violeta, Ca = verde-escuro, Mg = verde, N = azul, O = vermelho, H = branco).

# %% code title="Preparar o ambiente (só precisa rodar uma vez)"
import importlib.util, subprocess, sys

def _tem(modulo):
    return importlib.util.find_spec(modulo) is not None

for pacote in ("rdkit", "py3Dmol"):
    if not _tem(pacote):
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", pacote], check=True)

import rdkit, py3Dmol
print("rdkit", rdkit.__version__, "| py3Dmol", py3Dmol.__version__)

# %% code title="Funções auxiliares de desenho (3D interativo e 2D com risquinhos)"
import io, unicodedata
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from IPython.display import SVG
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")        # silencia os avisos de campo de força dos íons
rdDepictor.SetPreferCoordGen(True)
SEMENTE = 0xC0FFEE

ESTILOS_3D = {
    "ball-and-stick": {"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                       "sphere": {"scale": 0.28, "colorscheme": "Jmol"}},
    "stick": {"stick": {"radius": 0.15, "colorscheme": "Jmol"}},
    "space-filling": {"sphere": {"scale": 0.95, "colorscheme": "Jmol"}},
}

def mol3d(smiles, semente=SEMENTE, campo="auto"):
    """SMILES -> molécula com H explícitos, geometria 3D reproduzível e otimizada.

    campo: "auto" (MMFF se existir, senão UFF), "MMFF", "UFF" ou "nenhum".
    Nos óxidos de enxofre o MMFF erra os ângulos; ali use campo="UFF".
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    mol = Chem.AddHs(mol)
    parametros = AllChem.ETKDGv3()
    parametros.randomSeed = semente          # mesmo desenho em toda execução
    if AllChem.EmbedMolecule(mol, parametros) != 0:
        raise RuntimeError(f"não foi possível gerar a geometria 3D de {smiles!r}")
    if campo == "auto":
        campo = "MMFF" if AllChem.MMFFHasAllMoleculeParams(mol) else "UFF"
    if campo == "MMFF" and AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    elif campo == "UFF" and AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol)
    # íons metálicos não têm parâmetros de campo de força: fica só a geometria
    return mol

def molecula_3d(smiles, largura=620, altura=400, estilo="ball-and-stick", rotulo=None,
                campo="auto"):
    """Molécula 3D colorida e INTERATIVA: arraste para girar 360°, roda = zoom."""
    mol = mol3d(smiles, campo=campo)
    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle(ESTILOS_3D[estilo])
    visor.setBackgroundColor("#ffffff")
    if rotulo:
        visor.addLabel(rotulo, {"fontColor": "black", "backgroundColor": "white",
                                "backgroundOpacity": 0.75, "fontSize": 13})
    visor.zoomTo()
    return visor

def geometria_3d(smiles, atomo=0, largura=620, altura=420, cor="#1f6feb", campo="auto"):
    """3D interativo com setas do átomo central aos vizinhos e o ângulo medido."""
    mol = mol3d(smiles, campo=campo)
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
                        "radius": 0.045, "radiusRatio": 1.8, "mid": 0.85, "color": cor})
    a, b = vizinhos[0], vizinhos[1]
    va = np.array(conf.GetAtomPosition(a)) - centro
    vb = np.array(conf.GetAtomPosition(b)) - centro
    angulo = float(np.degrees(np.arccos(np.dot(va, vb) /
                                        (np.linalg.norm(va) * np.linalg.norm(vb)))))
    meio = centro + 0.6 * (va + vb) / 2
    visor.addLabel(f"{angulo:.1f}°", {"position": {"x": float(meio[0]),
                                                   "y": float(meio[1]),
                                                   "z": float(meio[2])},
                                      "fontColor": cor, "fontSize": 14,
                                      "backgroundOpacity": 0.0})
    visor.zoomTo()
    return visor, angulo

def legenda_segura(texto):
    """ASCII para o desenho do RDKit (ele corrompe acento e subscrito)."""
    tabela = str.maketrans("₀₁₂₃₄₅₆₇₈₉₊₋°", "0123456789+-o")
    texto = unicodedata.normalize("NFKD", texto.translate(tabela))
    return texto.encode("ascii", "ignore").decode("ascii")

def estrutura_2d(smiles, legenda="", largura=520, altura=360, mostrar_h=False):
    """Estrutura 2D com os risquinhos das ligações (1, 2 ou 3 traços)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if mostrar_h:
        mol = Chem.AddHs(mol)
    for atomo in mol.GetAtoms():        # idem: "C" só quando não tem C vizinho
        if atomo.GetSymbol() == "C" and not any(
                vizinho.GetSymbol() == "C" for vizinho in atomo.GetNeighbors()):
            atomo.SetProp("atomLabel", "C")
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    opcoes = desenho.drawOptions()
    opcoes.bondLineWidth = 2.2
    opcoes.padding = 0.10
    desenho.DrawMolecule(mol, legend=legenda_segura(legenda))
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

def ligacao_ionica(mol):
    """Liga os fragmentos iônicos com uma ligação de ordem ZERO (pontilhado).

    Sal e óxido metálico não têm ligação covalente entre cátion e ânion: existe
    atração eletrostática. O RDKit não desenha isso sozinho e, sem nada, os íons
    ficam "flutuando"; a ligação de ordem zero sai pontilhada, que é a convenção
    usada aqui (mesma das ligações de hidrogênio) para dizer "isto é atração,
    não covalência".
    """
    fragmentos = Chem.GetMolFrags(mol)
    if len(fragmentos) < 2:
        return mol
    cargas = {i: sum(mol.GetAtomWithIdx(a).GetFormalCharge() for a in fr)
              for i, fr in enumerate(fragmentos)}
    cation = max(cargas, key=cargas.get)
    positivos = [a for a in fragmentos[cation] if mol.GetAtomWithIdx(a).GetFormalCharge() > 0]
    ancora_cation = positivos[0] if positivos else fragmentos[cation][0]

    rw = Chem.RWMol(mol)
    for i, fragmento in enumerate(fragmentos):
        if i == cation:
            continue
        negativos = [a for a in fragmento if mol.GetAtomWithIdx(a).GetFormalCharge() < 0]
        ancora = negativos[0] if negativos else fragmento[0]
        rw.AddBond(int(ancora_cation), int(ancora), Chem.BondType.ZERO)
    novo = rw.GetMol()
    novo.UpdatePropertyCache(strict=False)
    for atomo in novo.GetAtoms():        # a ligação extra não deve virar "radical" (ponto)
        atomo.SetNumRadicalElectrons(0)
    return novo

def desenho_png(smiles, tamanho=(420, 340), mostrar_h=False, ionico=False,
                letra_carbono=True):
    """Molécula -> imagem (array) para colocar dentro de uma figura matplotlib.

    letra_carbono=True mostra a letra C (o padrão do RDKit esconde o carbono do
    esqueleto, e em CO2/CO/esqueleto de sal o desenho fica "O=O", errado).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if mostrar_h:
        mol = Chem.AddHs(mol)
    if letra_carbono:
        # mostra o "C" só quando ele é o único carbono em volta (CO2, CO, CO3(2-),
        # HCO3-): num esqueleto normal (CH3COOH) a convenção é não escrever o C
        for atomo in mol.GetAtoms():
            if atomo.GetSymbol() == "C" and not any(
                    vizinho.GetSymbol() == "C" for vizinho in atomo.GetNeighbors()):
                atomo.SetProp("atomLabel", "C")
    if ionico:
        mol = ligacao_ionica(mol)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DCairo(*tamanho)
    opcoes = desenho.drawOptions()
    opcoes.bondLineWidth = 2.4
    opcoes.padding = 0.06
    # dá um respiro em volta do rótulo: sem isso o traço da ligação atravessa a letra
    opcoes.additionalAtomLabelPadding = 0.09
    desenho.DrawMolecule(mol)
    desenho.FinishDrawing()
    return np.asarray(Image.open(io.BytesIO(desenho.GetDrawingText())))

def galeria_2d(smiles, titulos, por_linha=3, ionico=False, mostrar_h=False):
    """Grade de estruturas 2D com títulos em matplotlib (acento e subscrito certos).

    ionico=True liga cátion e ânion com pontilhado (ver `ligacao_ionica`).
    mostrar_h=True deixa os hidrogênios visíveis — nos hidróxidos isso mostra
    o "H—O" de verdade, que é o que interessa em base.
    """
    linhas = (len(smiles) + por_linha - 1) // por_linha
    figura, eixos = plt.subplots(linhas, por_linha, figsize=(3.5 * por_linha, 3.0 * linhas))
    eixos = np.atleast_1d(eixos).ravel()
    for eixo in eixos[len(smiles):]:
        eixo.set_axis_off()
    for eixo, smi, titulo in zip(eixos, smiles, titulos):
        eixo.imshow(desenho_png(smi, ionico=ionico, mostrar_h=mostrar_h))
        eixo.set_axis_off()
        eixo.set_title(titulo, fontsize=10, pad=6)
    figura.tight_layout()
    return figura

print("funções de desenho carregadas")

# %% md
## 1. Conceito principal

**Base** é a substância que, em água, libera íons $\mathrm{OH^-}$ (hidroxila). Essa é a definição de **Arrhenius**, e vale para os hidróxidos metálicos. Na visão de **Brønsted–Lowry**, base é quem **recebe** um próton — o que explica a amônia, que não tem $\mathrm{OH^-}$ na fórmula mas produz $\mathrm{OH^-}$ ao reagir com a água.

$$
\boxed{\mathrm{BOH_{(aq)} \rightarrow B^+_{(aq)} + OH^-_{(aq)}}}
\qquad
\mathrm{NH_3 + H_2O \rightleftharpoons NH_4^+ + OH^-}
$$

Para base fraca, escreve-se a constante de basicidade:

$$
K_b = \frac{[\mathrm{NH_4^+}]\,[\mathrm{OH^-}]}{[\mathrm{NH_3}]}
\qquad
[\mathrm{OH^-}] \approx \sqrt{K_b\,c}
$$

E a ponte com o pH, a $25\ ^\circ\mathrm{C}$:

$$
\boxed{\mathrm{pH} + \mathrm{pOH} = 14}
\qquad
\mathrm{pOH} = -\log[\mathrm{OH^-}]
$$

Bases **fortes** (NaOH, KOH, LiOH, Ca(OH)$_2$, Ba(OH)$_2$) dissociam praticamente 100%: a conta do pOH é direta. Bases **fracas** (NH$_3$, aminas, Mg(OH)$_2$, Al(OH)$_3$) precisam do $K_b$.

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** uma solução de $\mathrm{NaOH}$ a $0{,}010\ \mathrm{mol/L}$ tem qual pOH e qual pH? E uma solução de amônia de mesma concentração ($K_b = 1{,}8\times10^{-5}$)? Qual a massa de $\mathrm{NaOH}$ em $400\ \mathrm{mL}$ dessa solução?

**Passo 1 — base forte, dissociação total:**

$$
[\mathrm{OH^-}] = 0{,}010\ \mathrm{mol/L}
\qquad
\mathrm{pOH} = -\log_{10}(0{,}010) = 2{,}00
$$

$$
\mathrm{pH} = 14 - 2{,}00 = \boxed{12{,}00}
$$

**Passo 2 — base fraca (amônia):**

$$
[\mathrm{OH^-}] \approx \sqrt{1{,}8\times10^{-5}\cdot 0{,}010} = 4{,}24\times10^{-4}\ \mathrm{mol/L}
$$

$$
\mathrm{pOH} = 3{,}37 \qquad \mathrm{pH} = 14 - 3{,}37 = \boxed{10{,}63}
$$

**Passo 3 — massa** ($M_{\mathrm{NaOH}} = 40{,}00\ \mathrm{g/mol}$), com $n = cV = 0{,}010 \cdot 0{,}400 = 4{,}0\times10^{-3}\ \mathrm{mol}$:

$$
m = n\,M = 4{,}0\times10^{-3} \cdot 40{,}00 = \boxed{0{,}16\ \mathrm{g}}
$$

# %% md
## 3. Conferindo o cálculo com Python

# %% code title="Contas do exemplo"
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

C = 0.010        # mol/L
V = 0.400        # L
KB_NH3 = 1.8e-5

massa_molar = lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s))
formula = lambda s: rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))

M_naoh = massa_molar("[Na+].[OH-]")
n = C * V
m = n * M_naoh

pOH_forte = -math.log10(C)
pH_forte = 14 - pOH_forte

OH_fraco = math.sqrt(KB_NH3 * C)
pOH_fraco = -math.log10(OH_fraco)
pH_fraco = 14 - pOH_fraco

print(f"Fórmulas:            {formula('[Na+].[OH-]')} e {formula('N')}")
print(f"Massa molar NaOH:    {M_naoh:.2f} g/mol")
print(f"Quantidade:          {n:.4f} mol")
print(f"Massa de NaOH:       {m:.3f} g")
print(f"pOH do NaOH:         {pOH_forte:.2f}")
print(f"pH do NaOH:          {pH_forte:.2f}")
print(f"[OH-] da amônia:     {OH_fraco:.2e} mol/L")
print(f"pH da amônia:        {pH_fraco:.2f}")

# %% md
## 4. Visualização

### 4.1 A amônia em 3D — e o ângulo que quase fecha o tetraedro

O $\mathrm{NH_3}$ tem **3** ligações e **1 par de elétrons não ligante** no nitrogênio (azul). A repulsão desse par empurra as ligações N–H e o ângulo fica menor que o de um tetraedro perfeito ($109{,}5^\circ$). Gire a figura e olhe o rótulo.

> A tabela de livro dá $107{,}3^\circ$. A estrutura aqui é otimizada por campo de força (MMFF), que costuma medir perto de $106^\circ$–$108^\circ$: o importante é o **valor abaixo de $109{,}5^\circ$**, causado pelo par de elétrons livre.

# %% code title="Amônia: geometria medida"
visor, angulo = geometria_3d("N", atomo=0, cor="#3050f8")
print(f"ângulo H–N–H medido: {angulo:.1f}°   (livro: 107,3° | tetraedro ideal: 109,5°)")
visor.show()

# %% md
### 4.2 O íon hidróxido: a razão de toda base de Arrhenius

O $\mathrm{OH^-}$ aparece incluso nos hidróxidos metálicos. Aqui, o ânion sozinho em 3D (a carga negativa está no oxigênio, vermelho).

# %% code title="Íon hidróxido em 3D"
molecula_3d("[OH-]", rotulo="OH-", largura=460, altura=320)

# %% md
### 4.3 Hidróxidos metálicos em 3D

Nos hidróxidos, o cátion metálico fica **separado** dos ânions $\mathrm{OH^-}$ — não há ligação covalente entre eles, é atração iônica. Por isso o desenho mostra as espécies lado a lado: é exatamente essa separação que acontece quando o composto se dissolve.

# %% code title="NaOH em 3D interativo"
molecula_3d("[Na+].[OH-]", rotulo="NaOH", altura=330)

# %% code title="Ca(OH)2 em 3D interativo"
molecula_3d("[Ca+2].[OH-].[OH-]", rotulo="Ca(OH)2", altura=330)

# %% code title="Mg(OH)2 (leite de magnésia) em 3D interativo"
molecula_3d("[Mg+2].[OH-].[OH-]", rotulo="Mg(OH)2", altura=330)

# %% md
### 4.4 Os hidróxidos em 2D, com os risquinhos

# %% code title="Galeria 2D das bases"
galeria_2d(
    ["[Na+].[OH-]", "[K+].[OH-]", "[Ca+2].[OH-].[OH-]",
     "[Mg+2].[OH-].[OH-]", "[Al+3].[OH-].[OH-].[OH-]", "N"],
    ["NaOH\n(soda cáustica)", "KOH\n(potassa cáustica)", "Ca(OH)$_2$\n(cal hidratada)",
     "Mg(OH)$_2$\n(leite de magnésia)", "Al(OH)$_3$\n(antiácido)", "NH$_3$\n(amônia)"],
    ionico=True,      # pontilhado = atração eletrostática (não é ligação covalente)
    mostrar_h=True,   # aqui o H do "H—O" é o assunto: tem de aparecer
)

# %% md
### 4.5 A amônia em 2D — quando a base não tem $\mathrm{OH^-}$

No desenho 2D, o $\mathrm{NH_3}$ é só o nitrogênio com três hidrogênios; nenhuma hidroxila. A base aparece na reação com a água, que "rouba" um próton do $\mathrm{NH_3}$ e deixa um $\mathrm{OH^-}$ livre.

# %% code title="Amônia em 2D com hidrogênios explícitos"
estrutura_2d("N", legenda="amonia (NH3)", largura=460, altura=340, mostrar_h=True)

# %% code title="A reação com a água, do lado dos produtos"
galeria_2d(["[NH4+]", "[OH-]"], ["NH$_4^+$ (íon amônio)", "OH$^-$ (hidroxila)"], por_linha=2)

# %% md
### 4.6 Força das bases pelo $\mathrm{p}K_b$

Mesma lógica do $\mathrm{p}K_a$: **quanto menor o $\mathrm{p}K_b$, mais forte a base**. Os hidróxidos de metais alcalinos são bases fortes; a amônia e as aminas são fracas. (Valores aproximados, medidos em água; para os hidróxidos solúveis a dissociação é total e o $\mathrm{p}K_b$ só serve de comparação.)

# %% code title="pKb das bases comuns"
nomes = ["NaOH", "KOH", "Ca(OH)$_2$", "Ba(OH)$_2$", "NH$_3$", "CH$_3$NH$_2$", "Mg(OH)$_2$"]
pKbs = [-0.6, -0.6, -1.4, -0.7, 4.75, 3.36, 2.60]
cores = ["#c44e52" if p < 0 else "#4c72b0" for p in pKbs]

figura, eixo = plt.subplots(figsize=(8, 4.4))
eixo.barh(nomes[::-1], pKbs[::-1], color=cores[::-1])
for y, p in enumerate(pKbs[::-1]):
    # rótulo DENTRO da barra, na ponta: não encosta no eixo dos nomes nem no 0
    eixo.annotate(f"{p:+.2f}".replace(".", ","), (p, y),
                  xytext=(6 if p < 0 else -6, 0), textcoords="offset points",
                  va="center", ha="left" if p < 0 else "right",
                  fontsize=9, color="white")
eixo.axvline(0, color="#333", linewidth=1)
eixo.set_xlabel("p$K_b$  (menor = mais forte)")
eixo.set_title("Força das bases: vermelho = forte (dissocia 100%), azul = fraca")
eixo.grid(axis="x", alpha=0.3)
figura.tight_layout()

# %% md
### 4.7 pH conforme a concentração: base forte × base fraca

# %% code title="pH x concentração (bases)"
concentracoes = np.logspace(-4, 0, 80)
pOH_forte = -np.log10(concentracoes)
pOH_fraco = -np.log10(np.sqrt(KB_NH3 * concentracoes))

plt.figure(figsize=(7.5, 4.5))
plt.plot(concentracoes, 14 - pOH_forte, color="#c44e52", label="NaOH (forte)")
plt.plot(concentracoes, 14 - pOH_fraco, color="#4c72b0", label="NH$_3$ (fraca)")
plt.axhline(7, color="#888", linestyle=":", label="pH 7 = neutro")
plt.xscale("log")
plt.xlabel("Concentração (mol/L)")
plt.ylabel("pH")
plt.title("Bases: quanto mais diluída, mais perto do neutro")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.show()

# %% md
## 5. Como interpretar o resultado

- **Base forte é conta direta:** $[\mathrm{OH^-}]$ é a concentração da base **vezes o número de hidroxilas** — em $0{,}010\ \mathrm{mol/L}$ de $\mathrm{Ca(OH)_2}$ há $0{,}020\ \mathrm{mol/L}$ de $\mathrm{OH^-}$ (pOH $1{,}70$).
- **pOH primeiro, pH depois.** Trocar a ordem é o erro mais comum: quem dá o pH direto é $[\mathrm{H^+}]$, não $[\mathrm{OH^-}]$.
- **Base fraca tem pH mais baixo** que base forte de mesma concentração — no gráfico 4.7 a curva azul fica sempre abaixo.
- **A amônia é o caso que derruba a definição de Arrhenius:** não tem $\mathrm{OH^-}$ na fórmula e ainda assim é base. Brønsted–Lowry resolve.
- **No 3D**, o $\mathrm{NH_3}$ mede **abaixo de** $109{,}5^\circ$ (a tabela diz $107{,}3^\circ$) por causa do par de elétrons livre — o campo de força pode dar $106^\circ$–$108^\circ$, mas o efeito é o mesmo.
- Volume de solução sempre em **litros** em $n = cV$.

# %% md
## 6. Exercício de fixação

Um frasco tem $300\ \mathrm{mL}$ de solução de $\mathrm{Ca(OH)_2}$ a $0{,}020\ \mathrm{mol/L}$.

1. Qual a concentração de $\mathrm{OH^-}$?
2. Qual o pOH e o pH?
3. Qual a massa de $\mathrm{Ca(OH)_2}$ no frasco?
4. Se o frasco fosse de $\mathrm{NaOH}$ a $0{,}020\ \mathrm{mol/L}$, o pH seria maior ou menor?

Tente prever antes de executar a célula seguinte.

# %% code title="Exercício"
C_EX, V_EX = 0.020, 0.300
M_EX = massa_molar("[Ca+2].[OH-].[OH-]")

OH_EX = 2 * C_EX                      # duas hidroxilas por fórmula
pOH_EX = -math.log10(OH_EX)
pH_EX = 14 - pOH_EX
n_EX = C_EX * V_EX
m_EX = n_EX * M_EX
pH_naoh = 14 + math.log10(C_EX)

print(f"1) [OH-]:                 {OH_EX:.3f} mol/L")
print(f"2) pOH = {pOH_EX:.2f}  ->  pH = {pH_EX:.2f}")
print(f"3) massa:                 {m_EX:.3f} g   (n = {n_EX:.4f} mol)")
print(f"4) pH se fosse NaOH:      {pH_naoh:.2f} -> "
      f"{'maior' if pH_naoh > pH_EX else 'menor'} "
      "(NaOH tem só 1 hidroxila por fórmula)")

# %% md
## 7. Resumo

- Base, em água, libera $\mathrm{OH^-}$; doa par de elétrons / recebe próton (Brønsted–Lowry).
- $\mathrm{pOH} = -\log[\mathrm{OH^-}]$ e $\mathrm{pH} = 14 - \mathrm{pOH}$ (a $25\ ^\circ\mathrm{C}$).
- **Base forte:** $\mathrm{OH^-}$ = concentração × nº de hidroxilas. **Base fraca:** $[\mathrm{OH^-}] \approx \sqrt{K_b c}$.
- Fortes comuns: NaOH, KOH, LiOH, $\mathrm{Ca(OH)_2}$, $\mathrm{Ba(OH)_2}$.
- Fracas comuns: $\mathrm{NH_3}$, aminas, $\mathrm{Mg(OH)_2}$, $\mathrm{Al(OH)_3}$ (usados como antiácidos).
- Base + ácido → sal + água (neutralização) — veja o notebook de **sais** na pasta irmã.

## Referências

- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- FELTRE, R. *Química Geral*. Moderna.
