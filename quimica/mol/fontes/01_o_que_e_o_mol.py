# -*- coding: utf-8 -*-
"""Fonte do notebook 01_o_que_e_o_mol.ipynb (série mol).

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/mol/fontes/01_o_que_e_o_mol.py quimica/mol/01_o_que_e_o_mol.ipynb
"""

# %% md
# O que é o mol? Número de Avogadro e massa molar
### Mol, massa molar e a ponte entre o mundo microscópico e o macroscópico

**Objetivo:** Entender o mol como unidade de quantidade de matéria, calcular a massa molar a partir da fórmula molecular e converter massa em número de partículas (e vice-versa).

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. As cores seguem a convenção dos químicos (C = cinza, H = branco, O = vermelho, N = azul).

# %% code title="Preparar o ambiente (só precisa rodar uma vez)"
import importlib.util, subprocess, sys

def _tem(modulo):
    return importlib.util.find_spec(modulo) is not None

for pacote in ("rdkit", "py3Dmol"):
    if not _tem(pacote):
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", pacote], check=True)

import rdkit, py3Dmol
print("rdkit", rdkit.__version__, "| py3Dmol", py3Dmol.__version__)

# %% md
## 1. Conceito principal

Químicos não contam moléculas uma a uma: elas são numerosas demais. O **mol** é a unidade que resolve isso — é uma "dúzia química":

$$\boxed{1\ \mathrm{mol} = N_A\ \text{entidades} = 6{,}022\times10^{23}\ \mathrm{mol^{-1}}}$$

$N_A$ é o **número de Avogadro**. Desde 2019 ele tem valor exato, $N_A = 6{,}02214076\times10^{23}\ \mathrm{mol^{-1}}$; nas contas de escola, $6{,}02\times10^{23}$ já basta.

A **massa molar** $M$ (em $\mathrm{g/mol}$) é a massa de 1 mol da substância. Ela sai da soma das massas atômicas da fórmula. Para a água, $\mathrm{H_2O}$:

$$\boxed{M(\mathrm{H_2O}) = 2\cdot 1{,}008 + 16{,}00 = 18{,}02\ \mathrm{g/mol}}$$

Com $M$ e $N_A$ na mão, duas conversões resolvem quase tudo:

$$n = \frac{m}{M} \qquad\text{e}\qquad N = n\cdot N_A$$

onde $n$ é a quantidade de matéria (mol), $m$ a massa (g) e $N$ o número de partículas.

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** quantos mols e quantas moléculas há em $36{,}0\ \mathrm{g}$ de água?

**Passo 1 — massa molar.** Já calculada acima: $M = 18{,}02\ \mathrm{g/mol}$.

**Passo 2 — de gramas a mols.**

$$n = \frac{m}{M} = \frac{36{,}0}{18{,}02} = \boxed{2{,}00\ \mathrm{mol}}$$

**Passo 3 — de mols a moléculas.**

$$N = n\cdot N_A = 2{,}00 \cdot 6{,}022\times10^{23} = \boxed{1{,}20\times10^{24}\ \text{moléculas}}$$

Repare na ordem de grandeza: $36\ \mathrm{g}$ de água (um copo pequeno) já contém mais de **um trilhão de trilhões** de moléculas. É por isso que o mol existe.

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula lê a massa molar direto da fórmula (via RDKit, que usa as massas atômicas oficiais) e repete as duas conversões.

# %% code title="Contas do exemplo"
from rdkit import Chem
from rdkit.Chem import Descriptors

SMILES = "O"              # água: H2O
NA = 6.022e23             # número de Avogadro (mol^-1)
m_agua = 36.0             # g

mol = Chem.MolFromSmiles(SMILES)
M = Descriptors.MolWt(mol)          # g/mol
n = m_agua / M                      # mol
N = n * NA                          # moléculas

print(f"Fórmula molecular: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Massa molar:       {M:.2f} g/mol")
print(f"Quantidade:        {n:.2f} mol")
print(f"Número de moléculas: {N:.3e}")

# %% md
## 4. Visualização

### 4.1 A molécula de água em 3D — arraste para girar

A figura abaixo é interativa. Arraste para girar em 360°, use a roda do mouse para aproximar. Repare que a água **não é uma molécula plana**: o oxigênio (vermelho) fica no vértice de um ângulo de cerca de $104{,}5^\circ$, com os dois hidrogênios (brancos) apontando para baixo.

# %% code title="Água em 3D interativo"
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def molecula_3d(smiles, largura=640, altura=420):
    """SMILES -> visualização 3D interativa (arrastar = girar, rolar = zoom)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))   # H explícitos: aparecem na figura
    parametros = AllChem.ETKDGv3()
    parametros.randomSeed = 0xC0FFEE               # mesmo desenho em toda execução
    if AllChem.EmbedMolecule(mol, parametros) != 0:
        raise RuntimeError("não foi possível gerar a geometria 3D")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)          # relaxa as distâncias de ligação
    else:
        AllChem.UFFOptimizeMolecule(mol)

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.28, "colorscheme": "Jmol"}})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    return visor

molecula_3d(SMILES)

# %% md
### 4.2 A mesma molécula em 2D, com os risquinhos das ligações

No desenho 2D cada traço é uma ligação: **um** traço = simples, **dois** = dupla. A água tem duas ligações simples (O–H), e os hidrogênios aparecem aqui porque são o assunto da aula.

# %% code title="Estrutura 2D"
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

rdDepictor.SetPreferCoordGen(True)

def estrutura_2d(smiles, legenda="", largura=520, altura=360, mostrar_h=False):
    mol = Chem.MolFromSmiles(smiles)
    if mostrar_h:
        mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.drawOptions().addStereoAnnotation = True
    desenho.DrawMolecule(mol, legend=legenda)
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

estrutura_2d(SMILES, legenda="Agua (H2O)", mostrar_h=True)   # legenda do RDKit: sem acento

# %% md
### 4.3 Um mol de cada: compare as moléculas

> **Detalhe técnico que vale a pena saber:** o desenho do RDKit corrompe acentos e subscritos na legenda (`água` sai `gua`). Por isso as legendas dentro das figuras do RDKit ficam em ASCII; o texto explicativo com acento vai no markdown.

Cada molécula abaixo representa **1 mol** da substância — ou seja, $6{,}02\times10^{23}$ cópias dela. O número de moléculas é o mesmo; o que muda é a massa (seção 4.4).

# %% code title="Galeria 2D"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

nomes = ["agua (H2O)", "metano (CH4)", "amonia (NH3)",
         "dioxido de carbono (CO2)", "oxigenio (O2)", "etanol (C2H6O)"]
smiles = ["O", "C", "N", "O=C=O", "O=O", "CCO"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in smiles],
                              molsPerRow=3, subImgSize=(320, 260), legends=nomes, useSVG=True)
# dentro do notebook o RDKit já devolve um objeto SVG; fora dele, devolve texto
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 Massa de 1 mol de cada substância

# %% code title="Massa molar das substâncias"
import matplotlib.pyplot as plt

substancias = ["metano\n(CH₄)", "água\n(H₂O)", "amônia\n(NH₃)",
               "oxigênio\n(O₂)", "gás carbônico\n(CO₂)", "etanol\n(C₂H₆O)", "glicose\n(C₆H₁₂O₆)"]
massas = [16.04, 18.02, 17.03, 32.00, 44.01, 46.07, 180.16]   # g/mol (RDKit)

plt.figure(figsize=(8, 4.5))
barras = plt.bar(substancias, massas, color="#4c72b0")
for barra, m in zip(barras, massas):
    plt.annotate(f"{m:.2f}", (barra.get_x() + barra.get_width() / 2, m),
                 xytext=(0, 4), textcoords="offset points", ha="center")
plt.ylabel("Massa de 1 mol (g/mol)")
plt.title("Um mol de cada: mesmo número de moléculas, massas diferentes")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# %% md
## 5. Como interpretar o resultado

- A massa molar calculada pelo RDKit ($18{,}02\ \mathrm{g/mol}$) coincide com a do cálculo à mão: a fórmula química é a mesma conta.
- O mol é uma **contagem**, não uma massa: 1 mol de qualquer substância tem $6{,}02\times10^{23}$ partículas, mas a massa varia (de $16\ \mathrm{g}$ no metano a $180\ \mathrm{g}$ na glicose).
- No 3D, a água é angular, não linear: gire a figura e veja que os dois hidrogênios não estão em linha reta com o oxigênio.
- O erro mais comum em prova é confundir massa atômica (de um átomo) com massa molar (de 1 mol da substância): o valor numérico é o mesmo, mas a unidade muda de $\mathrm{u}$ para $\mathrm{g/mol}$.

# %% md
## 6. Exercício de fixação

Quantos mols e quantas moléculas há em $88{,}0\ \mathrm{g}$ de gás carbônico ($\mathrm{CO_2}$)?

Tente prever o resultado antes de executar a célula seguinte.

# %% code title="Exercício"
from rdkit import Chem
from rdkit.Chem import Descriptors

M_CO2 = Descriptors.MolWt(Chem.MolFromSmiles("O=C=O"))   # g/mol
n_ex = 88.0 / M_CO2
N_ex = n_ex * NA
print(f"M(CO2) = {M_CO2:.2f} g/mol")
print(f"n = {n_ex:.2f} mol")
print(f"N = {N_ex:.3e} moléculas")

# %% md
## 7. Resumo

- 1 mol = $6{,}02\times10^{23}$ entidades (número de Avogadro, $N_A$).
- Massa molar $M$ = soma das massas atômicas da fórmula, em $\mathrm{g/mol}$.
- $n = m/M$ (gramas → mols) e $N = n\cdot N_A$ (mols → partículas).
- Visualização 3D interativa mostra a geometria; 2D com risquinhos mostra as ligações; o gráfico mostra que 1 mol de substâncias diferentes tem massas diferentes.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book) — verbetes "mole" e "Avogadro constant".
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- BROWN, T. et al. *Química: a Ciência Central*. Pearson.
