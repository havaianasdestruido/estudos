# -*- coding: utf-8 -*-
"""Fonte do notebook 02_conversoes_massa_mol_particulas.ipynb (série mol).

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/mol/fontes/02_conversoes_massa_mol_particulas.py \
        quimica/mol/02_conversoes_massa_mol_particulas.ipynb
"""

# %% md
# Conversões: massa, mol e número de partículas
### O triângulo do mol e como ir de gramas a moléculas

**Objetivo:** Converter entre massa (g), quantidade de matéria (mol) e número de partículas nos dois sentidos, usando a massa molar $M$ e o número de Avogadro $N_A$.

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

Toda conversão do mol passa por um "triângulo" de três grandezas: **massa** $m$ (g), **quantidade de matéria** $n$ (mol) e **número de partículas** $N$. As pontes são a massa molar $M$ e o número de Avogadro $N_A$:

$$\boxed{n = \frac{m}{M}} \qquad \boxed{N = n\cdot N_A} \qquad \boxed{m = \frac{N\cdot M}{N_A}}$$

- Para ir de **gramas a mols**, divida pela massa molar.
- Para ir de **mols a partículas**, multiplique por $N_A = 6{,}02\times10^{23}\ \mathrm{mol^{-1}}$.
- Para voltar (**partículas a gramas**), divida por $N_A$ e multiplique por $M$.

Um cuidado de linguagem: "partículas" pode ser átomos, moléculas ou fórmulas unitárias. No sal de cozinha ($\mathrm{NaCl}$), 1 mol tem $6{,}02\times10^{23}$ **fórmulas** de NaCl — não moléculas, porque o sal é uma rede iônica, não um conjunto de moléculas soltas.

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** quantos mols e quantas moléculas há em $8{,}02\ \mathrm{g}$ de metano ($\mathrm{CH_4}$)?

**Passo 1 — massa molar do metano.**

$$M(\mathrm{CH_4}) = 12{,}01 + 4\cdot 1{,}008 = \boxed{16{,}04\ \mathrm{g/mol}}$$

**Passo 2 — gramas a mols.**

$$n = \frac{m}{M} = \frac{8{,}02}{16{,}04} = \boxed{0{,}500\ \mathrm{mol}}$$

**Passo 3 — mols a moléculas.**

$$N = n\cdot N_A = 0{,}500 \cdot 6{,}022\times10^{23} = \boxed{3{,}01\times10^{23}\ \text{moléculas}}$$

**E no sentido contrário?** Se tivéssemos $1{,}00\times10^{23}$ moléculas de $\mathrm{CO_2}$:

$$n = \frac{N}{N_A} = \frac{1{,}00\times10^{23}}{6{,}022\times10^{23}} = 0{,}166\ \mathrm{mol}, \qquad m = n\cdot M = 0{,}166\cdot 44{,}01 = 7{,}31\ \mathrm{g}$$

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula faz as conversões nos dois sentidos: gramas → mols → moléculas, e moléculas → mols → gramas.

# %% code title="Contas do exemplo"
from rdkit import Chem
from rdkit.Chem import Descriptors

NA = 6.022e23                       # número de Avogadro (mol^-1)

# Ida: 8,02 g de metano -> mol -> moléculas
M_CH4 = Descriptors.MolWt(Chem.MolFromSmiles("C"))
n_ch4 = 8.02 / M_CH4
N_ch4 = n_ch4 * NA
print(f"M(CH4) = {M_CH4:.2f} g/mol")
print(f"n = {n_ch4:.3f} mol")
print(f"N = {N_ch4:.3e} moléculas")

# Volta: 1,00e23 moléculas de CO2 -> mol -> gramas
M_CO2 = Descriptors.MolWt(Chem.MolFromSmiles("O=C=O"))
n_co2 = 1.00e23 / NA
m_co2 = n_co2 * M_CO2
print(f"\nM(CO2) = {M_CO2:.2f} g/mol")
print(f"n = {n_co2:.3f} mol")
print(f"m = {m_co2:.2f} g")

# %% md
## 4. Visualização

### 4.1 A molécula de metano em 3D — arraste para girar

A figura abaixo é interativa. Arraste para girar em 360°, use a roda do mouse para aproximar. O metano é **tetraédrico**: o carbono (cinza) fica no centro e os quatro hidrogênios (brancos) apontam para os vértices de um tetraedro, com ângulos de $109{,}5^\circ$.

# %% code title="Metano em 3D interativo"
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

molecula_3d("C")   # metano: CH4

# %% md
### 4.2 A mesma molécula em 2D, com os risquinhos das ligações

No desenho 2D cada traço é uma ligação simples. O metano tem **quatro** ligações C–H — por isso a fórmula é $\mathrm{CH_4}$.

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

estrutura_2d("C", legenda="Metano (CH4)", mostrar_h=True)   # legenda do RDKit: sem acento

# %% md
### 4.3 Outras substâncias para praticar as conversões

> **Detalhe técnico que vale a pena saber:** o desenho do RDKit corrompe acentos e subscritos na legenda. Por isso as legendas dentro das figuras do RDKit ficam em ASCII; o texto explicativo com acento vai no markdown.

# %% code title="Galeria 2D"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

nomes = ["metano (CH4)", "agua (H2O)", "amonia (NH3)",
         "dioxido de carbono (CO2)", "oxigenio (O2)", "etanol (C2H6O)"]
smiles = ["C", "O", "N", "O=C=O", "O=O", "CCO"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in smiles],
                              molsPerRow=3, subImgSize=(320, 260), legends=nomes, useSVG=True)
# dentro do notebook o RDKit já devolve um objeto SVG; fora dele, devolve texto
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 Massa de $1{,}00\times10^{23}$ partículas de cada substância

# %% code title="Massa de 1,00e23 partículas"
import matplotlib.pyplot as plt

substancias = ["metano\n(CH₄)", "água\n(H₂O)", "amônia\n(NH₃)",
               "oxigênio\n(O₂)", "gás carbônico\n(CO₂)", "etanol\n(C₂H₆O)"]
massas = [16.04, 18.02, 17.03, 32.00, 44.01, 46.07]   # g/mol
N_ALVO = 1.00e23
# m = N * M / NA : mesmo N, massas diferentes
m_alvo = [N_ALVO * M / NA for M in massas]

plt.figure(figsize=(8, 4.5))
barras = plt.bar(substancias, m_alvo, color="#4c72b0")
for barra, m in zip(barras, m_alvo):
    plt.annotate(f"{m:.2f} g", (barra.get_x() + barra.get_width() / 2, m),
                 xytext=(0, 4), textcoords="offset points", ha="center")
plt.ylabel("Massa (g)")
plt.title("Mesmo número de partículas (1,00×10²³), massas diferentes")
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# %% md
## 5. Como interpretar o resultado

- Converter é só multiplicar ou dividir por $M$ e por $N_A$ — a ordem importa: gramas → mols (divide por $M$), mols → partículas (multiplica por $N_A$).
- O gráfico mostra o ponto central do mol: fixar o **número** de partículas não fixa a **massa**. Quem fixa a massa é a massa molar.
- No 3D, o metano é tetraédrico ($109{,}5^\circ$): gire a figura e conte os quatro hidrogênios ao redor do carbono.
- Confira sempre a ordem de grandeza: $N_A \approx 6\times10^{23}$ significa que $0{,}1\ \mathrm{mol}$ já tem $\approx 6\times10^{22}$ partículas — nunca $6\times10^{24}$.

# %% md
## 6. Exercício de fixação

**a)** Quantos mols e quantas moléculas há em $4{,}25\ \mathrm{g}$ de amônia ($\mathrm{NH_3}$)?

**b)** Quantos mols e quantos gramas há em $2{,}00\times10^{23}$ moléculas de água?

Tente prever os resultados antes de executar a célula seguinte.

# %% code title="Exercício"
from rdkit import Chem
from rdkit.Chem import Descriptors

M_NH3 = Descriptors.MolWt(Chem.MolFromSmiles("N"))
n_a = 4.25 / M_NH3
N_a = n_a * NA
print(f"a) M(NH3) = {M_NH3:.2f} g/mol")
print(f"   n = {n_a:.3f} mol")
print(f"   N = {N_a:.3e} moléculas")

M_H2O = Descriptors.MolWt(Chem.MolFromSmiles("O"))
n_b = 2.00e23 / NA
m_b = n_b * M_H2O
print(f"b) M(H2O) = {M_H2O:.2f} g/mol")
print(f"   n = {n_b:.3f} mol")
print(f"   m = {m_b:.2f} g")

# %% md
## 7. Resumo

- $n = m/M$ (gramas ↔ mols) e $N = n\cdot N_A$ (mols ↔ partículas); $m = N\cdot M/N_A$ faz o caminho de volta.
- "Partículas" pode ser átomos, moléculas ou fórmulas unitárias (ex.: NaCl).
- Mesmo número de partículas não significa mesma massa: quem manda na massa é a massa molar.
- Visualização 3D interativa mostra a geometria; 2D com risquinhos mostra as ligações; o gráfico compara massas para um mesmo $N$.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book) — verbetes "mole" e "Avogadro constant".
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- BROWN, T. et al. *Química: a Ciência Central*. Pearson.
