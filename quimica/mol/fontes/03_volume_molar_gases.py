# -*- coding: utf-8 -*-
"""Fonte do notebook 03_volume_molar_gases.ipynb (série mol).

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/mol/fontes/03_volume_molar_gases.py \
        quimica/mol/03_volume_molar_gases.ipynb
"""

# %% md
# Volume molar: o mol dos gases
### Lei de Avogadro, CNTP e a equação PV = nRT

**Objetivo:** Calcular o volume de um gás a partir do número de mols (e vice-versa) usando o volume molar e a equação do gás ideal.

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

A **lei de Avogadro** diz que volumes iguais de gases ideais, na mesma temperatura e pressão, contêm o mesmo número de moléculas. Consequência direta: 1 mol de **qualquer** gás ideal ocupa o mesmo volume, chamado de **volume molar** $V_m$.

Nas **CNTP** (condições normais de temperatura e pressão: $0^\circ\mathrm{C} = 273{,}15\ \mathrm{K}$ e $1\ \mathrm{atm}$):

$$\boxed{V_m = 22{,}4\ \mathrm{L/mol}} \qquad\Longrightarrow\qquad \boxed{V = n\cdot V_m}$$

A $25^\circ\mathrm{C}$ e $1\ \mathrm{atm}$ (condições ambiente), o volume molar é maior: $V_m \approx 24{,}5\ \mathrm{L/mol}$.

O volume molar vem da equação do gás ideal, $PV = nRT$, com $R = 0{,}08206\ \mathrm{L\cdot atm\cdot mol^{-1}\cdot K^{-1}}$:

$$V = \frac{nRT}{P}$$

Repare que $V_m = RT/P$ **não depende do gás**: $\mathrm{CO_2}$, $\mathrm{O_2}$ e $\mathrm{H_2}$ ocupam o mesmo volume por mol (enquanto forem ideais).

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** qual o volume de $0{,}500\ \mathrm{mol}$ de $\mathrm{CO_2}$ nas CNTP?

**Pelo volume molar:**

$$V = n\cdot V_m = 0{,}500 \cdot 22{,}4 = \boxed{11{,}2\ \mathrm{L}}$$

**Pela equação do gás ideal (conferência):**

$$V = \frac{nRT}{P} = \frac{0{,}500 \cdot 0{,}08206 \cdot 273{,}15}{1} = \boxed{11{,}2\ \mathrm{L}}$$

Os dois caminhos concordam — como deve ser, porque o volume molar é a equação do gás ideal "pronta" para as CNTP.

**E no sentido contrário?** $44{,}8\ \mathrm{L}$ de $\mathrm{O_2}$ nas CNTP correspondem a $n = 44{,}8/22{,}4 = 2{,}00\ \mathrm{mol}$, ou seja, $2{,}00 \cdot 32{,}00 = 64{,}0\ \mathrm{g}$ de oxigênio.

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula calcula o volume pelos dois caminhos (volume molar e $PV = nRT$) e faz a conversão inversa.

# %% code title="Contas do exemplo"
R = 0.08206        # L·atm/(mol·K)
T_CNTP = 273.15    # K (0 °C)
P = 1.0            # atm
VM_CNTP = 22.4     # L/mol nas CNTP

n = 0.500          # mol de CO2
V_por_vm = n * VM_CNTP
V_por_gas_ideal = n * R * T_CNTP / P
print(f"Volume por Vm:        {V_por_vm:.2f} L")
print(f"Volume por PV = nRT:  {V_por_gas_ideal:.2f} L")

# Inverso: 44,8 L de O2 nas CNTP -> mol -> gramas
n_o2 = 44.8 / VM_CNTP
m_o2 = n_o2 * 32.00
print(f"\n44,8 L de O2 = {n_o2:.2f} mol = {m_o2:.1f} g")

# %% md
## 4. Visualização

### 4.1 A molécula de CO₂ em 3D — arraste para girar

A figura abaixo é interativa. Arraste para girar em 360°, use a roda do mouse para aproximar. O $\mathrm{CO_2}$ é **linear**: o carbono (cinza) fica no meio e os dois oxigênios (vermelhos) em lados opostos, com ângulo de $180^\circ$. É essa molécula que ocupa $11{,}2\ \mathrm{L}$ no exemplo acima.

# %% code title="CO2 em 3D interativo"
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

molecula_3d("O=C=O")   # gás carbônico: CO2

# %% md
### 4.2 A mesma molécula em 2D, com os risquinhos das ligações

No desenho 2D cada traço é uma ligação: **dois** traços = ligação dupla. O $\mathrm{CO_2}$ tem duas ligações duplas (C=O), o que o torna uma molécula linear.

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

estrutura_2d("O=C=O", legenda="Dioxido de carbono (CO2)")   # legenda do RDKit: sem acento

# %% md
### 4.3 Gases comuns: todos com o mesmo volume molar

> **Detalhe técnico que vale a pena saber:** o desenho do RDKit corrompe acentos e subscritos na legenda. Por isso as legendas dentro das figuras do RDKit ficam em ASCII; o texto explicativo com acento vai no markdown.

Cada molécula abaixo representa 1 mol do gás — e, nas CNTP, cada um desses mols ocupa $22{,}4\ \mathrm{L}$, independentemente da identidade do gás.

# %% code title="Galeria 2D"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

nomes = ["gas carbonico (CO2)", "oxigenio (O2)", "hidrogenio (H2)",
         "amonia (NH3)", "metano (CH4)", "nitrogenio (N2)"]
smiles = ["O=C=O", "O=O", "[H][H]", "N", "C", "[N]#[N]"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in smiles],
                              molsPerRow=3, subImgSize=(320, 260), legends=nomes, useSVG=True)
# dentro do notebook o RDKit já devolve um objeto SVG; fora dele, devolve texto
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 Volume em função do número de mols

# %% code title="Volume x mols"
import numpy as np
import matplotlib.pyplot as plt

mols = np.linspace(0, 3, 100)
V_cntp = mols * 22.4        # CNTP: 273,15 K, 1 atm
V_25c = mols * 24.5         # 25 °C, 1 atm (aprox.)

plt.figure(figsize=(7.5, 4.5))
plt.plot(mols, V_cntp, label="CNTP (0 °C, 1 atm): Vm = 22,4 L/mol")
plt.plot(mols, V_25c, "--", label="25 °C, 1 atm: Vm ≈ 24,5 L/mol")
plt.xlabel("Quantidade de matéria (mol)")
plt.ylabel("Volume (L)")
plt.title("Volume molar: a reta que passa pela origem")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# %% md
## 5. Como interpretar o resultado

- O volume molar é a **inclinação** da reta $V \times n$: cada mol adiciona $22{,}4\ \mathrm{L}$ nas CNTP, qualquer que seja o gás.
- Aquecer o gás aumenta o volume molar ($22{,}4 \to 24{,}5\ \mathrm{L/mol}$ de $0^\circ\mathrm{C}$ a $25^\circ\mathrm{C}$): é a equação $V = nRT/P$ em ação.
- No 3D, o $\mathrm{CO_2}$ é linear ($180^\circ$): gire a figura e veja os dois oxigênios em lados opostos do carbono.
- O erro mais comum é usar $22{,}4\ \mathrm{L/mol}$ fora das CNTP. Se a temperatura ou a pressão forem diferentes, use $V = nRT/P$.

# %% md
## 6. Exercício de fixação

**a)** Qual o volume de $1{,}50\ \mathrm{mol}$ de $\mathrm{H_2}$ nas CNTP? E a massa correspondente?

**b)** Quantos mols há em $24{,}5\ \mathrm{L}$ de $\mathrm{N_2}$ a $25^\circ\mathrm{C}$ e $1\ \mathrm{atm}$?

Tente prever os resultados antes de executar a célula seguinte.

# %% code title="Exercício"
from rdkit import Chem
from rdkit.Chem import Descriptors

# a) 1,50 mol de H2 nas CNTP
V_a = 1.50 * VM_CNTP
M_H2 = Descriptors.MolWt(Chem.MolFromSmiles("[H][H]"))
m_a = 1.50 * M_H2
print(f"a) V = {V_a:.1f} L")
print(f"   M(H2) = {M_H2:.2f} g/mol -> m = {m_a:.2f} g")

# b) 24,5 L de N2 a 25 °C, 1 atm
n_b = 24.5 / 24.5
print(f"b) n = {n_b:.2f} mol")

# %% md
## 7. Resumo

- Volume molar: $V_m = 22{,}4\ \mathrm{L/mol}$ nas CNTP; $V = n\cdot V_m$.
- Equação do gás ideal: $V = nRT/P$, com $R = 0{,}08206\ \mathrm{L\cdot atm\cdot mol^{-1}\cdot K^{-1}}$.
- O volume molar não depende da identidade do gás (ideal); depende de $T$ e $P$.
- Visualização 3D interativa mostra a geometria; 2D com risquinhos mostra as ligações duplas; o gráfico mostra $V \times n$ como uma reta.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book) — verbetes "molar volume" e "Avogadro law".
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- BROWN, T. et al. *Química: a Ciência Central*. Pearson.
