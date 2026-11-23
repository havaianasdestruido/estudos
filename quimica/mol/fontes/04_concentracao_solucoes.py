# -*- coding: utf-8 -*-
"""Fonte do notebook 04_concentracao_solucoes.ipynb (série mol).

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/mol/fontes/04_concentracao_solucoes.py \
        quimica/mol/04_concentracao_solucoes.ipynb
"""

# %% md
# Concentração em mol/L e diluição
### Molaridade: quantos mols por litro de solução

**Objetivo:** Preparar soluções em mol/L (calcular a massa a pesar e o volume) e diluir soluções usando a conservação dos mols.

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

A **concentração em mol/L** (molaridade, símbolo $c$) diz quantos mols de soluto há em cada litro de **solução** (não de solvente!):

$$\boxed{c = \frac{n}{V}} \qquad\Longrightarrow\qquad n = c\cdot V \qquad\text{e}\qquad m = c\cdot V\cdot M$$

onde $V$ é o volume da solução em litros, $n$ os mols de soluto e $M$ a massa molar (g/mol).

Na **diluição**, adicionamos solvente à solução: o volume aumenta, mas os mols de soluto **não mudam**. Logo:

$$\boxed{c_1\cdot V_1 = c_2\cdot V_2}$$

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** quantos gramas de ureia ($\mathrm{CO(NH_2)_2}$) são necessários para preparar $250\ \mathrm{mL}$ de solução a $0{,}200\ \mathrm{mol/L}$?

**Passo 1 — massa molar da ureia.**

$$M = 12{,}01 + 16{,}00 + 2\cdot 14{,}01 + 4\cdot 1{,}008 = \boxed{60{,}06\ \mathrm{g/mol}}$$

**Passo 2 — mols necessários** (cuidado: $250\ \mathrm{mL} = 0{,}250\ \mathrm{L}$).

$$n = c\cdot V = 0{,}200 \cdot 0{,}250 = \boxed{0{,}0500\ \mathrm{mol}}$$

**Passo 3 — massa a pesar.**

$$m = n\cdot M = 0{,}0500 \cdot 60{,}06 = \boxed{3{,}00\ \mathrm{g}}$$

**Diluição:** para obter $200\ \mathrm{mL}$ de solução $0{,}100\ \mathrm{mol/L}$ a partir de uma solução $1{,}00\ \mathrm{mol/L}$, quanto da solução concentrada é preciso?

$$V_1 = \frac{c_2\cdot V_2}{c_1} = \frac{0{,}100\cdot 200}{1{,}00} = \boxed{20{,}0\ \mathrm{mL}}$$

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula calcula a massa a pesar e o volume de diluição, com a massa molar vinda do RDKit.

# %% code title="Contas do exemplo"
from rdkit import Chem
from rdkit.Chem import Descriptors

SMILES = "NC(=O)N"        # ureia: CO(NH2)2
c = 0.200                 # mol/L
V = 0.250                 # L (250 mL)

M = Descriptors.MolWt(Chem.MolFromSmiles(SMILES))   # g/mol
n = c * V
m = n * M
print(f"Fórmula molecular: {Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(SMILES))}")
print(f"Massa molar:       {M:.2f} g/mol")
print(f"Quantidade:        {n:.4f} mol")
print(f"Massa a pesar:     {m:.2f} g")

# Diluição: 1,00 mol/L -> 0,100 mol/L, volume final 200 mL
V1 = 0.100 * 200 / 1.00
print(f"\nDiluição: pipetar {V1:.1f} mL da solução 1,00 mol/L e completar até 200 mL")

# %% md
## 4. Visualização

### 4.1 A molécula de ureia em 3D — arraste para girar

A figura abaixo é interativa. Arraste para girar em 360°, use a roda do mouse para aproximar. A ureia tem um carbono central (cinza) ligado a um oxigênio (vermelho, ligação dupla) e a dois grupos $\mathrm{-NH_2}$ (nitrogênios azuis com hidrogênios brancos).

# %% code title="Ureia em 3D interativo"
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

No desenho 2D cada traço é uma ligação: **dois** traços = dupla (C=O), **um** traço = simples (C–N). Os hidrogênios dos grupos $\mathrm{-NH_2}$ aparecem porque fazem parte da fórmula.

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

estrutura_2d(SMILES, legenda="Ureia (CH4N2O)", mostrar_h=True)   # legenda do RDKit: sem acento

# %% md
### 4.3 Outros solutos comuns em laboratório

> **Detalhe técnico que vale a pena saber:** o desenho do RDKit corrompe acentos e subscritos na legenda. Por isso as legendas dentro das figuras do RDKit ficam em ASCII; o texto explicativo com acento vai no markdown.

# %% code title="Galeria 2D"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

nomes = ["ureia (CH4N2O)", "etanol (C2H6O)", "acido acetico (C2H4O2)", "glicose (C6H12O6)"]
smiles = ["NC(=O)N", "CCO", "CC(=O)O", "C(C(C(C(C(C=O)O)O)O)O)O"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in smiles],
                              molsPerRow=2, subImgSize=(320, 260), legends=nomes, useSVG=True)
# dentro do notebook o RDKit já devolve um objeto SVG; fora dele, devolve texto
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 Massa a pesar em função da concentração

# %% code title="Massa x concentração"
import numpy as np
import matplotlib.pyplot as plt

concentracoes = np.linspace(0, 1, 100)      # mol/L
V_L = 0.250                                  # 250 mL
m_ureia = concentracoes * V_L * M            # g

plt.figure(figsize=(7.5, 4.5))
plt.plot(concentracoes, m_ureia, label="ureia, 250 mL (M = 60,06 g/mol)")
plt.xlabel("Concentração (mol/L)")
plt.ylabel("Massa a pesar (g)")
plt.title("Massa a pesar cresce com a concentração")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# %% md
## 5. Como interpretar o resultado

- A massa molar calculada pelo RDKit ($60{,}06\ \mathrm{g/mol}$) coincide com a do cálculo à mão: a fórmula química é a mesma conta.
- O gráfico é uma reta que passa pela origem: dobrar a concentração dobra a massa a pesar (para volume fixo).
- No 3D, a ureia é **plana** ao redor do carbono central: gire a figura e veja que os átomos ficam (quase) no mesmo plano.
- Os dois erros clássicos: esquecer de converter mL para L ($250\ \mathrm{mL} = 0{,}250\ \mathrm{L}$) e usar o volume do solvente em vez do volume da solução.

# %% md
## 6. Exercício de fixação

**a)** Quantos gramas de $\mathrm{NaCl}$ são necessários para preparar $500\ \mathrm{mL}$ de solução a $0{,}100\ \mathrm{mol/L}$? ($M = 58{,}44\ \mathrm{g/mol}$.)

**b)** Quanto da solução $2{,}00\ \mathrm{mol/L}$ é preciso pipetar para obter $100\ \mathrm{mL}$ de solução $0{,}200\ \mathrm{mol/L}$?

Tente prever os resultados antes de executar a célula seguinte.

# %% code title="Exercício"
M_NACL = 58.44

# a) 500 mL de NaCl 0,100 mol/L
n_a = 0.100 * 0.500
m_a = n_a * M_NACL
print(f"a) n = {n_a:.4f} mol -> m = {m_a:.2f} g de NaCl")

# b) diluição: c1*V1 = c2*V2
V1_b = 0.200 * 100 / 2.00
print(f"b) pipetar {V1_b:.1f} mL da solução 2,00 mol/L e completar até 100 mL")

# %% md
## 7. Resumo

- $c = n/V$ (mol/L); $n = cV$; $m = cVM$.
- Diluição conserva os mols: $c_1V_1 = c_2V_2$.
- Sempre converta mL para L antes de calcular.
- Visualização 3D interativa mostra a geometria; 2D com risquinhos mostra as ligações; o gráfico mostra massa × concentração como uma reta.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book) — verbete "molarity".
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- BROWN, T. et al. *Química: a Ciência Central*. Pearson.
