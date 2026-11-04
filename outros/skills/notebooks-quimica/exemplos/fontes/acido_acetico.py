# -*- coding: utf-8 -*-
"""Rascunho-fonte do notebook exemplo_acido_acetico.ipynb.

Gere o .ipynb com:
    python scripts/criar_notebook.py rascunho_acido_acetico.py ../exemplos/exemplo_acido_acetico.ipynb
"""

# %% md
# Ácido Acético: massa molar, soluções e pH
### Funções orgânicas, mol e concentração

**Objetivo:** Calcular a massa molar do ácido acético a partir da fórmula molecular, preparar uma solução em mol/L e estimar o pH, olhando a molécula em 3D interativo e em 2D com as ligações.

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. As cores seguem a convenção dos químicos (C = cinza, H = branco, O = vermelho).

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

O ácido acético é um ácido carboxílico de fórmula molecular $\mathrm{C_2H_4O_2}$ e fórmula estrutural $\mathrm{CH_3COOH}$. A massa molar vem da soma das massas atômicas:

$$
\boxed{M = 2\cdot 12{,}01 + 4\cdot 1{,}008 + 2\cdot 16{,}00 = 60{,}05\ \mathrm{g/mol}}
$$

Em solução aquosa ele se ioniza parcialmente:

$$
\mathrm{CH_3COOH \rightleftharpoons H^+ + CH_3COO^-}, \qquad K_a = 1{,}8\times10^{-5}
$$

Por ser um ácido fraco, o pH não é o de um ácido forte da mesma concentração. Para uma concentração $c$, a aproximação clássica é

$$
[\mathrm{H^+}] \approx \sqrt{K_a\,c}, \qquad \mathrm{pH} = -\log_{10}[\mathrm{H^+}]
$$

# %% md
## 2. Exemplo resolvido passo a passo

Preparo de $250\ \mathrm{mL}$ de solução de ácido acético a $0{,}100\ \mathrm{mol/L}$ (é o vinagre "de laboratório", bem mais concentrado que o vinagre de cozinha, que tem cerca de $0{,}8\ \mathrm{mol/L}$).

**Passo 1 — massa necessária.**

$$
n = c \cdot V = 0{,}100 \cdot 0{,}250 = 0{,}0250\ \mathrm{mol}
$$
$$
m = n \cdot M = 0{,}0250 \cdot 60{,}05 = \boxed{1{,}50\ \mathrm{g}}
$$

**Passo 2 — pH.**

$$
[\mathrm{H^+}] \approx \sqrt{1{,}8\times10^{-5} \cdot 0{,}100} = 1{,}34\times10^{-3}\ \mathrm{mol/L}
$$
$$
\mathrm{pH} = -\log_{10}(1{,}34\times10^{-3}) = \boxed{2{,}87}
$$

Esse valor é típico de ácido fraco: se fosse forte, o pH seria $1{,}00$.

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula implementa as duas fórmulas com as grandezas e unidades do exemplo.

# %% code title="Contas do exemplo"
SMILES = "CC(=O)O"          # ácido acético: CH3-C(=O)OH
C_USADA = 0.100             # mol/L
V_USADO = 0.250             # L
KA = 1.8e-5

# As massas atômicas oficiais vêm do RDKit (iguais à tabela periódica)
from rdkit import Chem
from rdkit.Chem import Descriptors
mol = Chem.MolFromSmiles(SMILES)
M = Descriptors.MolWt(mol)                     # g/mol
n = C_USADA * V_USADO                          # mol
m = n * M                                      # g
H = (KA * C_USADA) ** 0.5                      # mol/L
pH = -__import__("math").log10(H)

print(f"Fórmula molecular: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Massa molar:       {M:.2f} g/mol")
print(f"Quantidade:        {n:.4f} mol")
print(f"Massa a pesar:     {m:.2f} g")
print(f"[H+] estimado:     {H:.2e} mol/L")
print(f"pH estimado:       {pH:.2f}")

# %% md
## 4. Visualização

### 4.1 A molécula em 3D — arraste para girar

A figura abaixo é interativa. Arraste para girar em 360°, use a roda do mouse para aproximar. Repare no grupo carboxila ($\mathrm{-COOH}$): o oxigênio da hidroxila (branco ligado ao vermelho) é o que perde o $\mathrm{H^+}$.

# %% code title="Ácido acético em 3D interativo"
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
    AllChem.MMFFOptimizeMolecule(mol)              # relaxa as distâncias de ligação

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

No desenho 2D cada traço é uma ligação: **um** traço = simples, **dois** = dupla. O hidrogênio do $-COOH$ aparece porque ele é justamente o que sai na ionização; os hidrogênios do carbono já estão implícitos, como manda a convenção do esqueleto.

# %% code title="Estrutura 2D"
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

rdDepictor.SetPreferCoordGen(True)

def estrutura_2d(smiles, legenda="", largura=520, altura=360):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.drawOptions().addStereoAnnotation = True
    desenho.DrawMolecule(mol, legend=legenda)
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

estrutura_2d(SMILES, legenda="Acido acetico (CH3COOH)")   # legenda do RDKit: sem acento

# %% md
### 4.3 Comparando com outros ácidos carboxílicos

> **Detalhe técnico que vale a pena saber:** o desenho do RDKit corrompe acentos e subscritos na legenda (`ácido` sai `cido`). Por isso as legendas dentro das figuras do RDKit ficam em ASCII; o texto explicativo com acento vai no markdown, e quando precisa de título bonito usamos matplotlib (veja a figura 4.5).

# %% code title="Galeria 2D"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

nomes = ["acido metanoico (formico)", "acido etanoico (acetico)", "acido propanoico (propionico)"]
smiles = ["C(=O)O", "CC(=O)O", "CCC(=O)O"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in smiles],
                              molsPerRow=3, subImgSize=(320, 260), legends=nomes, useSVG=True)
# dentro do notebook o RDKit já devolve um objeto SVG; fora dele, devolve texto
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 O pH conforme a concentração

# %% code title="pH x concentração"
import numpy as np
import matplotlib.pyplot as plt

concentracoes = np.logspace(-3, 0, 60)          # de 0,001 a 1 mol/L
ph_fraco = -np.log10((KA * concentracoes) ** 0.5)
ph_forte = -np.log10(concentracoes)

plt.figure(figsize=(7.5, 4.5))
plt.plot(concentracoes, ph_fraco, label="ácido acético (fraco, Ka = 1,8×10⁻⁵)")
plt.plot(concentracoes, ph_forte, "--", label="ácido forte de mesma concentração")
plt.xscale("log")
plt.xlabel("Concentração (mol/L)")
plt.ylabel("pH")
plt.title("Por que o ácido fraco tem pH mais alto")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# %% md
## 5. Como interpretar o resultado

- A massa molar calculada pelo RDKit ($60{,}05\ \mathrm{g/mol}$) coincide com a do cálculo à mão: a fórmula química é a mesma conta.
- No 3D, a molécula **não é plana**: o $\mathrm{CH_3}$ é tetraédrico e o carbono da carboxila é trigonal plano. Gire a figura para ver os dois planos diferentes.
- A curva de pH mostra que a diferença entre ácido forte e fraco é de cerca de 1,5 unidade de pH nas concentrações usuais — e que a diferença diminui quando a solução fica muito diluída.
- O item que mais erra em prova é esquecer de converter o volume para litros: $250\ \mathrm{mL} = 0{,}250\ \mathrm{L}$, nunca $250$.

# %% md
## 6. Exercício de fixação

Qual é a massa de ácido acético necessária para preparar $500\ \mathrm{mL}$ de solução a $0{,}050\ \mathrm{mol/L}$? E qual o pH esperado?

Tente prever o resultado antes de executar a célula seguinte.

# %% code title="Exercício"
C_EX, V_EX = 0.050, 0.500
n_ex = C_EX * V_EX
m_ex = n_ex * M
pH_ex = -__import__("math").log10((KA * C_EX) ** 0.5)
print(f"n = {n_ex:.4f} mol")
print(f"massa = {m_ex:.2f} g")
print(f"pH esperado = {pH_ex:.2f}")

# %% md
## 7. Resumo

- Ácido acético: $\mathrm{CH_3COOH}$, massa molar $60{,}05\ \mathrm{g/mol}$.
- $n = cV$ e $m = nM$ resolvem o preparo de soluções.
- Ácido fraco: $[\mathrm{H^+}] \approx \sqrt{K_a c}$ e $\mathrm{pH} = -\log[\mathrm{H^+}]$.
- Visualização 3D interativa mostra geometria; 2D com risquinhos mostra a estrutura e as ligações duplas.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
