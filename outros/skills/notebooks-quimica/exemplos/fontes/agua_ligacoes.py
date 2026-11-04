# -*- coding: utf-8 -*-
"""Fonte do notebook exemplo_agua_ligacoes.ipynb.

Gere o .ipynb com:
    python scripts/criar_notebook.py fontes/agua_ligacoes.py ../exemplos/exemplo_agua_ligacoes.ipynb
"""

# %% md
# Água: geometria angular e ligações de hidrogênio
### Ligações químicas e forças intermoleculares

**Objetivo:** Entender por que a água é angular (104,5°) e como as ligações de hidrogênio explicam o ponto de ebulição alto dela em comparação com moléculas parecidas.

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** arraste com o mouse para girar em 360°. Na primeira figura, as setas mostram as direções das ligações O–H e o ângulo medido aparece sobre a molécula.

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

A água tem fórmula $\mathrm{H_2O}$: um oxigênio com dois pares de elétrons não ligantes. Esses pares "empurram" as ligações O–H, e o resultado é uma molécula **angular**, com ângulo próximo ao tetraédrico:

$$
\boxed{\theta_{\mathrm{H-O-H}} \approx 104{,}5^\circ}
$$

A diferença entre os átomos gera um dipolo elétrico: o lado do oxigênio fica parcialmente negativo ($\delta^-$) e o lado dos hidrogênios, positivo ($\delta^+$). Por isso uma molécula de água atrai o oxigênio da vizinha — a **ligação de hidrogênio**:

$$
\mathrm{O-H \cdots O} \qquad \text{(energia típica: } 20\ \mathrm{kJ/mol}\text{)}
$$

Compare com a ligação covalente O–H, de cerca de $460\ \mathrm{kJ/mol}$: a ligação de hidrogênio é umas 20 vezes mais fraca, mas como cada molécula pode fazer até quatro, o efeito coletivo é enorme.

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** por que a água ferve a $100\ ^\circ\mathrm{C}$ e o $\mathrm{H_2S}$ (mesma família, enxofre logo abaixo do oxigênio) ferve a $-60\ ^\circ\mathrm{C}$?

**Passo 1 — massas molares.** $\mathrm{H_2O}$: $18\ \mathrm{g/mol}$; $\mathrm{H_2S}$: $34\ \mathrm{g/mol}$. O $\mathrm{H_2S}$ é mais pesado, então *deveria* ferver mais alto.

**Passo 2 — forças intermoleculares.** O enxofre é maior e menos eletronegativo: o $\mathrm{H_2S}$ é quase apolar em comparação e não forma ligações de hidrogênio. A água forma.

**Passo 3 — energia necessária para ferver.**

$$
\Delta H_{\text{vap}}(\mathrm{H_2O}) \approx 40{,}7\ \mathrm{kJ/mol} \qquad
\Delta H_{\text{vap}}(\mathrm{H_2S}) \approx 18{,}7\ \mathrm{kJ/mol}
$$

A razão é $\dfrac{40{,}7}{18{,}7} \approx \boxed{2{,}2}$: a água precisa de mais que o dobro de energia por mol para passar ao vapor, apesar de ser mais leve.

# %% md
## 3. Conferindo o cálculo com Python

# %% code title="Comparação numérica"
dados = {
    "H2O": {"massa": 18.02, "vaporizacao": 40.7, "ebulicao": 100.0},
    "H2S": {"massa": 34.08, "vaporizacao": 18.7, "ebulicao": -60.0},
}
for nome, d in dados.items():
    print(f"{nome}: massa {d['massa']:.2f} g/mol | "
          f"ΔHvap {d['vaporizacao']:.1f} kJ/mol | Teb {d['ebulicao']:.0f} °C")
razao = dados["H2O"]["vaporizacao"] / dados["H2S"]["vaporizacao"]
print(f"razão das energias de vaporização: {razao:.2f}")

# %% md
## 4. Visualização

### 4.1 Geometria da água em 3D — gire e confira o ângulo

Arraste a figura. As setas azuis saem do oxigênio em direção aos dois hidrogênios; o número mostra o ângulo medido na geometria otimizada (deve dar perto de $104{,}5^\circ$).

# %% code title="Água em 3D com ângulo medido"
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import py3Dmol

def agua_com_angulo(largura=640, altura=460):
    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    parametros = AllChem.ETKDGv3()
    parametros.randomSeed = 42
    AllChem.EmbedMolecule(mol, parametros)
    AllChem.MMFFOptimizeMolecule(mol)

    conformacao = mol.GetConformer()
    posicoes = np.array([list(conformacao.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    o, h1, h2 = 0, 1, 2

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.28, "colorscheme": "Jmol"}})
    inicio = {"x": float(posicoes[o][0]), "y": float(posicoes[o][1]), "z": float(posicoes[o][2])}
    for h in (h1, h2):
        visor.addArrow({"start": inicio,
                        "end": {"x": float(posicoes[h][0]), "y": float(posicoes[h][1]),
                                "z": float(posicoes[h][2])},
                        "radius": 0.045, "radiusRatio": 1.8, "mid": 0.85,
                        "color": "#1f6feb"})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    return visor, mol

visor, mol_agua = agua_com_angulo()
angulo = Chem.rdMolTransforms.GetAngleDeg(mol_agua.GetConformer(), 1, 0, 2)
print(f"ângulo H–O–H = {angulo:.1f}°")
visor.show()

# %% md
### 4.2 As duas moléculas em 2D e a ligação de hidrogênio com risquinhos

O pontilhado entre o H de uma água e o O da outra representa a ligação de hidrogênio. No código, esse pontilhado é uma ligação de **ordem zero** (`Chem.BondType.ZERO`), que o RDKit desenha tracejada.

# %% code title="Dímero de água com ligação de hidrogênio"
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

def agua_dimer():
    mol = Chem.AddHs(Chem.MolFromSmiles("O.O"))     # dois oxigênios + hidrogênios explícitos
    oxigenios = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]

    # hidrogênio que pertence ao segundo oxigênio
    h_do_segundo = [v.GetIdx() for v in mol.GetAtomWithIdx(oxigenios[1]).GetNeighbors()
                    if v.GetSymbol() == "H"][0]

    rw = Chem.RWMol(mol)
    rw.AddBond(h_do_segundo, oxigenios[0], Chem.BondType.ZERO)   # pontilhado = ligação de H
    dimero = rw.GetMol()
    dimero.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(dimero, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)

    rdDepictor.Compute2DCoords(dimero)
    desenho = rdMolDraw2D.MolDraw2DSVG(600, 340)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.DrawMolecule(dimero, legend="duas aguas: ligacao de hidrogenio pontilhada")
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

agua_dimer()

# %% md
### 4.3 Ligações covalentes simples: um, dois e três risquinhos

# %% code title="Ordens de ligação"
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import SVG

exemplos = ["CC", "C=C", "C#N", "c1ccccc1"]
# legenda dentro do desenho do RDKit precisa ser ASCII (ele corrompe acentos)
legendas = ["etano (C-C, 1 traco)", "eteno (C=C, 2 tracos)",
            "HCN (C#N, 3 tracos)", "benzeno (aromatico)"]
imagem = Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in exemplos], molsPerRow=4,
                              subImgSize=(330, 300), legends=legendas, useSVG=True)
imagem if hasattr(imagem, "data") else SVG(imagem)

# %% md
### 4.4 Título corrigido com matplotlib (acento e subscrito de verdade)

A figura abaixo é a mesma estrutura, mas o título vem do matplotlib: aí `ácido`,
`ç` e o subscrito `₃` aparecem corretamente. Use este padrão para as figuras que
vão para o material final.

# %% code title="Estrutura 2D com título correto"
import io
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

def figura_com_legenda(smiles, titulo, largura=520, altura=360):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DCairo(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.drawOptions().padding = 0.12
    desenho.DrawMolecule(mol)                 # legenda fora do RDKit
    desenho.FinishDrawing()
    imagem = np.asarray(Image.open(io.BytesIO(desenho.GetDrawingText())))

    figura, eixo = plt.subplots(figsize=(largura / 100, (altura + 46) / 100))
    eixo.imshow(imagem)
    eixo.set_axis_off()
    eixo.set_title(titulo, fontsize=11, pad=8)
    figura.tight_layout()
    return figura

figura_com_legenda("O", "Água (H$_2$O) — geometria angular")   # $...$ = subscrito

# %% md
### 4.5 Ponto de ebulição da família 16

# %% code title="Ponto de ebulição x massa molar"
import matplotlib.pyplot as plt

compostos = ["H₂O", "H₂S", "H₂Se", "H₂Te"]
massa = [18.0, 34.1, 81.0, 129.6]
ebulicao = [100.0, -60.0, -41.0, -2.0]

plt.figure(figsize=(7.5, 4.5))
plt.plot(massa, ebulicao, "o-", color="#4c72b0", linewidth=2, label="dados reais")
plt.annotate("a água \"foge\" da tendência\n(ligações de hidrogênio)",
             xy=(massa[0], ebulicao[0]), xytext=(45, 55),
             arrowprops={"arrowstyle": "->", "color": "#c44e52"}, color="#c44e52")
plt.xlabel("Massa molar (g/mol)")
plt.ylabel("Ponto de ebulição (°C)")
plt.title("Tendência da família 16: a água é a exceção")
plt.grid(alpha=0.3)
plt.legend(loc="lower right")
plt.show()

# %% md
## 5. Como interpretar o resultado

- O ângulo medido no 3D (perto de $104{,}5^\circ$) **não** é o tetraédrico puro ($109{,}5^\circ$): os dois pares de elétrons livres comprimem a molécula.
- No desenho 2D do dímero, o pontilhado é a ligação de hidrogênio; ele fica entre o H (positivo) de uma água e o O (negativo) da outra. Note que a ligação covalente O–H continua com um traço cheio.
- No gráfico, a tendência normal é o ponto de ebulição subir com a massa; a água está muito acima da reta. Essa é a "assinatura" da ligação de hidrogênio.
- No 3D, gire a molécula e observe que a água é **plana** (os três átomos no mesmo plano), diferente do metano, tridimensional.

# %% md
## 6. Exercício de fixação

O $\mathrm{NH_3}$ (amônia) também faz ligações de hidrogênio. Estime o ângulo H–N–H com o mesmo procedimento e compare com o da água.

Tente prever antes de executar.

# %% code title="Exercício: ângulo da amônia"
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def angulo_hxh(smiles, central, largura=640, altura=420):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parametros = AllChem.ETKDGv3()
    parametros.randomSeed = 7
    AllChem.EmbedMolecule(mol, parametros)
    AllChem.MMFFOptimizeMolecule(mol)

    c = mol.GetConformer()
    vizinhos = [n.GetIdx() for n in mol.GetAtomWithIdx(central).GetNeighbors()]
    angulos = [Chem.rdMolTransforms.GetAngleDeg(c, vizinhos[i], central, vizinhos[i + 1])
               for i in range(len(vizinhos) - 1)]

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.28, "colorscheme": "Jmol"}})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    return visor, angulos

visor_nh3, angulos_nh3 = angulo_hxh("N", 0)
print("ângulos H–N–H:", ", ".join(f"{a:.1f}°" for a in angulos_nh3))
visor_nh3.show()

# %% md
## 7. Resumo

- $\mathrm{H_2O}$ é angular ($\approx 104{,}5^\circ$), polar e plana.
- Cada água pode doar dois e aceitar dois hidrogênios $\rightarrow$ até quatro ligações de hidrogênio.
- Ligação de hidrogênio ($\sim 20\ \mathrm{kJ/mol}$) é fraca sozinha, forte no conjunto: explica o ponto de ebulição anômalo.
- No notebook: 3D interativo para geometria, 2D com pontilhado para as forças intermoleculares, gráfico para a tendência.

## Referências

- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
