# -*- coding: utf-8 -*-
"""Rascunho-fonte de quimica/sais/01_sais_neutralizacao_solubilidade_hidrolise.ipynb.

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/fontes/sais.py quimica/sais/01_sais_neutralizacao_solubilidade_hidrolise.ipynb
"""

# %% md
# Sais: neutralização, solubilidade e hidrólise
### Funções inorgânicas: os sais

**Objetivo:** Reconhecer a função sal (NaCl, CaCO$_3$, NaHCO$_3$, NH$_4$Cl, KNO$_3$, CuSO$_4$), escrever a reação de neutralização que produz cada um, calcular concentração comum ($\%\,m/v$), quantidade de matéria e massa, e entender por que a solução de alguns sais não fica com pH 7.

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. Nas figuras de sais, o **cátion** e o **ânion** aparecem como íons separados: é assim que eles estão na solução.

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

def massa_molar(smiles):
    """Massa molar (g/mol) a partir do SMILES do composto (ou de seus íons)."""
    return Descriptors.MolWt(Chem.MolFromSmiles(smiles))

def formula(smiles):
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))

print("funções de desenho carregadas")

# %% md
## 1. Conceito principal

**Sal** é o composto iônico formado, em geral, pela reação de **neutralização** entre um ácido e uma base:

$$
\boxed{\mathrm{ácido + base \rightarrow sal + água}}
$$

$$
\mathrm{HCl + NaOH \rightarrow NaCl + H_2O}
\qquad
\mathrm{H_2CO_3 + Ca(OH)_2 \rightarrow CaCO_3 + 2\,H_2O}
$$

Todo sal é formado por um **cátion** (vindo da base) e um **ânion** (vindo do ácido). O nome segue a regra *ânion de ácido + cátion da base*: cloreto de sódio, carbonato de cálcio, nitrato de potássio, sulfato de cobre(II).

Sais podem ser **neutros** (NaCl), **ácidos** (NH$_4$Cl, um sal que deixa a solução ácida), **básicos** (CaCO$_3$, NaHCO$_3$) ou **hidratados** (CuSO$_4\cdot5\mathrm{H_2O}$, o "azul da pedra"). E o pH da solução **não é sempre 7**: depende de quem veio do ácido/base forte ou fraco — detalhe na seção 4.6.

Concentração comum, muito usada na prática (soro, soro fisiológico, água sanitária):

$$
\boxed{c\ (\%\,m/v) = \frac{m\ \mathrm{(g)}}{V\ \mathrm{(mL)}}\times 100}
\qquad
c\ \mathrm{(mol/L)} = \frac{m}{M\cdot V_{L}}
$$

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** o **soro fisiológico** é uma solução de $\mathrm{NaCl}$ a $0{,}9\ \%\,m/v$. Qual sua concentração em mol/L? Quantos gramas de $\mathrm{NaCl}$ há num frasco de $500\ \mathrm{mL}$? E qual é a massa de $\mathrm{CaCO_3}$ em $0{,}25\ \mathrm{mol}$ de calcário?

**Passo 1 — de $\%\,m/v$ para mol/L.** $0{,}9\ \%\,m/v$ significa $0{,}9\ \mathrm{g}$ em $100\ \mathrm{mL}$, ou seja $9{,}0\ \mathrm{g/L}$:

$$
c = \frac{m}{M\,V} = \frac{9{,}0}{58{,}44\cdot 1{,}00} = \boxed{0{,}154\ \mathrm{mol/L}}
$$

**Passo 2 — massa no frasco de $500\ \mathrm{mL}$** ($V = 0{,}500\ \mathrm{L}$):

$$
n = c\,V = 0{,}154 \cdot 0{,}500 = 0{,}0770\ \mathrm{mol}
\qquad
m = n\,M = 0{,}0770 \cdot 58{,}44 = \boxed{4{,}5\ \mathrm{g}}
$$

*(O atalho: $0{,}9\ \%\,m/v$ em $500\ \mathrm{mL}$ → $0{,}9 \cdot 5 = 4{,}5\ \mathrm{g}$.)*

**Passo 3 — massa de $\mathrm{CaCO_3}$** ($M = 100{,}09\ \mathrm{g/mol}$):

$$
m = n\,M = 0{,}25 \cdot 100{,}09 = \boxed{25{,}0\ \mathrm{g}}
$$

# %% md
## 3. Conferindo o cálculo com Python

# %% code title="Contas do exemplo"
from rdkit.Chem import Descriptors, rdMolDescriptors

SMILES_NACL = "[Na+].[Cl-]"
SMILES_CACO3 = "[Ca+2].O=C([O-])[O-]"

PORCENTO_MV = 0.9      # % m/v do soro fisiológico
V_ML = 500.0           # mL do frasco
N_CARBONATO = 0.25     # mol

M_nacl = Descriptors.MolWt(Chem.MolFromSmiles(SMILES_NACL))
M_caco3 = Descriptors.MolWt(Chem.MolFromSmiles(SMILES_CACO3))

gramas_por_litro = PORCENTO_MV * 10          # 0,9 % m/v = 9,0 g/L
c_nacl = gramas_por_litro / M_nacl           # mol/L
m_frasco = PORCENTO_MV / 100 * V_ML          # % m/v * mL = gramas
n_frasco = m_frasco / M_nacl
m_caco3 = N_CARBONATO * M_caco3

print(f"Massas molares:  NaCl = {M_nacl:.2f} g/mol | CaCO3 = {M_caco3:.2f} g/mol")
print(f"Concentração:    {gramas_por_litro:.1f} g/L = {c_nacl:.3f} mol/L")
print(f"Massa no frasco: {m_frasco:.1f} g  ({n_frasco:.4f} mol em {V_ML:.0f} mL)")
print(f"CaCO3:           {m_caco3:.1f} g em {N_CARBONATO} mol")

# %% md
## 4. Visualização

### 4.1 O carbonato de cálcio em 3D — o íon carbonato é plano

O $\mathrm{CaCO_3}$ (calcita, mármore, casca de ovo) é formado pelo cátion $\mathrm{Ca^{2+}}$ e pelo ânion $\mathrm{CO_3^{2-}}$. Gire a figura e repare: o **carbono no centro com três oxigênios num plano**, ângulos de $120^\circ$ — a assinatura do carbono $sp^2$.

# %% code title="Carbonato de cálcio em 3D com ângulo medido"
visor, angulo = geometria_3d("[Ca+2].O=C([O-])[O-]", atomo=2, cor="#c44e52")
print(f"ângulo O–C–O medido: {angulo:.1f}°   (esperado: 120°)")
visor.show()

# %% md
### 4.2 O bicarbonato de sódio e o sulfato de cobre

O $\mathrm{NaHCO_3}$ é o bicarbonato de sódio (fermento químico, antiácido). O $\mathrm{CuSO_4}$ anidro é branco, mas com 5 águas de hidratação fica azul: é o teste clássico de presença de água.

> **Atenção a uma pegadinha de prova:** o bicarbonato **não** é chamado de "sal ácido" por ter H na fórmula; a classificação vem da **hidrólise** — ele vem de ácido fraco + base forte, então sua solução é **básica** (pH $\approx 8{,}3$, veja 4.6). Quem deixa a solução ácida é o $\mathrm{NH_4Cl}$.

# %% code title="Bicarbonato de sódio em 3D"
molecula_3d("[Na+].OC([O-])=O", rotulo="NaHCO3", altura=330)

# %% code title="Sulfato de cobre(II) em 3D"
molecula_3d("[Cu+2].[O-]S(=O)(=O)[O-]", rotulo="CuSO4", altura=330)

# %% code title="Cloreto de amônio em 3D (sal de ácido forte + base fraca)"
molecula_3d("[NH4+].[Cl-]", rotulo="NH4Cl", altura=330)

# %% md
### 4.3 Os sais em 2D, com os risquinhos

# %% code title="Galeria 2D dos sais"
galeria_2d(
    ["[Na+].[Cl-]", "[Ca+2].O=C([O-])[O-]", "[Na+].OC([O-])=O",
     "[NH4+].[Cl-]", "[K+].[O-][N+](=O)[O-]", "[Cu+2].[O-]S(=O)(=O)[O-]"],
    ["NaCl\n(cloreto de sódio)", "CaCO$_3$\n(carbonato de cálcio)",
     "NaHCO$_3$\n(bicarbonato de sódio)", "NH$_4$Cl\n(cloreto de amônio)",
     "KNO$_3$\n(nitrato de potássio)", "CuSO$_4$\n(sulfato de cobre II)"],
    ionico=True,     # pontilhado = atração eletrostática entre os íons
)

# %% md
### 4.4 As peças soltas: cátions e ânions

Em solução, o sal não existe como "molécula": existe como **íons separados**, rodeados de água. Estes são os íons mais comuns.

# %% code title="Cátions e ânions em 2D"
galeria_2d(["[Na+]", "[Ca+2]", "[NH4+]", "[Cl-]", "O=C([O-])[O-]", "[O-][N+](=O)[O-]"],
           ["Na$^+$", "Ca$^{2+}$", "NH$_4^+$", "Cl$^-$", "CO$_3^{2-}$", "NO$_3^-$"])

# %% md
### 4.5 Solubilidade em água: por que o nitrato "some" e o carbonato "afunda"

Solubilidade em $\mathrm{g}$ de sal por $100\ \mathrm{g}$ de água, em função da temperatura. Repare que o $\mathrm{KNO_3}$ dispara com o aquecimento (é a base do salitre e da cristalização fracionada), enquanto o $\mathrm{NaCl}$ é quase constante — é isso que permite separar os dois.

# %% code title="Solubilidade x temperatura"
temperaturas = np.array([0, 20, 40, 60, 80, 100])
solubilidade = {
    "NaCl":  np.array([35.7, 36.0, 36.6, 37.3, 38.4, 39.8]),
    "KCl":   np.array([27.6, 34.0, 40.0, 45.5, 51.1, 56.7]),
    "KNO$_3$": np.array([13.3, 31.6, 63.9, 110.0, 169.0, 246.0]),
}

plt.figure(figsize=(7.5, 4.6))
for (nome, valores), cor in zip(solubilidade.items(), ["#4c72b0", "#55a868", "#c44e52"]):
    plt.plot(temperaturas, valores, "o-", color=cor, label=nome)
plt.xlabel("Temperatura (°C)")
plt.ylabel("Solubilidade (g por 100 g de água)")
plt.title("Solubilidade de três sais: o nitrato é o sensível à temperatura")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# %% md
### 4.6 Por que a solução de sal não tem sempre pH 7

O sal "herda" a força do ácido e da base que o geraram:

| Sal | Vem de | pH da solução |
|---|---|---|
| NaCl | ácido forte + base forte | $\approx 7$ (neutro) |
| NH$_4$Cl | ácido forte + base fraca | $< 7$ (ácido) |
| NaHCO$_3$ | ácido fraco + base forte | $> 7$ (básico) |
| CH$_3$COONa | ácido fraco + base forte | $> 7$ (básico) |

# %% code title="pH de soluções salinas 0,1 mol/L"
sais = ["NaCl", "NH$_4$Cl", "CH$_3$COONa", "NaHCO$_3$", "KNO$_3$"]
phs = [7.00, 5.13, 8.87, 8.31, 7.00]
cores = ["#888888" if abs(p - 7) < 0.05 else ("#c44e52" if p < 7 else "#4c72b0")
         for p in phs]

figura, eixo = plt.subplots(figsize=(7.8, 4.3))
barras = eixo.bar(sais, phs, color=cores)
eixo.axhline(7, color="#333", linestyle="--", linewidth=1)
eixo.annotate("neutro", (4.35, 7.08), fontsize=9, color="#333")
for barra, p in zip(barras, phs):
    eixo.annotate(f"{p:.2f}".replace(".", ","), (barra.get_x() + barra.get_width() / 2, p),
                  xytext=(0, 4), textcoords="offset points", ha="center", fontsize=9)
eixo.set_ylim(0, 11)
eixo.set_ylabel("pH de solução 0,1 mol/L")
eixo.set_title("Hidrólise do sal: vermelho = ácido, azul = básico, cinza = neutro")
eixo.grid(axis="y", alpha=0.3)
figura.tight_layout()

# %% md
## 5. Como interpretar o resultado

- **A reação de neutralização conserva os átomos:** o cátion vem da base e o ânion vem do ácido. Confira sempre os coeficientes — com ácido/base dipróticos aparece o fator 2.
- **$\%\,m/v$ é atalho:** $x\ \%\,m/v$ equivale a $x\cdot 10\ \mathrm{g/L}$. Para mol/L, divida pela massa molar. É assim que se lê rótulo de soro ($0{,}9\%$) e de água sanitária ($2{,}0$–$2{,}5\%$ de $\mathrm{NaClO}$).
- **Solubilidade não é "regra de decoreba" pura:** nitratos e sais de sódio/potássio/amônio são solúveis; carbonatos, fosfatos e hidróxidos de metais pesados são pouco solúveis. O gráfico 4.5 mostra o efeito da temperatura, que a regra não mostra.
- **Sal não quer dizer pH 7:** a hidrólise do íon que veio do fraco "puxa" o pH para o lado do fraco. $\mathrm{NH_4^+}$ doa próton (solução ácida); $\mathrm{HCO_3^-}$ e $\mathrm{CH_3COO^-}$ recebem próton (solução básica).
- No 3D, o íon carbonato é **plano**, com $120^\circ$ — se você mediu perto disso, a geometria está certa.

# %% md
## 6. Exercício de fixação

Uma receita de **soro caseiro** usa $3{,}5\ \mathrm{g}$ de sal de cozinha ($\mathrm{NaCl}$) e $20\ \mathrm{g}$ de açúcar em $1{,}0\ \mathrm{L}$ de água.

1. Qual a concentração de $\mathrm{NaCl}$ em $\%\,m/v$?
2. Qual a concentração de $\mathrm{NaCl}$ em mol/L?
3. Qual a quantidade de matéria de $\mathrm{NaCl}$ em um copo de $200\ \mathrm{mL}$ desse soro?
4. Esse soro é hipertônico, isotônico ou hipotônico em relação ao soro fisiológico ($0{,}9\ \%\,m/v$)?

Tente prever antes de executar a célula seguinte.

# %% code title="Exercício"
MASSA_SAL, VOLUME_L = 3.5, 1.0
COPO_ML = 200.0

porcento_ex = MASSA_SAL / (VOLUME_L * 1000) * 100      # % m/v
c_ex = MASSA_SAL / (M_nacl * VOLUME_L)                 # mol/L
n_copo = c_ex * (COPO_ML / 1000)
comparacao = ("hipertônico" if porcento_ex > PORCENTO_MV else
              "isotônico" if abs(porcento_ex - PORCENTO_MV) < 0.05 else "hipotônico")

print(f"1) concentração:      {porcento_ex:.2f} % m/v")
print(f"2) em mol/L:          {c_ex:.3f} mol/L")
print(f"3) em 200 mL:         {n_copo:.4f} mol  ({n_copo * M_nacl:.2f} g)")
print(f"4) em relação ao soro fisiológico (0,9 %): {comparacao}")

# %% md
## 7. Resumo

- Sal = composto iônico de cátion (da base) + ânion (do ácido); nasce da neutralização ácido + base.
- Nomenclatura: ânion + "de" + cátion (cloreto de sódio, nitrato de potássio, carbonato de cálcio).
- $\%\,m/v \leftrightarrow$ mol/L: $x\ \%\,m/v = 10x\ \mathrm{g/L}$; $c = m/(M V_L)$.
- Solubilidade depende do sal **e** da temperatura (nitratos muito sensíveis, cloretos quase não).
- Hidrólise: forte+forte → neutro; forte+fraco → pH do lado do forte (ácido ou básico).
- Sais aparecem em casa: sal de cozinha, bicarbonato, calcário, fertilizante (KNO$_3$, $\mathrm{NH_4Cl}$).

## Referências

- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- FELTRE, R. *Química Geral*. Moderna.
