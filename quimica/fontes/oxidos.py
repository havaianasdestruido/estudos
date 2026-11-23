# -*- coding: utf-8 -*-
"""Rascunho-fonte de quimica/oxidos/01_oxidos_acidicos_basicos_anfotericos.ipynb.

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/fontes/oxidos.py quimica/oxidos/01_oxidos_acidicos_basicos_anfotericos.ipynb
"""

# %% md
# Óxidos: ácidos, básicos, anfóteros e neutros
### Funções inorgânicas: os óxidos

**Objetivo:** Reconhecer a função óxido (CO$_2$, CO, SO$_2$, CaO, Na$_2$O, Al$_2$O$_3$, ZnO, NO), classificá-los em ácidos, básicos, anfóteros e neutros, prever o produto da reação com água e calcular quanto óxido se forma na combustão (estequiometria com massa molar).

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. As cores seguem a convenção dos químicos (C = cinza, O = vermelho, S = amarelo, N = azul, metais em cores variadas).

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
from rdkit.Chem import AllChem, Descriptors, rdDepictor, rdMolDescriptors
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
    """Massa molar (g/mol) a partir do SMILES do composto."""
    return Descriptors.MolWt(Chem.MolFromSmiles(smiles))

def formula(smiles):
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))

print("funções de desenho carregadas")

# %% md
## 1. Conceito principal

**Óxido** é todo composto binário em que o oxigênio é o elemento mais **eletronegativo** — ou seja, o oxigênio vem sempre em segundo lugar na fórmula: $\mathrm{CO_2}$, $\mathrm{CaO}$, $\mathrm{Al_2O_3}$, $\mathrm{SO_3}$.

A classificação depende do **comportamento em água**, e isso define tudo em prova:

| Classe | Comportamento | Exemplos | Reação com água |
|---|---|---|---|
| **Ácido** | óxido de ametal; forma ácido | $\mathrm{CO_2}$, $\mathrm{SO_2}$, $\mathrm{SO_3}$, $\mathrm{NO_2}$ | $\mathrm{CO_2 + H_2O \rightarrow H_2CO_3}$ |
| **Básico** | óxido de metal; forma base | $\mathrm{Na_2O}$, $\mathrm{CaO}$, $\mathrm{MgO}$, $\mathrm{Fe_2O_3}$ | $\mathrm{CaO + H_2O \rightarrow Ca(OH)_2}$ |
| **Anfótero** | reage com ácido **e** com base | $\mathrm{Al_2O_3}$, $\mathrm{ZnO}$, $\mathrm{PbO}$ | não reage de forma simples com a água |
| **Neutro** | não reage com água nem neutraliza | $\mathrm{CO}$, $\mathrm{NO}$, $\mathrm{N_2O}$ | — |

A ponte com as outras funções:

$$
\boxed{\text{óxido básico} + \text{água} \rightarrow \text{base}}
\qquad
\boxed{\text{óxido ácido} + \text{água} \rightarrow \text{ácido}}
$$

$$
\mathrm{óxido\ básico + óxido\ ácido \rightarrow sal}
\qquad
\mathrm{CaO + CO_2 \rightarrow CaCO_3}
$$

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** queimando $1{,}0\ \mathrm{kg}$ de carvão (considere carbono puro), quantos gramas de $\mathrm{CO_2}$ se formam? E se o carvão fosse de enxofre, quanto $\mathrm{SO_2}$ sairia?

**Passo 1 — a reação é 1 para 1 em mol:**

$$
\mathrm{C + O_2 \rightarrow CO_2}
\qquad
n_{\mathrm{CO_2}} = n_{\mathrm{C}}
$$

**Passo 2 — quantidade de matéria de carbono em $1{,}0\ \mathrm{kg}$:**

$$
n = \frac{m}{M} = \frac{1000}{12{,}01} = 83{,}3\ \mathrm{mol}
$$

**Passo 3 — massa de $\mathrm{CO_2}$** ($M = 44{,}01\ \mathrm{g/mol}$):

$$
m = n\,M = 83{,}3 \cdot 44{,}01 = \boxed{3{,}66\times10^{3}\ \mathrm{g} \approx 3{,}7\ \mathrm{kg}}
$$

**Passo 4 — se fosse enxofre** ($M_{\mathrm{S}} = 32{,}06$; $M_{\mathrm{SO_2}} = 64{,}06$): a razão de massas é $64{,}06/32{,}06 \approx 2{,}00$, então

$$
m_{\mathrm{SO_2}} = 1000 \cdot \frac{64{,}06}{32{,}06} = \boxed{2{,}00\times10^{3}\ \mathrm{g} = 2{,}0\ \mathrm{kg}}
$$

Repare: a massa **aumenta** porque o oxigênio do ar entra no produto — nada se perde, só se combina.

# %% md
## 3. Conferindo o cálculo com Python

# %% code title="Contas do exemplo"
MASSA_CARVAO_G = 1000.0     # g

M_C = massa_molar("[C]")           # átomo isolado: "C" sozinho vira CH4 no RDKit
M_CO2 = massa_molar("O=C=O")
M_S = massa_molar("[S]")
M_SO2 = massa_molar("O=S=O")

n_C = MASSA_CARVAO_G / M_C
m_CO2 = n_C * M_CO2
m_SO2 = MASSA_CARVAO_G * (M_SO2 / M_S)

print(f"Massas molares:   C = {M_C:.2f} | CO2 = {M_CO2:.2f} | S = {M_S:.2f} | SO2 = {M_SO2:.2f} g/mol")
print(f"carbono:          n = {n_C:.1f} mol")
print(f"CO2 formado:      {m_CO2:.0f} g = {m_CO2 / 1000:.2f} kg")
print(f"SO2 (se enxofre): {m_SO2:.0f} g = {m_SO2 / 1000:.2f} kg")
print(f"razão de massas SO2/S: {M_SO2 / M_S:.2f}")

# %% md
## 4. Visualização

### 4.1 O dióxido de carbono em 3D — e por que ele não tem ângulo "de água"

O $\mathrm{CO_2}$ é **linear**: o carbono (cinza) no centro, dois oxigênios (vermelho) em lados opostos, $180^\circ$. Compare com a água, angular, do notebook de ligações: mesma família de molécula pequena, geometria completamente diferente.

# %% code title="CO2: ângulo O–C–O medido"
visor, angulo = geometria_3d("O=C=O", atomo=1, cor="#c44e52")
print(f"ângulo O–C–O medido: {angulo:.1f}°   (esperado: 180°)")
visor.show()

# %% md
### 4.2 O dióxido de enxofre: o "V" que anuncia o $sp^2$ com par livre

O $\mathrm{SO_2}$ (o gás do cheiro de fósforo queimado, e o que forma a chuva ácida) é **angular**, com ângulo perto de $119^\circ$: o enxofre tem um par de elétrons não ligante além das duas ligações. Gire e confira o rótulo.

# %% code title="SO2: ângulo O–S–O medido"
visor, angulo = geometria_3d("O=S=O", atomo=1, cor="#e6a700", campo="UFF")
print(f"ângulo O–S–O medido: {angulo:.1f}°   (esperado: ~119°)")
visor.show()

# %% md
### 4.3 Óxidos ácidos em 3D

O $\mathrm{SO_3}$ é o óxido que, com água, forma o ácido sulfúrico — o vilão da chuva ácida. Note o enxofre rodeado por **três** oxigênios, plano, a $120^\circ$.

> **Detalhe técnico:** a geometria aqui é otimizada com o campo de força **UFF** (`campo="UFF"`), porque o MMFF (o padrão) achata o $\mathrm{SO_3}$ e distorce o ângulo do $\mathrm{SO_2}$ — com o UFF o desenho sai com os $120^\circ$ esperados.

# %% code title="SO3 em 3D interativo"
molecula_3d("O=S(=O)=O", rotulo="SO3", altura=340, campo="UFF")

# %% code title="NO2 em 3D interativo"
molecula_3d("[O-][N+]=O", rotulo="NO2", altura=340, campo="UFF")

# No desenho 2D do NO2 aparece um pontinho no nitrogênio: é o elétron não
# emparelhado (o NO2 é um radical). É por isso que ele é tão reativo no ar.

# %% md
### 4.4 Óxidos básicos em 3D

Nos óxidos metálicos o oxigênio aparece como ânion $\mathrm{O^{2-}}$: a ligação é **iônica**, não covalente. Por isso o desenho mostra os íons como esferas separadas — é o retrato do sólido iônico.

# %% code title="CaO (cal virgem) em 3D interativo"
molecula_3d("[Ca+2].[O-2]", rotulo="CaO", largura=480, altura=340)

# %% code title="Na2O em 3D interativo"
molecula_3d("[Na+].[Na+].[O-2]", rotulo="Na2O", altura=340)

# %% code title="Al2O3 (alumina) em 3D interativo"
molecula_3d("[Al+3].[Al+3].[O-2].[O-2].[O-2]", rotulo="Al2O3", altura=340)

# %% md
### 4.5 Os óxidos em 2D, com os risquinhos

Aqui aparece bem a diferença entre **óxido ácido** (ligações covalentes: traços entre C e O) e **óxido básico** (íons). Nos óxidos metálicos, o traço **pontilhado** entre o metal e o oxigênio representa a **atração eletrostática**, não uma ligação covalente — é o mesmo código visual da ligação de hidrogênio.

# %% code title="Galeria 2D dos óxidos"
# mistura covalentes e iônicos: cada célula é desenhada com o estilo certo
figura, eixos = plt.subplots(2, 3, figsize=(10.5, 6.4))
itens = [("O=C=O", "CO$_2$\n(ácido)", False),
         ("O=S=O", "SO$_2$\n(ácido)", False),
         ("O=S(=O)=O", "SO$_3$\n(ácido)", False),
         ("[Ca+2].[O-2]", "CaO\n(básico)", True),
         ("[Na+].[Na+].[O-2]", "Na$_2$O\n(básico)", True),
         ("[C-]#[O+]", "CO\n(neutro)", False)]
for eixo, (smi, titulo, ionico) in zip(eixos.ravel(), itens):
    eixo.imshow(desenho_png(smi, ionico=ionico))
    eixo.set_axis_off()
    eixo.set_title(titulo, fontsize=10, pad=6)
figura.tight_layout()

# %% code title="O monóxido de carbono: o óxido neutro"
estrutura_2d("[C-]#[O+]", legenda="monoxido de carbono (CO)", largura=460, altura=330)

# %% md
### 4.6 Anfoterismo: o mesmo óxido reagindo com ácido **e** com base

Óxidos anfóteros como $\mathrm{Al_2O_3}$ e $\mathrm{ZnO}$ formam sal tanto com ácido quanto com base. É o motivo de a alumina resistir a quase tudo — e de o alumínio ser "protegido" pela camada de óxido.

# %% code title="pH depois de adicionar cada óxido à água"
# valores medidos/típicos de solução saturada ou suspensão aquosa
nomes = ["CaO", "Na$_2$O", "MgO", "Al$_2$O$_3$*", "ZnO*", "CO$_2$", "SO$_2$", "SO$_3$", "CO"]
phs = [12.4, 12.0, 10.4, 7.0, 7.0, 3.9, 2.9, 1.2, 7.0]
familias = ["básico", "básico", "básico", "anfótero", "anfótero",
            "ácido", "ácido", "ácido", "neutro"]
cores = {"básico": "#4c72b0", "ácido": "#c44e52",
         "anfótero": "#dd8452", "neutro": "#888888"}

figura, eixo = plt.subplots(figsize=(8.4, 4.4))
barras = eixo.bar(nomes, phs, color=[cores[f] for f in familias])
eixo.axhline(7, color="#333", linestyle="--", linewidth=1)
for barra, p in zip(barras, phs):
    eixo.annotate(f"{p:.1f}".replace(".", ","),
                  (barra.get_x() + barra.get_width() / 2, p),
                  xytext=(0, 4), textcoords="offset points", ha="center", fontsize=9)
eixo.set_ylim(0, 14)
eixo.set_ylabel("pH da solução/suspensão")
eixo.set_title("Caráter do óxido em água: azul = básico, vermelho = ácido\n"
               "laranja = anfótero (*não dissolve: pH ~ neutro), cinza = neutro")
eixo.tick_params(axis="x", labelrotation=20)
eixo.grid(axis="y", alpha=0.3)
figura.tight_layout()

# %% md
### 4.7 Óxidos e o efeito estufa: concentração de $\mathrm{CO_2}$ na atmosfera

O $\mathrm{CO_2}$ é o óxido ácido mais famoso — e o número que se acompanha. Dados aproximados do Observatório de Mauna Loa (média anual de maio), ppm.

# %% code title="CO2 atmosférico (ppm)"
anos = [1960, 1970, 1980, 1990, 2000, 2010, 2020, 2024]
ppm = [317, 326, 339, 354, 369, 390, 414, 427]

plt.figure(figsize=(7.5, 4.5))
plt.plot(anos, ppm, "o-", color="#c44e52")
for ano, valor in zip(anos, ppm):
    plt.annotate(f"{valor}", (ano, valor), xytext=(0, 7), textcoords="offset points",
                 ha="center", fontsize=8)
plt.xlabel("Ano")
plt.ylabel("CO$_2$ (ppm, partes por milhão)")
plt.title("O óxido ácido que mais se mede: CO$_2$ em Mauna Loa")
plt.grid(alpha=0.3)
plt.show()

# %% md
## 5. Como interpretar o resultado

- **Classificar óxido é olhar quem é o outro elemento:** ametal (ou metal de alta valência como o Mn) → óxido **ácido**; metal alcalino/alcalinoterroso → **básico**; Al, Zn, Pb, Sn → **anfótero**; CO, NO, N$_2$O → **neutro**.
- **A massa "cresce" na combustão** porque o oxigênio entra: $1\ \mathrm{kg}$ de carvão vira $\approx 3{,}7\ \mathrm{kg}$ de $\mathrm{CO_2}$. Não é mágica, é conservação de massa com um reagente gasoso.
- **A proporção é sempre em mol**, depois convertida em massa com a massa molar. Carbono queima 1:1 com $\mathrm{CO_2}$; enxofre 1:1 com $\mathrm{SO_2}$.
- **No 3D**, a geometria conta a história: $\mathrm{CO_2}$ linear ($180^\circ$), $\mathrm{SO_2}$ angular ($\approx 119^\circ$), $\mathrm{SO_3}$ trigonal plano ($120^\circ$).
- **Anfótero não "não reage"**: reage com ácido e com base, mas por isso mesmo não muda o pH da água pura.
- **A chuva ácida** é a soma de $\mathrm{SO_2}$ e $\mathrm{NO_2}$ com água: $\mathrm{H_2SO_3}$, $\mathrm{H_2SO_4}$, $\mathrm{HNO_3}$ — todos óxidos ácidos da tabela.

# %% md
## 6. Exercício de fixação

Um carro a gasolina percorre $100\ \mathrm{km}$ queimando $8{,}0\ \mathrm{kg}$ de combustível. Suponha que o combustível seja **octano** ($\mathrm{C_8H_{18}}$) puro:

$$
\mathrm{2\,C_8H_{18} + 25\,O_2 \rightarrow 16\,CO_2 + 18\,H_2O}
$$

1. Quantos mol de octano há em $8{,}0\ \mathrm{kg}$?
2. Quantos mol de $\mathrm{CO_2}$ se formam (use a proporção da equação)?
3. Qual a massa de $\mathrm{CO_2}$ liberada por quilômetro?
4. Que tipo de óxido é o $\mathrm{CO_2}$? O que acontece com ele na atmosfera úmida?

Tente prever antes de executar a célula seguinte.

# %% code title="Exercício"
MASSA_OCTANO_G = 8000.0
KM = 100.0
SMILES_OCTANO = "CCCCCCCC"

M_oct = massa_molar(SMILES_OCTANO)
n_oct = MASSA_OCTANO_G / M_oct
n_co2 = n_oct * (16 / 2)          # 16 CO2 para cada 2 octanos
m_co2 = n_co2 * M_CO2

print(f"Octano (C8H18):   M = {M_oct:.2f} g/mol")
print(f"1) n de octano:   {n_oct:.2f} mol")
print(f"2) n de CO2:      {n_co2:.2f} mol  (8 x n_octano)")
print(f"3) massa de CO2:  {m_co2:.0f} g = {m_co2 / 1000:.1f} kg")
print(f"   por km:        {m_co2 / KM:.0f} g/km")
print("4) CO2 = óxido ácido: com água forma H2CO3 (chuva levemente ácida, pH ~5,6)")

# %% md
## 7. Resumo

- Óxido: composto binário com oxigênio como segundo elemento.
- **Ácido** (ametal): CO$_2$, SO$_2$, SO$_3$, NO$_2$ → com água dão ácidos.
- **Básico** (metal): Na$_2$O, CaO, MgO, Fe$_2$O$_3$ → com água dão bases.
- **Anfótero**: Al$_2$O$_3$, ZnO, PbO, SnO → reagem com ácido e com base.
- **Neutro**: CO, NO, N$_2$O → não reagem com água.
- Óxido básico + óxido ácido → sal ($\mathrm{CaO + CO_2 \rightarrow CaCO_3}$).
- Estequiometria: converta massa → mol → mol (proporção) → massa, com a massa molar.
- Geometria ajuda a classificar: linear (CO$_2$), angular (SO$_2$), plano trigonal (SO$_3$).

## Referências

- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- FELTRE, R. *Química Geral*. Moderna.
- NOAA Global Monitoring Laboratory. *Trends in Atmospheric Carbon Dioxide* (Mauna Loa).
