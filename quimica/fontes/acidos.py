# -*- coding: utf-8 -*-
"""Rascunho-fonte de quimica/acidos/01_acidos_ionizacao_forca_ph.ipynb.

Gere o .ipynb com:
    python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
        quimica/fontes/acidos.py quimica/acidos/01_acidos_ionizacao_forca_ph.ipynb
"""

# %% md
# Ácidos: ionização, força e pH
### Funções inorgânicas: os ácidos

**Objetivo:** Reconhecer a função ácido (HCl, H$_2$SO$_4$, HNO$_3$, H$_3$PO$_4$, ácidos carboxílicos), escrever a ionização em água, comparar a força dos ácidos pelo $K_a$ e calcular o pH de uma solução a partir da concentração.

### Como usar este notebook

Execute as células na ordem, de cima para baixo. As células Markdown apresentam os conceitos; as células de código fazem as contas e geram as visualizações.

> **Moléculas em 3D:** clique na figura e **arraste o mouse** — ela gira em 360° em qualquer direção. Use a roda do mouse (ou dois dedos no celular) para dar zoom. As cores seguem a convenção dos químicos (C = cinza, H = branco, O = vermelho, N = azul, S = amarelo, Cl = verde, P = laranja).

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

**Ácido** é a substância que, em água, libera íons $\mathrm{H^+}$ (na verdade $\mathrm{H_3O^+}$, o hidrônio). Essa é a definição de **Arrhenius**. Brønsted e Lowry ampliaram a ideia: ácido é quem **doa** um próton.

A ionização é uma reação reversível. Para um ácido genérico $\mathrm{HA}$:

$$
\boxed{\mathrm{HA_{(aq)} \rightleftharpoons H^+_{(aq)} + A^-_{(aq)}}}
\qquad
K_a = \frac{[\mathrm{H^+}]\,[\mathrm{A^-}]}{[\mathrm{HA}]}
$$

Quanto **maior** o $K_a$, mais forte é o ácido. Por isso se usa a escala logarítmica

$$
\mathrm{p}K_a = -\log_{10} K_a
$$

e agora a regra se inverte: **quanto menor o $\mathrm{p}K_a$, mais forte o ácido**. Ácidos fortes ($\mathrm{HCl}$, $\mathrm{HNO_3}$, $\mathrm{H_2SO_4}$) ionizam praticamente 100%: a conta do pH fica direta, $\mathrm{pH} = -\log[\mathrm{ácido}]$. Ácidos fracos (acético, carbônico, hipocloroso) ionizam pouco e precisam do $K_a$.

# %% md
## 2. Exemplo resolvido passo a passo

**Pergunta:** qual é o pH de $250\ \mathrm{mL}$ de solução de $\mathrm{HCl}$ a $0{,}020\ \mathrm{mol/L}$? E o de uma solução de ácido acético de **mesma concentração**? Qual é a massa de $\mathrm{HCl}$ nessa solução?

**Passo 1 — quantidade de matéria.** (volume sempre em litros!)

$$
n = c\,V = 0{,}020 \cdot 0{,}250 = 5{,}0\times10^{-3}\ \mathrm{mol}
$$

**Passo 2 — massa**, com $M_{\mathrm{HCl}} = 36{,}46\ \mathrm{g/mol}$:

$$
m = n\,M = 5{,}0\times10^{-3} \cdot 36{,}46 = \boxed{0{,}18\ \mathrm{g}}
$$

**Passo 3 — pH do ácido forte.** O $\mathrm{HCl}$ ioniza totalmente, então $[\mathrm{H^+}] = 0{,}020$:

$$
\mathrm{pH} = -\log_{10}(0{,}020) = \boxed{1{,}70}
$$

**Passo 4 — pH do ácido fraco.** Para o acético ($K_a = 1{,}8\times10^{-5}$) vale a aproximação $[\mathrm{H^+}] \approx \sqrt{K_a\,c}$:

$$
[\mathrm{H^+}] = \sqrt{1{,}8\times10^{-5} \cdot 0{,}020} = 6{,}0\times10^{-4}\ \mathrm{mol/L}
\qquad
\mathrm{pH} = \boxed{3{,}22}
$$

Mesma concentração, quase **1,5 unidade** de diferença de pH (cerca de 30 vezes menos $\mathrm{H^+}$ no ácido fraco).

# %% md
## 3. Conferindo o cálculo com Python

A próxima célula repete as quatro contas com as grandezas e as unidades do exemplo.

# %% code title="Contas do exemplo"
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

SMILES_HCL = "Cl"
SMILES_ACETICO = "CC(=O)O"
C = 0.020        # mol/L
V = 0.250        # L
KA_ACETICO = 1.8e-5

massa_molar = lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s))
formula = lambda s: rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s))

M_hcl, M_ac = massa_molar(SMILES_HCL), massa_molar(SMILES_ACETICO)
n = C * V
m = n * M_hcl
pH_forte = -math.log10(C)
H_fraco = math.sqrt(KA_ACETICO * C)
pH_fraco = -math.log10(H_fraco)

print(f"Fórmulas:          {formula(SMILES_HCL)} e {formula(SMILES_ACETICO)}")
print(f"Massa molar HCl:   {M_hcl:.2f} g/mol")
print(f"Quantidade:        {n:.4f} mol")
print(f"Massa de HCl:      {m:.3f} g")
print(f"pH do HCl:         {pH_forte:.2f}")
print(f"[H+] do acético:   {H_fraco:.2e} mol/L")
print(f"pH do acético:     {pH_fraco:.2f}")
print(f"Razão [H+]:        {C / H_fraco:.0f} vezes mais H+ no ácido forte")

# %% md
## 4. Visualização

### 4.1 O ácido sulfúrico em 3D — arraste para girar

O $\mathrm{H_2SO_4}$ é o ácido mais usado na indústria (e o que mais aparece em prova). Repare no enxofre (amarelo) no centro de um **tetraedro**: dois oxigênios com ligação dupla e duas hidroxilas com ligação simples. São justamente as duas hidroxilas que perdem $\mathrm{H^+}$ — por isso ele é **diprótico**.

# %% code title="H2SO4 em 3D interativo"
molecula_3d("OS(=O)(=O)O", rotulo="H2SO4")

# %% md
### 4.2 A geometria medida: ângulo O–S–O

As setas mostram as ligações que saem do enxofre; o rótulo é o ângulo medido na estrutura otimizada. Se o ângulo ficar perto de $109{,}5^\circ$, aparece o tetraedro que a fórmula de Lewis prevê para 4 pares de elétrons.

> Nos compostos de enxofre o campo de força MMFF (o padrão) às vezes erra o ângulo por alguns graus. O valor do livro é $\approx 109^\circ$; o que importa aqui é ver o **tetraedro**, não o número exato.

# %% code title="Ângulo O–S–O medido"
visor, angulo = geometria_3d("OS(=O)(=O)O", atomo=1)
print(f"ângulo O–S–O medido: {angulo:.1f}°")
visor.show()

# %% md
### 4.3 O ácido nítrico e o fosfórico também

O ácido nítrico ($\mathrm{HNO_3}$) é monoprótico e forte; o fosfórico ($\mathrm{H_3PO_4}$) é triprótico e fraco (é o acidulante dos refrigerantes tipo cola). Gire cada um: no nítrico há um $\mathrm{N}$ central ligado a três oxigênios, um deles com o $\mathrm{H}$.

# %% code title="HNO3 em 3D interativo"
molecula_3d("O[N+](=O)[O-]", rotulo="HNO3", altura=330)

# %% code title="H3PO4 em 3D interativo"
molecula_3d("OP(=O)(O)O", rotulo="H3PO4", altura=330)

# %% md
### 4.4 Os ácidos em 2D, com os risquinhos das ligações

Cada traço é uma ligação: **um** traço = simples, **dois** = dupla. Os hidrogênios ácidos aparecem ligados ao oxigênio (é o $\mathrm{H}$ do $-OH$), porque é ele que sai na ionização; o $\mathrm{H}$ do carbono fica implícito, como manda a convenção do esqueleto.

# %% code title="Galeria 2D dos ácidos"
galeria_2d(
    ["Cl", "OS(=O)(=O)O", "O[N+](=O)[O-]",
     "OP(=O)(O)O", "CC(=O)O", "OC(=O)O"],
    ["HCl\n(clorídrico)", "H$_2$SO$_4$\n(sulfúrico)", "HNO$_3$\n(nítrico)",
     "H$_3$PO$_4$\n(fosfórico)", "CH$_3$COOH\n(acético)", "H$_2$CO$_3$\n(carbônico)"],
)

# %% md
### 4.5 O ácido acético em 2D, sozinho e maior

No ácido acético o grupo $-COOH$ (carboxila) é o responsável pela acidez. O $\mathrm{H}$ que sai é o do oxigênio, **não** o do carbono: quem sai leva o par de elétrons da ligação O–H.

# %% code title="Estrutura 2D do ácido acético"
estrutura_2d("CC(=O)O", legenda="acido acetico (CH3COOH)", largura=560, altura=380)

# %% md
### 4.6 Comparando a força dos ácidos pelo $\mathrm{p}K_a$

Neste gráfico, **barra para a esquerda = ácido mais forte**. Note que os ácidos "fortes" têm $\mathrm{p}K_a$ negativo: valores aproximados, medidos em água.

# %% code title="pKa dos ácidos comuns"
nomes = ["HCl", "H$_2$SO$_4$ (1º H)", "HNO$_3$", "H$_3$PO$_4$ (1º H)",
         "CH$_3$COOH", "H$_2$CO$_3$", "HClO"]
pkas = [-6.0, -3.0, -1.4, 2.15, 4.76, 6.35, 7.50]
cores = ["#c44e52" if p < 0 else "#4c72b0" for p in pkas]

figura, eixo = plt.subplots(figsize=(8, 4.4))
eixo.barh(nomes[::-1], pkas[::-1], color=cores[::-1])
for y, p in enumerate(pkas[::-1]):
    # rótulo DENTRO da barra, na ponta: não encosta no eixo dos nomes nem no 0
    eixo.annotate(f"{p:+.2f}".replace(".", ","), (p, y),
                  xytext=(6 if p < 0 else -6, 0), textcoords="offset points",
                  va="center", ha="left" if p < 0 else "right",
                  fontsize=9, color="white")
eixo.axvline(0, color="#333", linewidth=1)
eixo.set_xlabel("p$K_a$  (menor = mais forte)")
eixo.set_title("Força dos ácidos: vermelho = forte, azul = fraco")
eixo.grid(axis="x", alpha=0.3)
figura.tight_layout()

# %% md
### 4.7 pH conforme a concentração: forte × fraco

# %% code title="pH x concentração"
concentracoes = np.logspace(-4, 0, 80)              # de 0,0001 a 1 mol/L
ph_forte = -np.log10(concentracoes)
ph_fraco = -np.log10(np.sqrt(KA_ACETICO * concentracoes))

plt.figure(figsize=(7.5, 4.5))
plt.plot(concentracoes, ph_forte, color="#c44e52", label="HCl (forte: ioniza 100%)")
plt.plot(concentracoes, ph_fraco, color="#4c72b0", label="CH$_3$COOH (fraco: ioniza pouco)")
plt.scatter([C], [pH_forte] if C > 0.01 else [pH_forte], color="#c44e52", zorder=5)
plt.xscale("log")
plt.xlabel("Concentração (mol/L)")
plt.ylabel("pH")
plt.title("Mesma concentração, pH diferente: a força do ácido importa")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.show()

# %% md
## 5. Como interpretar o resultado

- **pH é logarítmico:** a diferença de $1{,}5$ unidade entre $\mathrm{HCl}$ e acético a $0{,}020\ \mathrm{mol/L}$ significa cerca de $30$ vezes menos $\mathrm{H^+}$ livre.
- **Ácido forte não tem $K_a$ útil:** a ionização é total, então a conta é só $-\log c$. Só use $\sqrt{K_a c}$ para ácido fraco.
- **Número de hidrogênios ácidos ≠ número de H da fórmula.** O acético tem 4 H e apenas 1 é ácido; o fosfórico tem 3 H e os três são ácidos (mas cada um mais fraco que o anterior).
- **No 3D**, a geometria do enxofre e do fósforo é tetraédrica, e a do carbono da carboxila é trigonal plana. Gire as figuras e compare.
- O erro clássico em prova é esquecer de converter o volume para litros: $250\ \mathrm{mL} = 0{,}250\ \mathrm{L}$.

# %% md
## 6. Exercício de fixação

Uma solução de ácido sulfúrico $\mathrm{H_2SO_4}$ foi preparada com $0{,}0050\ \mathrm{mol/L}$, considerando a **primeira** ionização completa.

1. Qual o pH dessa solução?
2. Qual a massa de $\mathrm{H_2SO_4}$ em $2{,}0\ \mathrm{L}$ dessa solução?
3. Se a solução fosse de ácido acético, o pH seria maior ou menor?

Tente prever antes de executar a célula seguinte.

# %% code title="Exercício"
C_EX, V_EX = 0.0050, 2.0
M_EX = massa_molar("OS(=O)(=O)O")

n_ex = C_EX * V_EX
m_ex = n_ex * M_EX
pH_ex = -math.log10(C_EX)
pH_ac_ex = -math.log10(math.sqrt(KA_ACETICO * C_EX))

print(f"1) pH do H2SO4:               {pH_ex:.2f}")
print(f"2) massa em 2,0 L:            {m_ex:.2f} g   (n = {n_ex:.4f} mol)")
print(f"3) pH se fosse acético:       {pH_ac_ex:.2f}  -> "
      f"{'maior' if pH_ac_ex > pH_ex else 'menor'} (ácido fraco libera menos H+)")

# %% md
## 7. Resumo

- Ácido, em água, libera $\mathrm{H^+}$; doa próton (Brønsted–Lowry).
- Ionização: $\mathrm{HA \rightleftharpoons H^+ + A^-}$, com $K_a$ e $\mathrm{p}K_a = -\log K_a$.
- **Ácido forte:** ioniza ~100% → $\mathrm{pH} = -\log c$. **Ácido fraco:** $[\mathrm{H^+}] \approx \sqrt{K_a c}$.
- Strong list comum: $\mathrm{HCl}$, $\mathrm{HBr}$, $\mathrm{HI}$, $\mathrm{HNO_3}$, $\mathrm{H_2SO_4}$ (1º H), $\mathrm{HClO_4}$.
- Fracos comuns: $\mathrm{CH_3COOH}$, $\mathrm{H_2CO_3}$, $\mathrm{H_3PO_4}$, $\mathrm{HClO}$, $\mathrm{HF}$.
- Volume de solução sempre em **litros** na fórmula $n = cV$.

## Referências

- ATKINS, P.; JONES, L. *Princípios de Química*. Bookman.
- IUPAC. *Compendium of Chemical Terminology* (Gold Book).
- FELTRE, R. *Química Geral*. Moderna.
