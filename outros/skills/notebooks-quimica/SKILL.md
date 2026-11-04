---
name: notebooks-quimica
description: >
  Como criar notebooks de química (.ipynb) para estudo, no mesmo estilo dos
  notebooks de física/matemática deste repositório. Define o padrão visual:
  toda molécula aparece em 3D colorida e INTERATIVA (arrastar o mouse gira em
  360°), mais estruturas 2D "com risquinhos" (ligações simples/duplas/triplas,
  anel aromático e ligação de hidrogênio pontilhada), legendas com acento,
  gráficos 2D de dados e figuras de energia. Inclui scripts prontos
  (rdkit + py3Dmol) e dois notebooks-exemplo completos e já executados.
---

# Skill: notebooks de química

Este skill descreve **como fazer um notebook de química** que combine com o resto do
repositório (`fisica/`, `matematica/`): markdown em português explicando o conceito,
células de código curtas com as contas e, no fim, visualizações.

A regra de ouro deste material: **molécula nunca é só uma imagem estática.**

| O que | Como aparece | Biblioteca |
|---|---|---|
| Molécula 3D | bola-e-vareta colorida, **girando com o mouse em 360°** (arrastar; roda do mouse = zoom) | `py3Dmol` (3Dmol.js) |
| Estrutura 2D | traço preto com "risquinhos": simples 1 traço, dupla 2 traços, tripla 3 traços, anel aromático com duplas alternadas, ligação de hidrogênio pontilhada | `rdkit` (`rdMolDraw2D`) |
| Legendas/rótulos | título em matplotlib (acento e subscrito corretos) | `matplotlib` |
| Gráficos 2D | eixos, pontos, barras, legendas — dados numéricos | `matplotlib` |
| Geometria/energia 3D | vetores e ângulos que também giram em 360° | `plotly` (WebGL) |

---

## 1. Estrutura obrigatória do notebook

Siga exatamente a espinha dorsal dos notebooks de física (veja
`fisica/energia/conservacao_energia/14_energia_mecanica.ipynb`):

1. **Célula 1 — título** (markdown): `# Nome da molécula/fenômeno`, subtítulo
   `### Tema maior`, `**Objetivo:** ...`, seção `### Como usar este notebook`
   explicando "execute as células em ordem; as de código geram as visualizações;
   nas moléculas 3D, arraste com o mouse para girar".
2. **`## 1. Conceito principal`** — fórmula em LaTeX `$$...$$` com `\boxed{}`.
3. **`## 2. Exemplo resolvido passo a passo`** — números reais, cada conta em
   um bloco `$$`.
4. **`## 3. Conferindo o cálculo com Python`** — célula curta que reproduz a
   conta e imprime com unidades (`print(f"... {valor:.2f} mol")`).
5. **`## 4. Visualização`** — aqui entram as moléculas: 3D interativa primeiro,
   depois 2D com risquinhos.
6. **`## 5. Como interpretar o resultado`** — leitura do desenho, o que olhar.
7. **`## 6. Exercício de fixação`** — enunciado em markdown, célula seguinte
   resolve; peça para o aluno prever antes de rodar.
8. **`## 7. Resumo`** — tópicos com `-`.
9. **`## Referências`** (opcional mas recomendado) — IUPAC / livros usados.

Metadados do `.ipynb` (copie o padrão já existente):

```json
"metadata": {
  "authors": [{"name": "Material didático ..."}],
  "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
  "language_info": {"name": "python", "version": "3.11.2", "...": "..."},
  "title": "Nome do notebook"
}
```

`nbformat: 4`, `nbformat_minor: 5`.

---

## 2. Instalação das dependências

O notebook precisa rodar no Colab **e** localmente. Use sempre esta célula de
preparação no topo (idempotente — não reinstala se já existir):

```python
# @title Preparar o ambiente (só precisa rodar uma vez)
import importlib.util, sys

def _tem(modulo):
    return importlib.util.find_spec(modulo) is not None

if not _tem("rdkit"):
    %pip install -q rdkit
if not _tem("py3Dmol"):
    %pip install -q py3Dmol
import rdkit, py3Dmol  # noqa: verifica se carregou
print("rdkit", rdkit.__version__, "| py3Dmol", py3Dmol.__version__)
```

Notas:
- No Colab, `%pip` funciona direto. Fora do Colab, troque por
  `subprocess.run([sys.executable, "-m", "pip", "install", "-q", "rdkit"])`.
- `py3Dmol` carrega o **3Dmol.js por CDN** (`cdn.jsdelivr.net`). Sem internet o
  bloco 3D mostra um aviso vermelho; por isso use `mostrar_molecula(...)` (seção
  3.5), que cai automaticamente para o 2D quando o `py3Dmol` não está instalado.
- Não use `mol.render()` de `rdkit.Chem.Draw.IPythonConsole`: ele gera uma imagem
  **estática**. A interatividade vem do `py3Dmol`.

---

## 3. Receita: molécula 3D colorida e interativa (o principal)

Código mínimo, testado:

```python
# @title Molécula 3D interativa — arraste com o mouse para girar 360°
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def molecula_3d(smiles, largura=640, altura=420, estilo="ball-and-stick"):
    """SMILES -> visualização 3D interativa (arrastar = girar, rolar = zoom)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    mol = Chem.AddHs(mol)                       # hidrogênios explícitos: aparecem na figura

    params = AllChem.ETKDGv3()                  # gera a geometria 3D
    params.randomSeed = 0xC0FFEE                # mesmo desenho toda vez que rodar
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("não foi possível gerar a geometria 3D")
    AllChem.MMFFOptimizeMolecule(mol)           # relaxa a estrutura (campo de força)

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.28, "colorscheme": "Jmol"}})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    return visor

molecula_3d("CC(=O)O")   # ácido acético
```

Por que cada linha:
- `AddHs` — sem os hidrogênios a geometria 3D fica "torta" (o ângulo do carbono
  não fecha) e o aluno não vê o que a fórmula mostra. Sempre adicione H.
- `ETKDGv3` + `randomSeed` fixo — reprodutibilidade: o notebook gera o mesmo
  desenho em qualquer máquina. **Sempre** fixe a semente.
- `MMFFOptimizeMolecule` — sem ele as ligações aparecem esticadas.
- `colorscheme: "Jmol"` — o esquema de cores por elemento que os químicos esperam:
  **C = cinza, H = branco, O = vermelho, N = azul, S = amarelo, P = laranja,
  Cl = verde, F = verde-claro, Br = marrom-escuro, I = roxo**.
- `stick.radius` 0.10–0.16 e `sphere.scale` 0.22–0.32 dão o clássico
  "ball-and-stick". `scale` muito alto vira "space-filling" e esconde as ligações.
- `zoomTo()` **antes** de `render`; se você adicionar rótulos/setas depois do
  `zoomTo`, chame `viewer.render()` de novo (o py3Dmol faz isso sozinho no último
  comando, mas ao montar a cena em etapas prefira `zoomTo` no final).

Interatividade já vem de graça no 3Dmol.js: **arrastar = rotação 360°**
(horizontal e vertical), **roda/pinça = zoom**, **botão direito ou 2 dedos = pan**.
Não precisa de JavaScript extra. Coisas úteis:

```python
visor.spin("y", 0.6)          # giro automático (efeito "vitrine"); use com moderação
visor.spin(False)             # desliga
visor.setStyle({"stick": {"radius": 0.13, "colorscheme": "Jmol"}})  # só varetas
visor.setStyle({"sphere": {"scale": 0.9, "colorscheme": "Jmol"}})   # só bolas
visor.addLabel("ácido acético", {"position": {"x": 0, "y": 0, "z": 0},
                                 "fontColor": "black", "backgroundColor": "white",
                                 "backgroundOpacity": 0.7})
```

### 3.1 Estilos nomeados (use estes, não invente)

| estilo | spec | quando usar |
|---|---|---|
| `"ball-and-stick"` | `stick 0.13` + `sphere 0.28` | padrão para moléculas pequenas: mostra ligações **e** volume atômico |
| `"stick"` | `stick 0.15` | moléculas grandes / muitas ligações |
| `"space-filling"` | `sphere 0.95` | volume molecular, "como a molécula ocupa espaço" |
| `"wireframe"` | `line 2.0` | visão esquemática rápida |

### 3.2 Comparar moléculas lado a lado (interativas, cada uma gira sozinha)

Use uma célula por molécula; nunca tente pôr duas moléculas no mesmo `view`
(mesmo com `viewergrid`, o zoom fica ruim). Para "galeria", prefira o 2D
(`MolsToGridImage`, seção 4).

```python
for nome, smiles in [("etanol", "CCO"),
                     ("ácido acético", "CC(=O)O"),
                     ("acetato de etila", "CC(=O)OCC")]:
    print(f"{nome} — {smiles}")
    display(molecula_3d(smiles, altura=360).show())   # `display` em célula de loop
```

### 3.3 Geometria molecular: ângulos e vetores em 3D interativo

Para explicar 109,5° / 120° / 104,5°, desenhe vetores e rótulos que giram junto:

```python
# @title Geometria do carbono (tetraédrico) em 3D interativo
import numpy as np, py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

def geometria_3d(smiles, atom_idx=0, largura=640, altura=460):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    centro = np.array(conf.GetAtomPosition(atom_idx))
    vizinhos = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx).GetNeighbors()]

    visor = py3Dmol.view(width=largura, height=altura)
    visor.addModel(Chem.MolToMolBlock(mol), "mol")
    visor.setStyle({"stick": {"radius": 0.11, "colorscheme": "Jmol"},
                    "sphere": {"scale": 0.25, "colorscheme": "Jmol"}})
    for i in vizinhos:                                   # seta do centro até cada vizinho
        p = np.array(conf.GetAtomPosition(i))
        visor.addArrow({"start": {"x": float(centro[0]), "y": float(centro[1]), "z": float(centro[2])},
                        "end":   {"x": float(p[0]),      "y": float(p[1]),      "z": float(p[2])},
                        "radius": 0.045, "radiusRatio": 1.8,
                        "color": "#1f6feb", "mid": 0.85})
    a, b = vizinhos[0], vizinhos[1]                      # rótulo do ângulo
    va, vb = (np.array(conf.GetAtomPosition(a)) - centro,
              np.array(conf.GetAtomPosition(b)) - centro)
    ang = np.degrees(np.arccos(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb))))
    meio = centro + 0.55 * (va + vb) / 2
    visor.addLabel(f"{ang:.1f}°", {"position": {"x": float(meio[0]), "y": float(meio[1]),
                                                "z": float(meio[2])},
                                   "fontColor": "#1f6feb", "fontSize": 14,
                                   "backgroundOpacity": 0.0})
    visor.setBackgroundColor("#ffffff")
    visor.zoomTo()
    return visor, ang

visor, angulo = geometria_3d("C", atom_idx=0)   # metano
print(f"ângulo H–C–H medido: {angulo:.2f}°")
visor.show()
```

> `addArrow` do 3Dmol.js recebe `radius` pequeno (0.03–0.06) para não parecer
> tubeira; `mid` entre 0 e 1 é onde a ponta começa.

### 3.4 Molécula "real" a partir de arquivo (SDF/PDB) em vez de SMILES

Para moléculas que já vêm de um arquivo (cristalografia, base de dados):

```python
mol = Chem.MolFromMolFile("cafeina.mol", removeHs=False)   # ou MolFromMolBlock(texto)
print(Chem.MolToSmiles(mol))                                # confira se leu certo!
visor = molecula_3d(Chem.MolToSmiles(mol))                  # reaproveita a receita
```

Sempre reimprima o SMILES: é o teste rápido de que o arquivo foi lido como você
esperava (ordem de átomos, hidrogênios, cargas).

### 3.5 Função "mostrar" com plano B (use esta no notebook final)

```python
def mostrar_molecula(smiles, largura=640, altura=420):
    """3D interativo se o py3Dmol existir; senão, o 2D do RDKit."""
    try:
        return molecula_3d(smiles, largura=largura, altura=altura)
    except ImportError:
        print("py3Dmol não instalado — mostrando a estrutura 2D no lugar do 3D.")
        return estrutura_2d(smiles, largura=largura, altura=altura)

mostrar_molecula("CC(=O)O")
```

Assim o notebook nunca "quebra" no meio da aula por causa de uma dependência.

---

## 4. Receita: estrutura 2D "com risquinhos"

"Risquinhos" na prática significa: **um traço por ligação simples, dois por
ligação dupla, três por tripla, pontilhado para ligação de hidrogênio / força
intermolecular, e duplas alternadas no anel aromático**. Quem desenha isso é o
RDKit (o anel aromático sai no estilo Kekulé — é o desenho padrão de vestibular;
não existe opção de "círculo no meio" na API de desenho, e o tracejado se faz com
uma ligação de ordem zero, seção 4.4).

### 4.1 Uma molécula, 2D, com legenda

```python
# @title Estrutura 2D
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

rdDepictor.SetPreferCoordGen(True)   # layout 2D mais "de livro didático"

def estrutura_2d(smiles, largura=520, altura=360, legenda=None, mostrar_h=False):
    mol = Chem.MolFromSmiles(smiles)
    if mostrar_h:
        mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    opts = desenho.drawOptions()
    opts.bondLineWidth = 2.2          # espessura do "risquinho"
    opts.addStereoAnnotation = True   # marca R/S quando houver carbono quiral
    desenho.DrawMolecule(mol, legend=legenda or "")
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())   # SVG = nítido em qualquer zoom

estrutura_2d("CC(=O)O", legenda="acido acetico")   # legenda do RDKit: sem acento (seção 4.5)
```

Regra do "risquinho": **não desenhe hidrogênio ligado a carbono** (é o padrão
esqueleto; o H do –OH sim, porque muda a função química). Use `mostrar_h=True`
só quando o hidrogênio for o assunto da aula (ligação de hidrogênio, quiralidade).

### 4.2 Galeria de moléculas (compare funções orgânicas)

```python
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG

mols = [Chem.MolFromSmiles(s) for s in
        ["CCO", "CC=O", "CC(=O)O", "CC(=O)OCC", "c1ccccc1", "CC(=O)N"]]
legendas = ["etanol (alcool)", "etanal (aldeido)", "acido acetico (acido)",
            "acetato de etila (ester)", "benzeno (aromatico)", "acetamida (amida)"]
imagem = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(320, 260),
                              legends=legendas, useSVG=True)
imagem if hasattr(imagem, "data") else SVG(imagem)   # funciona com e sem IPython
```

> ⚠️ **Atenção ao acento na legenda** (testado, é o erro mais comum deste skill):
> o motor de desenho do RDKit corrompe caracteres fora do ASCII — `"ácido acético"`
> sai como `"cido ac tico"` e `"H₂O"` perde o subscrito. Então:
> **legenda dentro do desenho do RDKit → ASCII** (`legenda_segura`);
> **título bonito, com acento e subscrito → matplotlib** (seção 4.6).

### 4.3 Ligação de hidrogênio com tracejado (os "risquinhos" entre moléculas)

O truque: criar uma ligação de **ordem zero** (`Chem.BondType.ZERO`) entre o H de
uma molécula e o O da outra. O RDKit desenha pontilhado — é exatamente a
convenção de ligação de hidrogênio.

```python
# @title Ligação de hidrogênio entre duas águas
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

def agua_dimer():
    rw = Chem.RWMol(Chem.MolFromSmiles("O.O"))     # dois O, H implícitos por enquanto
    rw = Chem.RWMol(Chem.AddHs(rw.GetMol()))       # agora os H existem para poder ligar
    # descobre índices de um O e do H ligado ao outro O
    oxigenios = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "O"]
    hidrogenios = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "H"]
    ligacao_h = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rw.GetBonds()]
    h_do_segundo_o = [h for o, h in ligacao_h if o == oxigenios[1]][0]
    rw.AddBond(h_do_segundo_o, oxigenios[0], Chem.BondType.ZERO)   # pontilhado
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    rdDepictor.Compute2DCoords(mol)

    desenho = rdMolDraw2D.MolDraw2DSVG(560, 300)
    desenho.drawOptions().bondLineWidth = 2.0
    desenho.DrawMolecule(mol, legend="ligação de hidrogênio (pontilhado)")
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

agua_dimer()
```

### 4.4 Realce de um grupo (destaque colorido no 2D)

Útil para "olhe o grupo funcional": `highlightAtoms` pinta o átomo e as ligações
que chegam nele.

```python
def estrutura_realcada(smiles, padrao_smarts, largura=520, altura=340):
    """2D com o grupo funcional destacado em laranja (padrão SMARTS)."""
    mol = Chem.MolFromSmiles(smiles)
    padrao = Chem.MolFromSmarts(padrao_smarts)
    corres = mol.GetSubstructMatch(padrao)
    rdDepictor.Compute2DCoords(mol)
    desenho = rdMolDraw2D.MolDraw2DSVG(largura, altura)
    desenho.DrawMolecule(mol, highlightAtoms=list(corres),
                         highlightAtomColors={i: (0.99, 0.75, 0.30) for i in corres})
    desenho.FinishDrawing()
    return SVG(desenho.GetDrawingText())

estrutura_realcada("CC(=O)OCC", "C(=O)O")   # destaca o grupo éster/carboxilato
```

### 4.5 Legenda correta: `legenda_segura` + figura com título em matplotlib

Como o desenho do RDKit corrompe acento/subscrito, o padrão do repositório é:

```python
import io, unicodedata
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

def legenda_segura(texto):
    """ASCII para o desenho do RDKit (que corrompe acentos e subscritos)."""
    tabela = str.maketrans("₀₁₂₃₄₅₆₇₈₉°", "0123456789o")
    texto = unicodedata.normalize("NFKD", texto.translate(tabela))
    return texto.encode("ascii", "ignore").decode("ascii")

def figura_com_legenda(smiles, titulo, largura=520, altura=360, mostrar_h=False):
    """Molécula 2D + título em matplotlib (aceita acento e mathtext)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if mostrar_h:
        mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol)

    desenho = rdMolDraw2D.MolDraw2DCairo(largura, altura)
    desenho.drawOptions().bondLineWidth = 2.2
    desenho.drawOptions().padding = 0.12          # não corta "-OH" na borda
    desenho.DrawMolecule(mol)                     # sem legenda: ela vem do matplotlib
    desenho.FinishDrawing()
    imagem = np.asarray(Image.open(io.BytesIO(desenho.GetDrawingText())))

    figura, eixo = plt.subplots(figsize=(largura / 100, (altura + 46) / 100))
    eixo.imshow(imagem)
    eixo.set_axis_off()
    eixo.set_title(titulo, fontsize=11, pad=8)
    figura.tight_layout()
    return figura

figura_com_legenda("CC(=O)O", "Ácido acético (CH$_3$COOH)")   # $...$ = mathtext
```

Resultado: imagem nítida (Cairo), legenda com **"Ácido acético (CH₃COOH)"** de
verdade, e `°`, `–`, `≡` funcionando no título. Se quiser a imagem em PNG/SVG
para o markdown, use `figura.savefig("acido.png", dpi=110, bbox_inches="tight")`.

Para uma grade com títulos corretos, use `galeria_com_legendas(...)` do
`scripts/rotulos.py` (mesma ideia: uma imagem RDKit por célula + título do
matplotlib).

---

## 5. Receita: gráficos 2D de dados (o padrão da casa)

Mesmo estilo dos notebooks de física: `matplotlib`, `figsize`, grade suave,
eixos com unidade, título em português, legenda só quando há mais de uma série.

```python
# @title Energia de ligação por molécula
import matplotlib.pyplot as plt

moleculas = ["H–H", "O=O", "N≡N", "C–H"]
energias = [436, 498, 946, 413]   # kJ/mol
plt.figure(figsize=(7.5, 4.5))
plt.bar(moleculas, energias, color="#4c72b0")
for x, e in zip(moleculas, energias):
    plt.annotate(f"{e}", (x, e), xytext=(0, 4), textcoords="offset points", ha="center")
plt.ylabel("Energia de ligação (kJ/mol)")
plt.title("Quanto mais forte a ligação, mais energia para quebrá-la")
plt.grid(axis="y", alpha=0.3)
plt.show()
```

Outros gráficos que combinam com química de ensino médio:
- **Curva de energia potencial × distância** de ligação (mostra o poço e a
  distância de equilíbrio) — `plt.plot(distancia, energia)`.
- **Titulação**: pH × volume de titulante.
- **Lei de Beer–Lambert**: absorbância × concentração (reta).
- **Cinética**: concentração × tempo para ordem 0, 1, 2 (compare com
  `matematica/PG` e `PA`, que já têm o padrão de séries).
- **Balanceamento**: barras empilhadas de átomos antes/depois.

### 5.1 Origem 3D com vetores (plotly) — quando quiser medir ângulo "de verdade"

`plotly` dá rotação 360° em WebGL e é a melhor opção para mostrar vetores de
ligação, orbitais ou planos de simetria sem depender do CDN do 3Dmol:

```python
# @title Vetores de ligação em 3D (arraste para girar)
import numpy as np, plotly.graph_objects as go

def vetores_3d(vetores, rotulos=None, titulo="Vetores de ligação"):
    fig = go.Figure()
    for i, (nome, v) in enumerate(zip(rotulos or [""] * len(vetores), vetores)):
        v = np.asarray(v, dtype=float)
        fig.add_trace(go.Scatter3d(x=[0, v[0]], y=[0, v[1]], z=[0, v[2]],
                                   mode="lines", line=dict(color="#4c72b0", width=6),
                                   hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter3d(x=[v[0]], y=[v[1]], z=[v[2]], mode="markers+text",
                                   marker=dict(size=4, color="#4c72b0"),
                                   text=[nome], textposition="top center",
                                   showlegend=False))
    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers",
                               marker=dict(size=6, color="#222"), showlegend=False))
    fig.update_layout(title=titulo, scene=dict(xaxis_title="x (Å)", yaxis_title="y (Å)",
                                                zaxis_title="z (Å)", aspectmode="cube"),
                      margin=dict(l=0, r=0, t=40, b=0))
    fig.show()

vetores_3d([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
           ["H1", "H2", "H3", "H4"], "Metano: quatro ligações a 109,5°")
```

Para salvar figuras em vez de exibir interativo (relatório/PDF):
`fig.write_image("...png")` (precisa de `kaleido`) ou
`fig.write_html("...html", include_plotlyjs="cdn")`.

---

## 6. Tabela de cores e convenções (siga sempre)

| Elemento | Cor Jmol | Hex aproximado | Observação |
|---|---|---|---|
| H | branco | `#FFFFFF` | em fundo branco use contorno/`stick`; bola branca "desaparece" |
| C | cinza-escuro | `#909090` | nunca pinte carbono de preto puro em 3D |
| N | azul | `#3050F8` | |
| O | vermelho | `#FF0D0D` | |
| S | amarelo | `#FFFF30` | |
| P | laranja | `#FF8000` | |
| Cl | verde | `#1FF01F` | |
| F | verde-claro | `#90E050` | |
| Br | marrom | `#A62929` | |
| I | roxo | `#940094` | |
| metais | variados | — | use `colorscheme: "Jmol"` e confira no desenho |

Outras convenções:
- **Fundo branco** nos dois tipos (`#ffffff`), para imprimir e para o notebook
  ficar igual aos gráficos de matplotlib.
- **Não use gradiente nem "comicMode"**: o material é didático, não é pôster.
- Legendas em português, com a **fórmula molecular** entre parênteses:
  `"ácido acético (CH₃COOH)"`.
- Unidades sempre: `Å` para distâncias, `°` para ângulos, `kJ/mol` para energia.
- Se a molécula tiver **carbono quiral**, ligue `opts.addStereoAnnotation = True`
  e escreva no markdown qual carbono é (ex.: "o carbono 2 é quiral").

---

## 7. Checklist antes de dizer que o notebook está pronto

- [ ] Primeira célula tem título, objetivo e "Como usar".
- [ ] A célula de preparação do ambiente roda do zero (Colab) sem erro.
- [ ] Toda molécula tem **uma** visualização 3D colorida e interativa
      (`py3Dmol`), **não** só imagem estática.
- [ ] As moléculas pequenas usam `ball-and-stick` (dá para ver ligação *e* bola).
- [ ] Há pelo menos uma estrutura **2D com risquinhos** da molécula-chave.
- [ ] Todo `EmbedMolecule` tem `randomSeed` fixo.
- [ ] Todo `MolFromSmiles` tem checagem de `None` ou `try/except` com mensagem clara.
- [ ] Há pelo menos um gráfico 2D de dados (não só desenhos de molécula).
- [ ] Exercício e resumo no fim; contas com unidade nos `print`.
- [ ] Texto em português do Brasil, com vírgula decimal nas contas
      (`7{,}5`) e ponto no Python (`7.5`) — cuidado para não misturar.
- [ ] O notebook executa de cima a baixo sem erro (`Kernel > Restart & Run All`).
- [ ] Não há saída gigante salva (SVG/PNG com mais de ~500 KB): reduza a
      `subImgSize` ou não salve a célula 3D executada (a cena 3D deve ser
      gerada, não salva em base64 no `.ipynb`).

---

## 8. Erros comuns e como resolver

| Sintoma | Causa | Solução |
|---|---|---|
| `ModuleNotFoundError: rdkit` | ambiente sem RDKit | rode a célula da seção 2 (`pip install rdkit`) |
| 3D mostra caixa vermelha "3Dmol.js failed to load" | sem internet / CDN bloqueado | avise no markdown e mantenha o 2D como fallback |
| Molécula 3D "torta"/estranha | faltou `AddHs` ou `MMFFOptimize` | `Chem.AddHs` antes do `EmbedMolecule` |
| Geometria muda a cada execução | semente aleatória | `params.randomSeed = 0xC0FFEE` |
| `ValueError: SMILES inválido` | erro de digitação (`C(=O` sem fechar) | teste o SMILES no 2D primeiro |
| 2D sem hidrogênios no –OH | RDKit omite H implícito | `Chem.AddHs(mol)` e `Compute2DCoords` de novo |
| Ligação dupla "colada" | molécula grande/pequena | ajuste `largura`/`altura` e `opts.fixedBondLength` |
| Legenda cortada no 2D | texto grande em imagem pequena | use SVG (não PNG) e `legenda` curta |
| `addStereoAnnotation` sem efeito | molécula sem estereoquímica definida | escreva o SMILES com `@`/`@@` |

---

## 9. Arquivos deste skill

```
outros/skills/notebooks-quimica/
├── SKILL.md                          # este guia
├── scripts/
│   ├── quimica.py                    # funções prontas: molecula_3d, estrutura_2d, geometria_3d...
│   ├── rotulos.py                    # legenda sem acento + título correto via matplotlib
│   ├── criar_notebook.py             # monta um .ipynb a partir de um script de células
│   └── validar_notebook.py           # confere a estrutura e roda o notebook de ponta a ponta
└── exemplos/
    ├── exemplo_acido_acetico.ipynb   # massa molar, solução, pH; 3D + 2D + gráfico
    ├── exemplo_agua_ligacoes.ipynb   # ângulo 104,5°, ligação de H pontilhada, ebulição
    └── fontes/                       # os .py que geraram os dois notebooks acima
```

- `scripts/quimica.py` — funções prontas (`molecula_3d`, `estrutura_2d`,
  `galeria_2d`, `geometria_3d`, `descricao`, `formula`, `massa_molar`). Copie o
  bloco de código da seção 3 e 4 para dentro do notebook (o material didático
  deve ser autocontido); o arquivo serve como fonte única e para testar antes de
  colar.
- `scripts/rotulos.py` — `legenda_segura()` e `figura_com_legenda()`; use quando
  o título precisar de acento/subscrito bonito (seção 4.5).
- `scripts/criar_notebook.py` — evita editar JSON à mão: você escreve as células
  num `.py` separado por marcadores e ele gera o `.ipynb` com os metadados
  corretos.
- `scripts/validar_notebook.py` — executa o notebook inteiro num kernel
  temporário e lista as células que falharam.

Uso:

```bash
python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
    minha_aula.py outro/notebook.ipynb
python outros/skills/notebooks-quimica/scripts/validar_notebook.py \
    outro/notebook.ipynb
```
