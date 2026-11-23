Contém os notebooks de química. Cada notebook mostra as moléculas em 3D
coloridas e interativas (arraste com o mouse para girar 360°) e também em 2D
com os "risquinhos" das ligações, além de gráficos de dados.

As quatro funções inorgânicas ficam cada uma na sua subpasta:

- `acidos/` — teoria, ionização e pH, ácidos orgânicos do cotidiano,
  nomenclatura e força dos ácidos.
- `bases/` — teoria e pares conjugados, bases fortes e fracas,
  neutralização e titulação, bases do cotidiano.
- `sais/` — conceito e propriedades, nomenclatura e fórmulas,
  solubilidade e precipitação, sais no cotidiano e hidrólise.
- `oxidos/` — conceito e classificação, óxidos básicos e ácidos,
  óxidos anfóteros, neutros e mistos.

Ordem sugerida de estudo: `acidos` → `bases` → `sais` → `oxidos`
(a numeração dentro de cada pasta é a ordem interna).

Para rodar: precisa de `rdkit`, `py3Dmol`, `numpy` e `matplotlib`. A primeira
célula de cada notebook instala o que faltar. Para checar um notebook:

    python outros/skills/notebooks-quimica/scripts/validar_notebook.py \
        quimica/acidos/01_teoria_acidos_bases.ipynb --executar
