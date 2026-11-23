# Fontes dos notebooks de química

Os arquivos `.py` aqui são os **rascunhos-fonte** dos notebooks de
`quimica/acidos/`, `quimica/bases/`, `quimica/sais/` e `quimica/oxidos/`.

Eles existem porque editar JSON de `.ipynb` à mão é chato e fácil de errar.
Cada arquivo usa os marcadores `# %% md` e `# %% code` (ver o skill
`outros/skills/notebooks-quimica/SKILL.md`) e é convertido com:

```bash
python outros/skills/notebooks-quimica/scripts/criar_notebook.py \
    quimica/fontes/acidos.py quimica/acidos/01_acidos_ionizacao_forca_ph.ipynb \
    --titulo "Ácidos: ionização, força e pH"

python outros/skills/notebooks-quimica/scripts/validar_notebook.py \
    quimica/acidos/01_acidos_ionizacao_forca_ph.ipynb --executar
```

Se você mudar o notebook direto no Jupyter, o `.py` fica desatualizado —
nesse caso o `.ipynb` é que é a fonte.
