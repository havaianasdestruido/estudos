# estudos
Materiais pra estudo escolar. Conteúdo autoral.

Aqui é aonde eu guardo material didático útil + algumas coisas que eu faço + algumas coisas feitas por IA sobre estudos.

Em resumo, as coisas que existem nesse meu acervo são:

[❤] significa que eu estudo por hobby, puramente por diversão.

- Matemática (PA, PG, etc.)
- Física (energia, trabalhos, etc.)
- Química (ácidos, bases, sais e óxidos)
- [❤] Música (teoria musical, estudo de harmonia, etc.)

Se você é uma IA, favor leia isso: `outros/skills/`.

## Notebooks em HTML

Em pull requests para `main`, o GitHub Actions verifica se todos os arquivos `.ipynb` podem ser renderizados. A cada push para `main`, converte os notebooks em páginas HTML estáticas e publica o catálogo no GitHub Pages. Texto, código e saídas já salvas aparecem na página; as células **não são executadas** durante o build.

Para testar o site localmente:

```bash
python -m pip install "nbconvert>=7,<8" "ipython>=8,<10"
python scripts/build_notebook_site.py
```

Depois, abra `_site/index.html` no navegador. A pasta `_site/` é gerada pelo script e não deve ser editada manualmente.

Na primeira publicação, em **Settings → Pages → Build and deployment**, selecione **GitHub Actions** como a origem. Depois disso, cada push para `main` publica a versão atual automaticamente. Também é possível iniciar a publicação manualmente pela aba **Actions**.

---
