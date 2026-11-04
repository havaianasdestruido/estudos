"""Monta um notebook .ipynb a partir de um arquivo .py com células.

Por que isso existe: editar JSON de notebook à mão é chato e fácil de errar
(vírgula, escape, saída gigante). Aqui você escreve Python normal, marca onde
começa cada célula e o script gera o .ipynb com os metadados do repositório.

Formato do arquivo de entrada
------------------------------

    # %% md
    # Título do notebook
    Texto em markdown, com $E=mc^2$ se quiser.

    # %% code
    print("olá")

    # %% code title="Preparar o ambiente"
    import rdkit

Regras:
  * `# %% md`   -> célula markdown (o conteúdo vai literalmente como está)
  * `# %% code` -> célula de código; `title="..."` vira o comentário `# @title`
  * linhas antes do primeiro marcador são ignoradas (use como nota/cabeçalho)
  * nada é executado: o notebook é gerado "limpo", sem saídas

Uso:
    python criar_notebook.py entrada.py saida.ipynb [--titulo "Meu notebook"]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

MARCADOR = re.compile(r"^#\s*%%\s*(md|code)\s*(.*)$", re.IGNORECASE)
TITLE = re.compile(r'title\s*=\s*"([^"]*)"')

METADADOS_BASE = {
    "authors": [{"name": "Material didático gerado a partir do script de células"}],
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py",
        "mimetype": "text/x-python",
        "name": "python",
        "nbconvert_exporter": "python",
        "pygments_lexer": "ipython3",
        "version": "3.11.2",
    },
}


def _id(indice: int) -> str:
    """id da célula (nbformat 4.5 exige; só precisa ser único)."""
    return f"celula-{indice + 1:02d}"


def ler_celulas(texto: str) -> list[dict]:
    """Divide o texto nos marcadores `# %% md` / `# %% code`."""
    celulas: list[dict] = []
    atual_tipo: str | None = None
    atual_linhas: list[str] = []
    atual_opcoes = ""

    def fechar():
        if atual_tipo is None:
            return
        corpo = "\n".join(atual_linhas).strip("\n")
        if atual_tipo == "md":
            if corpo.strip():
                celulas.append({"cell_type": "markdown", "id": _id(len(celulas)),
                                "metadata": {},
                                "source": corpo.splitlines(keepends=True)})
        else:
            if corpo.strip():
                titulo = TITLE.search(atual_opcoes)
                if titulo:
                    corpo = f'# @title {titulo.group(1)}\n' + corpo
                celulas.append({"cell_type": "code", "execution_count": None,
                                "id": _id(len(celulas)), "metadata": {}, "outputs": [],
                                "source": corpo.splitlines(keepends=True)})

    for linha in texto.splitlines():
        achado = MARCADOR.match(linha)
        if achado:
            fechar()
            atual_tipo = achado.group(1).lower()
            atual_opcoes = achado.group(2)
            atual_linhas = []
        elif atual_tipo is not None:
            atual_linhas.append(linha)
    fechar()
    return celulas


def montar_notebook(texto: str, titulo: str | None = None) -> dict:
    celulas = ler_celulas(texto)
    if not celulas:
        raise ValueError("nenhuma célula encontrada: use os marcadores '# %% md' e '# %% code'")
    metadados = json.loads(json.dumps(METADADOS_BASE))  # cópia
    if titulo:
        metadados["title"] = titulo
    return {"cells": celulas, "metadata": metadados, "nbformat": 4, "nbformat_minor": 5}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("entrada", type=Path, help="arquivo .py com as células")
    parser.add_argument("saida", type=Path, help="caminho do .ipynb a gerar")
    parser.add_argument("--titulo", default=None, help="metadado 'title' do notebook")
    args = parser.parse_args(argv)

    texto = args.entrada.read_text(encoding="utf-8")
    notebook = montar_notebook(texto, args.titulo or args.entrada.stem.replace("_", " ").title())

    args.saida.parent.mkdir(parents=True, exist_ok=True)
    args.saida.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
                          encoding="utf-8")
    print(f"{len(notebook['cells'])} células escritas em {args.saida}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
