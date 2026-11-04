"""Valida um notebook de química: estrutura, células e execução de cima a baixo.

Duas etapas (a segunda é opcional):

1. **Estrutura** (sempre roda, sem dependências):
   - JSON válido, `nbformat` 4, kernel Python;
   - pelo menos um título `# `, um objetivo e, em notebooks de química, uma
     visualização 3D (py3Dmol) e uma estrutura 2D (RDKit);
   - `EmbedMolecule` sempre acompanhado de `randomSeed` fixo (reprodutibilidade);
   - avisa sobre saídas grandes salvas (SVG/base64 gigante infla o .ipynb).

2. **Execução** (`--executar`): abre um kernel temporário e roda todas as células,
   reportando a célula que falhou. Precisa de `nbclient` + `ipykernel`
   (`pip install nbclient ipykernel`).

Uso:
    python validar_notebook.py notebook.ipynb
    python validar_notebook.py notebook.ipynb --executar
    python validar_notebook.py notebook.ipynb --executar --permitir-sem-3d
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

LIMITE_SAIDA_KB = 500


def texto_celula(celula: dict) -> str:
    fonte = celula.get("source", "")
    return "".join(fonte) if isinstance(fonte, list) else fonte


def validar_estrutura(caminho: Path, exigir_3d: bool = True) -> list[str]:
    problemas: list[str] = []
    try:
        notebook = json.loads(caminho.read_text(encoding="utf-8"))
    except json.JSONDecodeError as erro:
        return [f"JSON inválido: {erro}"]

    if notebook.get("nbformat") != 4:
        problemas.append("nbformat deve ser 4")
    if notebook.get("metadata", {}).get("kernelspec", {}).get("language") != "python":
        problemas.append("metadata.kernelspec.language deve ser 'python'")

    celulas = notebook.get("cells", [])
    if not celulas:
        return problemas + ["notebook sem células"]

    markdown = [texto_celula(c) for c in celulas if c.get("cell_type") == "markdown"]
    codigo = [texto_celula(c) for c in celulas if c.get("cell_type") == "code"]
    tudo_md = "\n".join(markdown)
    tudo_py = "\n".join(codigo)

    if not any(linha.lstrip().startswith("# ") for linha in tudo_md.splitlines()):
        problemas.append("falta o título '#' na primeira célula markdown")
    if "Objetivo" not in tudo_md:
        problemas.append("falta a linha '**Objetivo:**' (padrão do repositório)")
    if "Como usar" not in tudo_md:
        problemas.append("falta a seção '### Como usar este notebook'")

    # As exigências de química só valem quando o notebook é de química:
    # assim o mesmo validador serve para os notebooks de física/matemática.
    eh_quimica = any(pista in tudo_py + tudo_md
                     for pista in ("rdkit", "py3Dmol", "SMILES", "molécula", "molecula"))
    if eh_quimica:
        if exigir_3d and "py3Dmol" not in tudo_py:
            problemas.append("não há visualização 3D interativa (py3Dmol)")
        if "rdMolDraw2D" not in tudo_py and "MolsToGridImage" not in tudo_py:
            problemas.append("não há estrutura 2D com risquinhos (rdkit Draw)")
    if "matplotlib" not in tudo_py:
        problemas.append("não há gráfico 2D de dados (matplotlib) — recomendado")

    for i, codigo_celula in enumerate(codigo):
        if "EmbedMolecule" in codigo_celula and "randomSeed" not in codigo_celula:
            problemas.append(f"célula de código {i}: EmbedMolecule sem randomSeed fixo")
        for achado in re.finditer(r"MolFromSmiles\(([^)]*)\)", codigo_celula):
            argumento = achado.group(1).strip()
            if argumento.startswith(("f'", 'f"', "'", '"')) and "%" in argumento:
                problemas.append(f"célula de código {i}: SMILES montado com % — confira a validade")

    tamanho_saidas = sum(len(json.dumps(c.get("outputs", ""), ensure_ascii=False))
                         for c in celulas if c.get("cell_type") == "code")
    if tamanho_saidas > LIMITE_SAIDA_KB * 1024:
        problemas.append(f"saídas salvas somam {tamanho_saidas / 1024:.0f} KB "
                         f"(> {LIMITE_SAIDA_KB} KB): limpe as saídas 3D/grandes "
                         "com 'Cell > All Output > Clear'")
    return problemas


def executar(caminho: Path, tempo_limite: int = 600) -> list[str]:
    try:
        import nbformat
        from nbclient import NotebookClient
    except ImportError:
        return ["para executar, instale: pip install nbclient ipykernel"]

    notebook = nbformat.read(caminho, as_version=4)
    cliente = NotebookClient(notebook, timeout=tempo_limite,
                             kernel_name="python3", resources={"metadata": {"path": str(caminho.parent)}})
    try:
        cliente.execute()
    except Exception as erro:  # noqa: BLE001 - queremos a mensagem crua do kernel
        mensagem = str(erro)
        numero = re.search(r"Cell Execution Error.*", mensagem, re.S)
        return [f"falha ao executar: {mensagem[:1200]}"]
    return []


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("notebook", type=Path)
    parser.add_argument("--executar", action="store_true",
                        help="roda o notebook de ponta a ponta num kernel temporário")
    parser.add_argument("--permitir-sem-3d", action="store_true",
                        help="não exigir visualização 3D interativa")
    parser.add_argument("--mostrar", action="store_true",
                        help="imprime o relatório mesmo quando está tudo certo")
    args = parser.parse_args(argv)

    if not args.notebook.exists():
        print(f"arquivo não encontrado: {args.notebook}")
        return 2

    problemas = validar_estrutura(args.notebook, exigir_3d=not args.permitir_sem_3d)
    if args.executar:
        problemas += executar(args.notebook)

    if problemas:
        print(f"❌ {args.notebook}: {len(problemas)} problema(s)")
        for problema in problemas:
            print(f"  - {problema}")
        return 1

    print(f"✅ {args.notebook}: ok" + (" (executado)" if args.executar else
                                       " (só estrutura; use --executar para rodar)"))
    if args.mostrar:
        notebook = json.loads(args.notebook.read_text(encoding="utf-8"))
        codigo = [c for c in notebook.get("cells", []) if c.get("cell_type") == "code"]
        markdown = [c for c in notebook.get("cells", []) if c.get("cell_type") == "markdown"]
        print(f"  {len(markdown)} células de markdown, {len(codigo)} de código")
    return 0


if __name__ == "__main__":
    sys.exit(main())
