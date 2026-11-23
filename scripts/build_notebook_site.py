#!/usr/bin/env python3
"""Render every repository Jupyter notebook into a static, browsable site."""

from __future__ import annotations

import html
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

from nbconvert import HTMLExporter


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "_site"
SUBJECT_LABELS = {
    "fisica": "Física",
    "matematica": "Matemática",
    "quimica": "Química",
    "musica": "Música",
    "outros": "Outros",
}


INDEX_CSS = """
:root {
  color-scheme: light;
  --ink: #172b3a;
  --muted: #5c6b78;
  --line: #dce5eb;
  --accent: #087e8b;
  --paper: #ffffff;
  --background: #f3f7f9;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--background);
  color: var(--ink);
  font: 16px/1.55 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}
header {
  padding: 3.5rem 1.25rem 3rem;
  background: #123348;
  color: white;
}
.header-inner, main, footer { width: min(100% - 2rem, 1080px); margin-inline: auto; }
header h1 { margin: 0 0 .5rem; font-size: clamp(2rem, 5vw, 3rem); line-height: 1.1; }
header p { max-width: 680px; margin: .75rem 0 0; color: #d5e8ef; }
.count { display: inline-block; margin-top: 1rem; color: #a7e2df; font-size: .9rem; font-weight: 650; }
main { padding-block: 2rem 3rem; }
.subject {
  margin: 0 0 1.5rem;
  padding: 1.5rem;
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 14px;
  box-shadow: 0 4px 14px rgb(23 43 58 / 4%);
}
.subject h2 { margin: 0 0 1rem; font-size: 1.5rem; }
.folder { margin-top: 1.25rem; }
.folder h3 { margin: 0 0 .4rem; color: var(--muted); font-size: .95rem; font-weight: 650; }
ul { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 310px), 1fr)); gap: .5rem; margin: 0; padding: 0; list-style: none; }
li a {
  display: block;
  height: 100%;
  padding: .75rem .9rem;
  border: 1px solid var(--line);
  border-radius: 9px;
  color: var(--ink);
  text-decoration: none;
  transition: border-color .15s ease, background .15s ease;
}
li a:hover, li a:focus-visible { border-color: var(--accent); background: #f1fbfb; }
.note { display: block; margin-top: .2rem; color: var(--muted); font-size: .8rem; overflow-wrap: anywhere; }
footer { padding: 0 0 2rem; color: var(--muted); font-size: .85rem; }
@media (max-width: 520px) { header { padding-block: 2.5rem; } .subject { padding: 1rem; } }
"""


def read_notebook_title(path: Path) -> str:
    """Use notebook metadata or its first Markdown heading as the link title."""
    notebook = json.loads(path.read_text(encoding="utf-8"))
    title = notebook.get("metadata", {}).get("title")
    if isinstance(title, str) and title.strip():
        return title.strip()

    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "markdown":
            continue
        source = cell.get("source", "")
        if isinstance(source, list):
            source = "".join(source)
        for line in source.splitlines():
            heading = re.match(r"^\s{0,3}#\s+(.+?)\s*#*\s*$", line)
            if heading:
                return heading.group(1).strip()

    return path.stem.replace("_", " ").replace("-", " ").capitalize()


def find_notebooks() -> list[Path]:
    """Find notebooks while ignoring hidden folders and the generated site."""
    notebooks = []
    for path in ROOT.rglob("*.ipynb"):
        relative = path.relative_to(ROOT)
        if any(part.startswith(".") or part == "_site" for part in relative.parts[:-1]):
            continue
        notebooks.append(path)
    return sorted(notebooks, key=lambda path: path.relative_to(ROOT).as_posix().casefold())


def folder_label(parts: tuple[str, ...]) -> str:
    labels = []
    for part in parts:
        words = part.replace("_", " ").replace("-", " ")
        labels.append(words if part.isupper() else words[:1].upper() + words[1:])
    return " / ".join(labels)


def link_to_rendered_notebooks(rendered_html: str) -> str:
    """Point relative notebook links at their generated HTML pages."""
    pattern = re.compile(r"""(href=["'])([^"']+)\.ipynb(["'])""", re.IGNORECASE)

    def replace(match: re.Match[str]) -> str:
        url = match.group(2)
        if url.startswith(("/", "#")) or re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", url):
            return match.group(0)
        return f"{match.group(1)}{url}.html{match.group(3)}"

    return pattern.sub(replace, rendered_html)


def build_index(notebooks: list[Path]) -> str:
    grouped: dict[str, dict[tuple[str, ...], list[tuple[str, Path]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for notebook in notebooks:
        relative = notebook.relative_to(ROOT)
        subject = relative.parts[0] if len(relative.parts) > 1 else "outros"
        folders = relative.parent.parts[1:] if len(relative.parts) > 1 else ()
        grouped[subject][folders].append((read_notebook_title(notebook), notebook))

    sections = []
    for subject in sorted(grouped, key=lambda name: SUBJECT_LABELS.get(name, name).casefold()):
        folder_groups = grouped[subject]
        subject_name = SUBJECT_LABELS.get(subject, folder_label((subject,)))
        content = [
            '<section class="subject">',
            f"  <h2>{html.escape(subject_name)}</h2>",
        ]
        for folders in sorted(folder_groups, key=lambda parts: "/".join(parts).casefold()):
            items = folder_groups[folders]
            content.append('  <div class="folder">')
            if folders:
                content.append(f"    <h3>{html.escape(folder_label(folders))}</h3>")
            content.append("    <ul>")
            for title, notebook in items:
                relative_html = notebook.relative_to(ROOT).with_suffix(".html")
                href = quote(relative_html.as_posix(), safe="/")
                location = html.escape(notebook.parent.relative_to(ROOT).as_posix())
                content.extend(
                    [
                        "      <li>",
                        f'        <a href="{href}">{html.escape(title)}'
                        f'<span class="note">{location}</span></a>',
                        "      </li>",
                    ]
                )
            content.extend(["    </ul>", "  </div>"])
        content.extend(["</section>"])
        sections.append("\n".join(content))

    total = len(notebooks)
    noun = "caderno" if total == 1 else "cadernos"
    return f"""<!doctype html>
<html lang="pt-BR">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Cadernos de estudo em HTML, gerados a partir dos notebooks do repositório.">
  <title>Estudos — cadernos</title>
  <style>{INDEX_CSS}</style>
</head>
<body>
  <header>
    <div class="header-inner">
      <h1>Cadernos de estudo</h1>
      <p>Material didático em texto e código, pronto para ler no navegador. Escolha um caderno para abrir a versão estática.</p>
      <span class="count">{total} {noun}</span>
    </div>
  </header>
  <main>
    {''.join(sections)}
  </main>
  <footer>Gerado automaticamente dos arquivos <code>.ipynb</code> pelo GitHub Actions.</footer>
</body>
</html>
"""


def main() -> None:
    notebooks = find_notebooks()
    if not notebooks:
        raise SystemExit("Nenhum arquivo .ipynb encontrado.")

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    exporter = HTMLExporter(
        template_name="lab",
        exclude_input_prompt=True,
        exclude_output_prompt=True,
    )

    for notebook in notebooks:
        relative_html = notebook.relative_to(ROOT).with_suffix(".html")
        destination = OUTPUT_DIR / relative_html
        destination.parent.mkdir(parents=True, exist_ok=True)
        rendered_html, _resources = exporter.from_filename(str(notebook))
        rendered_html = link_to_rendered_notebooks(rendered_html)
        destination.write_text(rendered_html, encoding="utf-8")

    (OUTPUT_DIR / "index.html").write_text(build_index(notebooks), encoding="utf-8")
    (OUTPUT_DIR / ".nojekyll").write_text("", encoding="utf-8")
    print(f"Rendered {len(notebooks)} notebooks into {OUTPUT_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
