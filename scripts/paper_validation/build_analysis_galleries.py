from __future__ import annotations

import html
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent


def sort_key(path: Path) -> tuple:
    parts = [part.lower() for part in path.parts]
    return tuple(parts)


def iter_analysis_dirs() -> list[Path]:
    return sorted(
        [path for path in ROOT.iterdir() if path.is_dir() and not path.name.startswith("__")],
        key=lambda path: path.name.lower(),
    )


def collect_pngs(analysis_dir: Path) -> list[Path]:
    return sorted(
        [
            path
            for path in analysis_dir.rglob("*.png")
            if "__pycache__" not in path.parts
            and not any(part.startswith("_tmp") for part in path.parts)
            and path.name.lower() != "index.html"
        ],
        key=lambda path: sort_key(path.relative_to(analysis_dir)),
    )


def group_name(rel_path: Path) -> str:
    if len(rel_path.parts) <= 1:
        return "Top Level"
    first = rel_path.parts[0]
    return first.replace("_", " ").title()


def render_analysis_page(analysis_dir: Path, pngs: Iterable[Path]) -> str:
    rel_pngs = [path.relative_to(analysis_dir) for path in pngs]
    sections: list[tuple[str, list[Path]]] = []
    current_name = None
    current_paths: list[Path] = []
    for rel_path in rel_pngs:
        name = group_name(rel_path)
        if name != current_name:
            if current_paths:
                sections.append((current_name, current_paths))
            current_name = name
            current_paths = [rel_path]
        else:
            current_paths.append(rel_path)
    if current_paths:
        sections.append((current_name, current_paths))

    toc_items = []
    section_blocks = []
    for index, (name, paths) in enumerate(sections, start=1):
        anchor = f"section-{index}"
        toc_items.append(f'        <a href="#{anchor}">{html.escape(name)}</a>')
        cards = []
        for rel_path in paths:
            rel_text = rel_path.as_posix()
            alt = rel_path.stem.replace("_", " ")
            cards.append(
                "      <article class=\"image-card\">"
                f"<p class=\"caption\">{html.escape(rel_text)}</p>"
                f"<a class=\"image-link\" href=\"{html.escape(rel_text)}\">"
                f"<img src=\"{html.escape(rel_text)}\" alt=\"{html.escape(alt)}\"></a>"
                "</article>"
            )
        section_blocks.append(
            f"    <section id=\"{anchor}\">\n"
            f"      <h2>{html.escape(name)}</h2>\n"
            + "\n".join(cards)
            + "\n    </section>"
        )

    summary = f"Scrollable local gallery for {len(rel_pngs)} PNG outputs under {analysis_dir.name}."
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(analysis_dir.name.replace('_', ' '))} Gallery</title>
  <style>
    :root {{
      color-scheme: light;
      --page-max: 1200px;
      --nav-bg: rgba(255, 255, 255, 0.96);
      --border: #d9d9d9;
      --text: #1b1b1b;
      --muted: #5f5f5f;
      --accent: #0f5fbf;
    }}

    * {{
      box-sizing: border-box;
    }}

    html {{
      scroll-behavior: smooth;
    }}

    body {{
      margin: 0;
      font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
      color: var(--text);
      background: #ffffff;
      line-height: 1.45;
    }}

    header {{
      position: sticky;
      top: 0;
      z-index: 10;
      background: var(--nav-bg);
      border-bottom: 1px solid var(--border);
      backdrop-filter: blur(8px);
    }}

    .header-inner {{
      max-width: var(--page-max);
      margin: 0 auto;
      padding: 1rem 1.25rem 0.85rem;
    }}

    h1 {{
      margin: 0 0 0.35rem;
      font-size: 1.5rem;
      font-weight: 700;
    }}

    .summary {{
      margin: 0 0 0.8rem;
      color: var(--muted);
      font-size: 0.95rem;
    }}

    nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
    }}

    nav a {{
      text-decoration: none;
      color: var(--accent);
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 0.35rem 0.8rem;
      font-size: 0.92rem;
      background: #fff;
    }}

    .back-link {{
      display: inline-flex;
      align-items: center;
      text-decoration: none;
      color: var(--accent);
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 0.35rem 0.8rem;
      font-size: 0.92rem;
      background: #fff;
      margin-bottom: 0.8rem;
    }}

    main {{
      max-width: var(--page-max);
      margin: 0 auto;
      padding: 1.25rem;
    }}

    section {{
      margin: 0 0 2.5rem;
      padding-top: 0.25rem;
    }}

    section > h2 {{
      margin: 0 0 0.85rem;
      font-size: 1.35rem;
      border-bottom: 2px solid #efefef;
      padding-bottom: 0.4rem;
    }}

    .image-card {{
      margin: 0 0 1.5rem;
      padding: 0.85rem;
      border: 1px solid var(--border);
      border-radius: 12px;
      background: #fcfcfc;
    }}

    .caption {{
      margin: 0 0 0.6rem;
      font-size: 0.92rem;
      color: var(--muted);
      font-family: Consolas, "Courier New", monospace;
      word-break: break-word;
    }}

    .image-link {{
      display: block;
      text-decoration: none;
    }}

    img {{
      display: block;
      width: 100%;
      max-width: 100%;
      height: auto;
      border-radius: 8px;
      background: #fff;
    }}
  </style>
</head>
<body>
  <header>
    <div class="header-inner">
      <a class="back-link" href="../index.html">Back to Paper Validation Home</a>
      <h1>{html.escape(analysis_dir.name.replace('_', ' '))} Gallery</h1>
      <p class="summary">{html.escape(summary)}</p>
      <nav>
{chr(10).join(toc_items)}
      </nav>
    </div>
  </header>
  <main>
{chr(10).join(section_blocks)}
  </main>
</body>
</html>
"""


def render_root_page(analysis_dirs: Iterable[Path]) -> str:
    items = []
    for analysis_dir in analysis_dirs:
        png_count = len(collect_pngs(analysis_dir))
        items.append(
            "      <li>"
            f"<a href=\"{html.escape(analysis_dir.name)}/index.html\">{html.escape(analysis_dir.name.replace('_', ' '))}</a>"
            f"<span>{png_count} PNGs</span>"
            "</li>"
        )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Paper Validation Galleries</title>
  <style>
    body {{
      margin: 0;
      font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
      background: #ffffff;
      color: #1b1b1b;
    }}

    main {{
      max-width: 900px;
      margin: 0 auto;
      padding: 2rem 1.25rem 3rem;
    }}

    h1 {{
      margin: 0 0 0.5rem;
      font-size: 1.8rem;
    }}

    p {{
      margin: 0 0 1.5rem;
      color: #5f5f5f;
    }}

    ul {{
      list-style: none;
      padding: 0;
      margin: 0;
      display: grid;
      gap: 0.8rem;
    }}

    li {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 1rem;
      border: 1px solid #d9d9d9;
      border-radius: 12px;
      padding: 0.85rem 1rem;
      background: #fcfcfc;
    }}

    a {{
      color: #0f5fbf;
      text-decoration: none;
      font-weight: 600;
    }}

    span {{
      color: #5f5f5f;
      font-size: 0.92rem;
      white-space: nowrap;
    }}
  </style>
</head>
<body>
  <main>
    <h1>Paper Validation Galleries</h1>
    <p>Static local galleries for all paper-validation analysis folders.</p>
    <ul>
{chr(10).join(items)}
    </ul>
  </main>
</body>
</html>
"""


def main() -> None:
    analysis_dirs = iter_analysis_dirs()
    for analysis_dir in analysis_dirs:
        pngs = collect_pngs(analysis_dir)
        if not pngs:
            continue
        (analysis_dir / "index.html").write_text(
            render_analysis_page(analysis_dir, pngs),
            encoding="utf-8",
        )
    (ROOT / "index.html").write_text(render_root_page(analysis_dirs), encoding="utf-8")


if __name__ == "__main__":
    main()
