from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"


def should_include_png(path: Path) -> bool:
    parts_lower = tuple(part.lower() for part in path.relative_to(PLOTS_ROOT).parts)
    if "__pycache__" in parts_lower:
        return False
    if any(part.startswith("_tmp") for part in parts_lower):
        return False
    return True


def sort_key(path: Path) -> tuple[str, ...]:
    return tuple(part.lower() for part in path.parts)


def collect_pngs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted(
        [path for path in root.rglob("*.png") if should_include_png(path)],
        key=lambda path: sort_key(path.relative_to(root)),
    )


def direct_pngs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted(
        [path for path in root.glob("*.png") if should_include_png(path)],
        key=lambda path: sort_key(path.relative_to(root)),
    )


def iter_gallery_dirs() -> list[Path]:
    dirs = {PLOTS_ROOT}
    for png in collect_pngs(PLOTS_ROOT):
        dirs.update(png.parents)
    return sorted(
        [path for path in dirs if path == PLOTS_ROOT or PLOTS_ROOT in path.parents],
        key=lambda path: (len(path.relative_to(PLOTS_ROOT).parts), path.as_posix().lower()),
    )


def group_name(rel_path: Path) -> str:
    return "Plots"


def page_title(gallery_dir: Path) -> str:
    if gallery_dir == PLOTS_ROOT:
        return "ePC-SAFT Plot Gallery"
    rel = gallery_dir.relative_to(PLOTS_ROOT)
    return " / ".join(part.replace("_", " ").title() for part in rel.parts)


def back_href(gallery_dir: Path) -> str | None:
    if gallery_dir == PLOTS_ROOT:
        return None
    return "../index.html"


def child_gallery_links(gallery_dir: Path) -> list[tuple[str, str, int]]:
    children = []
    for child in sorted([p for p in gallery_dir.iterdir() if p.is_dir()], key=lambda p: p.name.lower()):
        count = len(collect_pngs(child))
        if count:
            children.append((child.name.replace("_", " ").title(), f"{child.name}/index.html", count))
    return children


def _folder_options(gallery_dir: Path) -> str:
    options = []
    if gallery_dir != PLOTS_ROOT:
        options.append('        <option value="../index.html">..</option>')
    for name, href, count in child_gallery_links(gallery_dir):
        options.append(f'        <option value="{html.escape(href)}">{html.escape(name)} ({count})</option>')
    if not options:
        return '        <option value="">No subfolders</option>'
    return "\n".join(options)


def render_gallery_page(gallery_dir: Path, pngs: Iterable[Path]) -> str:
    rel_pngs = [path.relative_to(gallery_dir) for path in pngs]
    child_links = child_gallery_links(gallery_dir)

    sections: list[tuple[str, list[Path]]] = []
    current_name = None
    current_paths: list[Path] = []
    for rel_path in rel_pngs:
        name = group_name(rel_path)
        if name != current_name:
            if current_paths:
                sections.append((current_name or "Top Level", current_paths))
            current_name = name
            current_paths = [rel_path]
        else:
            current_paths.append(rel_path)
    if current_paths:
        sections.append((current_name or "Top Level", current_paths))

    toc_items = []
    for index, (name, paths) in enumerate(sections, start=1):
        anchor = f"section-{index}"
        toc_items.append(
            f'        <button type="button" class="section-button" data-target="{html.escape(anchor)}">'
            f"{html.escape(name)} ({len(paths)})</button>"
        )

    child_items = [
        f'        <a class="child-card" href="{html.escape(href)}">'
        f"<span>{html.escape(name)}</span><strong>{count} PNGs</strong></a>"
        for name, href, count in child_links
    ]

    section_blocks = []
    for index, (name, paths) in enumerate(sections, start=1):
        anchor = f"section-{index}"
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
            f"    <section id=\"{anchor}\" class=\"gallery-section\">\n"
            f"      <h2>{html.escape(name)}</h2>\n"
            + "\n".join(cards)
            + "\n    </section>"
        )

    title = page_title(gallery_dir)
    rel_label = "." if gallery_dir == PLOTS_ROOT else gallery_dir.relative_to(PLOTS_ROOT).as_posix()
    summary = f"Browsable gallery for {len(rel_pngs)} PNG files under docs/plots/{rel_label}."
    default_section = ""
    back = back_href(gallery_dir)
    back_link = f'      <a class="back-link" href="{back}">Back</a>\n' if back else ""
    empty_message = ""
    if not section_blocks:
        empty_message = (
            '    <section class="empty-state">\n'
            "      <h2>No plots in this folder</h2>\n"
            "      <p>No plots are shown until you open a folder that contains PNG files.</p>\n"
            "    </section>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --page-max: 1200px;
      --nav-bg: rgba(255, 255, 255, 0.96);
      --border: #d9d9d9;
      --text: #1b1b1b;
      --muted: #5f5f5f;
      --accent: #0f5fbf;
      --surface: #f7f8fa;
    }}

    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
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
    h1 {{ margin: 0 0 0.35rem; font-size: 1.5rem; font-weight: 700; }}
    .summary {{ margin: 0 0 0.8rem; color: var(--muted); font-size: 0.95rem; }}
    .folder-nav, nav, .child-grid {{ display: flex; flex-wrap: wrap; gap: 0.5rem; align-items: center; }}
    .folder-nav {{ margin-bottom: 0.75rem; }}
    .folder-nav label {{ font-weight: 600; font-size: 0.92rem; }}
    .folder-select {{
      min-width: min(28rem, 100%);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 0.42rem 0.55rem;
      background: #fff;
      color: var(--text);
      font: inherit;
    }}
    .section-button, .back-link, .child-card {{
      appearance: none;
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 0.35rem 0.8rem;
      font-size: 0.92rem;
      background: #fff;
      color: var(--accent);
      cursor: pointer;
      text-decoration: none;
      transition: background 120ms ease, color 120ms ease, border-color 120ms ease;
    }}
    .section-button.is-active {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
    .child-card {{ display: inline-flex; gap: 0.5rem; align-items: center; }}
    .child-card strong {{ color: var(--muted); font-weight: 500; }}
    .back-link {{ display: inline-flex; margin-bottom: 0.8rem; }}
    main {{ max-width: var(--page-max); margin: 0 auto; padding: 1.25rem; }}
    .child-section {{ margin: 0 0 1.25rem; }}
    .child-section h2, .gallery-section > h2 {{
      margin: 0 0 0.85rem;
      font-size: 1.35rem;
      border-bottom: 2px solid #efefef;
      padding-bottom: 0.4rem;
    }}
    .gallery-section {{ display: none; margin: 0 0 2.5rem; padding-top: 0.25rem; }}
    .gallery-section.is-active {{ display: block; }}
    .empty-state {{
      max-width: 46rem;
      padding: 1rem;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: var(--surface);
    }}
    .empty-state p {{ margin: 0; color: var(--muted); }}
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
    .image-link {{ display: block; text-decoration: none; }}
    img {{ display: block; width: 100%; max-width: 100%; height: auto; border-radius: 8px; background: #fff; }}
  </style>
</head>
<body>
  <header>
    <div class="header-inner">
{back_link}\
      <h1>{html.escape(title)}</h1>
      <p class="summary">{html.escape(summary)}</p>
      <div class="folder-nav">
        <label for="folder-select">Subfolder</label>
        <select id="folder-select" class="folder-select">
          <option value="">Choose a folder</option>
{_folder_options(gallery_dir)}
        </select>
      </div>
      <nav>
{chr(10).join(toc_items)}
      </nav>
    </div>
  </header>
  <main>
    <section class="child-section">
      <h2>Subfolders</h2>
      <div class="child-grid">
{chr(10).join(child_items) if child_items else '        <span class="summary">No child galleries.</span>'}
      </div>
    </section>
{chr(10).join(section_blocks)}
{empty_message}
  </main>
  <script>
    (function() {{
      const folderSelect = document.getElementById('folder-select');
      const buttons = Array.from(document.querySelectorAll('.section-button'));
      const sections = Array.from(document.querySelectorAll('.gallery-section'));
      const defaultSection = {json.dumps(default_section)};

      function activateSection(sectionId) {{
        const targetId = sectionId || '';
        sections.forEach((section) => section.classList.toggle('is-active', section.id === targetId));
        buttons.forEach((button) => button.classList.toggle('is-active', button.dataset.target === targetId));
        if (targetId) {{
          if (window.location.hash !== '#' + targetId) {{
            history.replaceState(null, '', '#' + targetId);
          }}
          window.scrollTo({{ top: 0, behavior: 'auto' }});
        }}
      }}

      buttons.forEach((button) => button.addEventListener('click', () => activateSection(button.dataset.target)));
      if (folderSelect) {{
        folderSelect.addEventListener('change', () => {{
          if (folderSelect.value) {{
            window.location.href = folderSelect.value;
          }}
        }});
      }}
      const initialHash = window.location.hash ? window.location.hash.slice(1) : '';
      const hasInitialTarget = sections.some((section) => section.id === initialHash);
      activateSection(hasInitialTarget ? initialHash : '');
    }})();
  </script>
</body>
</html>
"""


def remove_stale_indexes(valid_dirs: set[Path]) -> None:
    for index_path in PLOTS_ROOT.rglob("index.html"):
        if index_path.parent not in valid_dirs:
            index_path.unlink()


def main() -> None:
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    gallery_dirs = iter_gallery_dirs()
    valid_dirs = set(gallery_dirs)
    remove_stale_indexes(valid_dirs)
    for gallery_dir in gallery_dirs:
        pngs = direct_pngs(gallery_dir)
        if not pngs and not child_gallery_links(gallery_dir) and gallery_dir != PLOTS_ROOT:
            continue
        (gallery_dir / "index.html").write_text(render_gallery_page(gallery_dir, pngs), encoding="utf-8")


if __name__ == "__main__":
    main()
