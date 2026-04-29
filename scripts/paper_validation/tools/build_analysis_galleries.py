from __future__ import annotations

import html
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
CSV_BACKFILL_MARKER = 'epcsaft-interactive-source="csv_backfill"'


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


def page_title() -> str:
    return "ePC-SAFT Plot Gallery"


def _source_path_for_output(rel_path: str) -> str:
    parts = rel_path.split("/")
    if not parts:
        return rel_path
    if parts[0] == "paper_validation" and len(parts) >= 2:
        source_study = parts[1] if parts[1].endswith("_analysis") else f"{parts[1]}_analysis"
        return "/".join(["scripts", "paper_validation", source_study, *parts[2:]])
    if parts[0] == "fits":
        return "/".join(["scripts", *parts])
    if parts[0] == "tests":
        return rel_path
    return "/".join(["docs", "plots", *parts])


def _folder_for(path: str) -> str:
    folder = str(Path(path).parent).replace("\\", "/")
    return "" if folder == "." else folder


def _interactive_source(html_path: Path) -> str:
    if not html_path.exists():
        return "static_only"
    try:
        head = html_path.read_text(encoding="utf-8", errors="ignore")[:2000]
    except OSError:
        return "native_plotly"
    is_csv_backfill = CSV_BACKFILL_MARKER in head or (
        "epcsaft-interactive-source" in head and 'content="csv_backfill"' in head
    )
    return "csv_backfill" if is_csv_backfill else "native_plotly"


def image_manifest(pngs: list[Path]) -> list[dict[str, str]]:
    manifest = []
    for png in pngs:
        output_path = png.relative_to(PLOTS_ROOT).as_posix()
        svg = png.with_suffix(".svg")
        interactive_html = png.with_suffix(".html")
        interactive_source = _interactive_source(interactive_html)
        data = png.parent / "data" / f"{png.stem}_plot_data.csv"
        source_path = _source_path_for_output(output_path)
        manifest.append(
            {
                "path": output_path,
                "folder": _folder_for(source_path),
                "output_path": output_path,
                "output_folder": _folder_for(output_path),
                "svg_path": svg.relative_to(PLOTS_ROOT).as_posix() if svg.exists() else "",
                "html_path": interactive_html.relative_to(PLOTS_ROOT).as_posix() if interactive_html.exists() else "",
                "interactive_source": interactive_source,
                "data_path": data.relative_to(PLOTS_ROOT).as_posix() if data.exists() else "",
                "source_path": source_path,
                "source_folder": _folder_for(source_path),
                "name": png.name,
                "title": png.stem.replace("_", " ").replace("-", " "),
            }
        )
    return manifest


_image_manifest = image_manifest


def _safe_json(data: object) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def render_gallery_page(root: Path, pngs: list[Path]) -> str:
    manifest = image_manifest(pngs)
    title = page_title()
    summary = f"{len(manifest)} PNG files under docs/plots. Select folders on the left to show every image in that subtree."

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f4f6f8;
      --panel: #ffffff;
      --panel-2: #f9fafb;
      --text: #17202a;
      --muted: #607086;
      --border: #d9e1ea;
      --border-strong: #b8c5d3;
      --blue: #195fb8;
      --blue-soft: #e8f1ff;
      --teal: #137d72;
      --shadow: 0 8px 24px rgba(25, 43, 65, 0.08);
      --sidebar: clamp(280px, 27vw, 420px);
      --sidebar-min: 220px;
      --sidebar-max: min(620px, 55vw);
      --sidebar-resizer: 8px;
      --radius: 8px;
    }}

    * {{ box-sizing: border-box; }}
    html, body {{ height: 100%; }}
    body {{
      margin: 0;
      overflow: hidden;
      font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
      color: var(--text);
      background: var(--bg);
      line-height: 1.4;
    }}
    button, input {{ font: inherit; }}
    body.sidebar-resizing {{
      cursor: col-resize;
      user-select: none;
    }}
    .app-shell {{
      display: grid;
      grid-template-columns: var(--sidebar) var(--sidebar-resizer) minmax(0, 1fr);
      height: 100vh;
    }}
    .sidebar {{
      display: flex;
      min-width: 0;
      flex-direction: column;
      background: var(--panel);
    }}
    .sidebar-resizer {{
      position: relative;
      width: var(--sidebar-resizer);
      border-inline: 1px solid transparent;
      background: linear-gradient(to right, transparent 0, transparent 3px, var(--border) 3px, var(--border) 4px, transparent 4px);
      cursor: col-resize;
      touch-action: none;
    }}
    .sidebar-resizer:hover,
    .sidebar-resizer:focus,
    body.sidebar-resizing .sidebar-resizer {{
      border-inline-color: rgba(25, 95, 184, 0.28);
      background: linear-gradient(to right, transparent 0, transparent 2px, var(--blue) 2px, var(--blue) 5px, transparent 5px);
      outline: none;
    }}
    .brand {{
      padding: 18px 18px 14px;
      border-bottom: 1px solid var(--border);
    }}
    h1 {{
      margin: 0;
      font-size: 1.18rem;
      line-height: 1.2;
      letter-spacing: 0;
    }}
    .summary {{
      margin: 7px 0 0;
      color: var(--muted);
      font-size: 0.88rem;
    }}
    .toolbar {{
      display: grid;
      gap: 10px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--border);
      background: var(--panel-2);
    }}
    .search {{
      width: 100%;
      min-height: 34px;
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 7px 10px;
      background: #fff;
      color: var(--text);
    }}
    .actions {{
      display: flex;
      gap: 8px;
      align-items: center;
    }}
    .mode-buttons {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
    }}
    .mode-button {{
      min-height: 32px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #fff;
      color: var(--muted);
      cursor: pointer;
      font-size: 0.82rem;
    }}
    .mode-button.is-active {{
      border-color: var(--blue);
      background: var(--blue-soft);
      color: var(--blue);
      font-weight: 650;
    }}
    .icon-button {{
      width: 34px;
      height: 34px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #fff;
      color: var(--text);
      cursor: pointer;
    }}
    .icon-button:hover {{ border-color: var(--blue); color: var(--blue); }}
    .tree-wrap {{
      min-height: 0;
      overflow: auto;
      padding: 10px 8px 18px;
    }}
    .tree {{
      margin: 0;
      padding: 0;
      list-style: none;
      font-size: 0.91rem;
    }}
    .tree ul {{
      margin: 0;
      padding: 0;
      list-style: none;
    }}
    .folder-row {{
      display: grid;
      grid-template-columns: 24px 22px minmax(0, 1fr) auto;
      align-items: center;
      min-height: 31px;
      gap: 4px;
      padding-right: 6px;
      border-radius: 6px;
    }}
    .folder-row:hover {{ background: #eef4fb; }}
    .folder-row.is-selected {{ background: var(--blue-soft); }}
    .toggle {{
      width: 22px;
      height: 22px;
      border: 0;
      border-radius: 4px;
      background: transparent;
      color: var(--muted);
      cursor: pointer;
    }}
    .toggle:disabled {{
      cursor: default;
      opacity: 0.2;
    }}
    .folder-check {{
      width: 15px;
      height: 15px;
      accent-color: var(--blue);
    }}
    .folder-name {{
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      cursor: pointer;
    }}
    .folder-count {{
      color: var(--muted);
      font-size: 0.78rem;
      font-variant-numeric: tabular-nums;
    }}
    .content {{
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      min-width: 0;
      height: 100vh;
    }}
    .content-header {{
      display: grid;
      gap: 12px;
      padding: 16px 20px;
      border-bottom: 1px solid var(--border);
      background: rgba(255, 255, 255, 0.9);
      backdrop-filter: blur(8px);
    }}
    .content-title-row {{
      display: flex;
      align-items: end;
      justify-content: space-between;
      gap: 16px;
    }}
    h2 {{
      margin: 0;
      font-size: 1.16rem;
      line-height: 1.2;
      letter-spacing: 0;
    }}
    .meta {{
      margin: 4px 0 0;
      color: var(--muted);
      font-size: 0.9rem;
    }}
    .view-tools {{
      display: flex;
      gap: 10px;
      align-items: center;
      color: var(--muted);
      font-size: 0.86rem;
      flex-wrap: wrap;
      justify-content: flex-end;
    }}
    .tile-control {{
      display: inline-flex;
      gap: 8px;
      align-items: center;
    }}
    .plot-mode-buttons {{
      display: inline-grid;
      grid-template-columns: 1fr 1fr;
      gap: 6px;
      min-width: 190px;
    }}
    .range {{
      width: 140px;
      accent-color: var(--teal);
    }}
    .chips {{
      display: flex;
      gap: 8px;
      min-height: 28px;
      overflow-x: auto;
      padding-bottom: 2px;
    }}
    .chip {{
      flex: 0 0 auto;
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 4px 10px;
      background: #fff;
      color: var(--muted);
      font-size: 0.82rem;
    }}
    .gallery-scroll {{
      min-height: 0;
      overflow: auto;
      padding: 20px;
    }}
    .gallery-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(var(--card-min, 280px), 1fr));
      gap: 16px;
      align-items: start;
    }}
    .image-card {{
      min-width: 0;
      overflow: hidden;
      border: 1px solid var(--border);
      border-radius: var(--radius);
      background: var(--panel);
      box-shadow: var(--shadow);
    }}
    .image-head {{
      display: grid;
      gap: 3px;
      padding: 10px 12px;
      border-bottom: 1px solid var(--border);
      background: var(--panel-2);
    }}
    .image-title {{
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      font-size: 0.91rem;
      font-weight: 650;
    }}
    .image-path {{
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      color: var(--muted);
      font-family: Consolas, "Courier New", monospace;
      font-size: 0.76rem;
    }}
    .image-actions {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      padding-top: 4px;
    }}
    .resource-link {{
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-height: 24px;
      border: 1px solid var(--border);
      border-radius: 5px;
      padding: 2px 7px;
      background: #fff;
      color: var(--blue);
      font-size: 0.74rem;
      font-weight: 650;
      text-decoration: none;
      cursor: pointer;
    }}
    .resource-link:hover {{
      border-color: var(--blue);
      background: var(--blue-soft);
    }}
    button.resource-link {{
      font: inherit;
    }}
    .asset-badge {{
      display: inline-flex;
      align-items: center;
      min-height: 24px;
      border: 1px solid #c9d5e2;
      border-radius: 5px;
      padding: 2px 7px;
      background: #eef3f8;
      color: var(--muted);
      font-size: 0.74rem;
      font-weight: 650;
    }}
    .plot-preview {{
      display: block;
      min-height: 220px;
      padding: 12px;
      background:
        linear-gradient(45deg, #f7f9fb 25%, transparent 25%),
        linear-gradient(-45deg, #f7f9fb 25%, transparent 25%),
        linear-gradient(45deg, transparent 75%, #f7f9fb 75%),
        linear-gradient(-45deg, transparent 75%, #f7f9fb 75%);
      background-color: #fff;
      background-position: 0 0, 0 8px, 8px -8px, -8px 0;
      background-size: 16px 16px;
    }}
    .plot-preview.is-interactive {{
      min-height: 430px;
      padding: 0;
      background: #fff;
    }}
    .interactive-frame {{
      display: block;
      width: 100%;
      height: clamp(430px, 54vw, 760px);
      min-height: 430px;
      border: 0;
      background: #fff;
    }}
    .interactive-placeholder {{
      display: grid;
      min-height: 430px;
      place-items: center;
      padding: 22px;
      background: #fff;
      color: var(--muted);
      text-align: center;
    }}
    img {{
      display: block;
      width: 100%;
      min-height: 180px;
      max-height: 620px;
      height: auto;
      object-fit: contain;
      border-radius: 6px;
      background: #fff;
    }}
    .empty-state {{
      max-width: 560px;
      border: 1px dashed var(--border-strong);
      border-radius: var(--radius);
      padding: 18px;
      background: #fff;
      color: var(--muted);
    }}
    .empty-state strong {{
      display: block;
      margin-bottom: 4px;
      color: var(--text);
    }}
    .asset-modal,
    .data-modal {{
      position: fixed;
      inset: 0;
      z-index: 20;
      display: grid;
      place-items: center;
      padding: 24px;
      background: rgba(23, 32, 42, 0.42);
    }}
    .asset-dialog,
    .data-dialog {{
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      width: min(1100px, 96vw);
      max-height: min(760px, 92vh);
      overflow: hidden;
      border: 1px solid var(--border-strong);
      border-radius: var(--radius);
      background: var(--panel);
      box-shadow: 0 22px 70px rgba(23, 32, 42, 0.28);
    }}
    .asset-dialog-head,
    .data-dialog-head {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 12px 14px;
      border-bottom: 1px solid var(--border);
      background: var(--panel-2);
    }}
    .asset-dialog-title,
    .data-dialog-title {{
      min-width: 0;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      font-weight: 650;
    }}
    .asset-dialog-body,
    .data-dialog-body {{
      min-height: 0;
      overflow: auto;
      padding: 12px 14px 16px;
    }}
    .asset-dialog-body {{
      display: grid;
      place-items: center;
      background: #fff;
    }}
    .asset-preview-image {{
      display: block;
      width: 100%;
      max-width: 100%;
      min-height: 0;
      max-height: 72vh;
      object-fit: contain;
      border-radius: 6px;
      background: #fff;
    }}
    .asset-preview-frame {{
      display: block;
      width: 100%;
      height: min(72vh, 760px);
      min-height: 420px;
      border: 0;
      background: #fff;
    }}
    .data-note {{
      margin: 0 0 10px;
      color: var(--muted);
      font-size: 0.84rem;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.78rem;
    }}
    .data-table th,
    .data-table td {{
      max-width: 260px;
      border: 1px solid var(--border);
      padding: 4px 6px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      text-align: left;
    }}
    .data-table th {{
      position: sticky;
      top: 0;
      background: var(--blue-soft);
      color: var(--text);
      z-index: 1;
    }}
    .hidden {{ display: none !important; }}

    @media (max-width: 860px) {{
      body {{ overflow: auto; }}
      .app-shell {{
        grid-template-columns: 1fr;
        height: auto;
        min-height: 100vh;
      }}
      .sidebar, .content {{ height: auto; }}
      .sidebar {{ border-right: 0; border-bottom: 1px solid var(--border); }}
      .sidebar-resizer {{ display: none; }}
      .tree-wrap {{ max-height: 42vh; }}
      .content {{ display: block; }}
      .gallery-scroll {{ overflow: visible; padding: 14px; }}
      .content-title-row {{ align-items: start; flex-direction: column; }}
      .view-tools {{ width: 100%; justify-content: space-between; }}
      .plot-mode-buttons {{ min-width: min(190px, 48vw); }}
      .range {{ width: min(180px, 48vw); }}
    }}
  </style>
</head>
<body>
  <div class="app-shell">
    <aside class="sidebar" aria-label="Plot folders">
      <div class="brand">
        <h1>{html.escape(title)}</h1>
        <p class="summary">{html.escape(summary)}</p>
      </div>
      <div class="toolbar">
        <input id="search" class="search" type="search" placeholder="Filter visible images" aria-label="Filter visible images">
        <div class="actions">
          <button id="expand-all" class="icon-button" type="button" title="Expand all folders" aria-label="Expand all folders">+</button>
          <button id="collapse-all" class="icon-button" type="button" title="Collapse all folders" aria-label="Collapse all folders">-</button>
          <button id="clear-selection" class="icon-button" type="button" title="Clear selected folders" aria-label="Clear selected folders">x</button>
        </div>
        <div class="mode-buttons" aria-label="Tree mode">
          <button id="source-mode" class="mode-button is-active" type="button">Source tree</button>
          <button id="output-mode" class="mode-button" type="button">Output tree</button>
        </div>
      </div>
      <div class="tree-wrap">
        <ul id="folder-tree" class="tree" data-testid="folder-tree"></ul>
      </div>
    </aside>
    <div
      id="sidebar-resizer"
      class="sidebar-resizer"
      role="separator"
      aria-label="Resize plot folder sidebar"
      aria-orientation="vertical"
      aria-valuemin="220"
      aria-valuemax="620"
      tabindex="0"
    ></div>
    <main class="content">
      <header class="content-header">
        <div class="content-title-row">
          <div>
            <h2 id="gallery-title">Selected plots</h2>
            <p id="gallery-meta" class="meta"></p>
          </div>
          <div class="view-tools">
            <div class="plot-mode-buttons" aria-label="Plot display mode">
              <button id="interactive-view" class="mode-button is-active" type="button">Interactive</button>
              <button id="static-view" class="mode-button" type="button">Static</button>
            </div>
            <label class="tile-control">
              Tile width
              <input id="tile-size" class="range" type="range" min="220" max="620" step="20" value="360">
            </label>
          </div>
        </div>
        <div id="selected-chips" class="chips" aria-label="Selected folders"></div>
      </header>
      <section class="gallery-scroll" aria-label="Plot images">
        <div id="gallery-grid" class="gallery-grid" data-testid="gallery-grid"></div>
        <div id="empty-state" class="empty-state hidden">
          <strong>No images to show</strong>
          <span>Select a folder checkbox on the left, or change the image filter.</span>
        </div>
      </section>
    </main>
  </div>
  <div id="data-modal" class="data-modal hidden" role="dialog" aria-modal="true" aria-labelledby="data-dialog-title">
    <div class="data-dialog">
      <div class="data-dialog-head">
        <div id="data-dialog-title" class="data-dialog-title">Plot data</div>
        <button id="data-close" class="icon-button" type="button" title="Close data table" aria-label="Close data table">x</button>
      </div>
      <div id="data-dialog-body" class="data-dialog-body"></div>
    </div>
  </div>
  <div id="asset-modal" class="asset-modal hidden" role="dialog" aria-modal="true" aria-labelledby="asset-dialog-title">
    <div class="asset-dialog">
      <div class="asset-dialog-head">
        <div id="asset-dialog-title" class="asset-dialog-title">Plot preview</div>
        <button id="asset-close" class="icon-button" type="button" title="Close plot preview" aria-label="Close plot preview">x</button>
      </div>
      <div id="asset-dialog-body" class="asset-dialog-body"></div>
    </div>
  </div>

  <script id="plot-data" type="application/json">{_safe_json(manifest)}</script>
  <script>
    (() => {{
      const images = JSON.parse(document.getElementById("plot-data").textContent);
      const treeEl = document.getElementById("folder-tree");
      const gridEl = document.getElementById("gallery-grid");
      const emptyEl = document.getElementById("empty-state");
      const searchEl = document.getElementById("search");
      const tileSizeEl = document.getElementById("tile-size");
      const titleEl = document.getElementById("gallery-title");
      const metaEl = document.getElementById("gallery-meta");
      const chipsEl = document.getElementById("selected-chips");
      const sourceModeEl = document.getElementById("source-mode");
      const outputModeEl = document.getElementById("output-mode");
      const interactiveViewEl = document.getElementById("interactive-view");
      const staticViewEl = document.getElementById("static-view");
      const sidebarResizerEl = document.getElementById("sidebar-resizer");
      const dataModalEl = document.getElementById("data-modal");
      const dataTitleEl = document.getElementById("data-dialog-title");
      const dataBodyEl = document.getElementById("data-dialog-body");
      const dataCloseEl = document.getElementById("data-close");
      const assetModalEl = document.getElementById("asset-modal");
      const assetTitleEl = document.getElementById("asset-dialog-title");
      const assetBodyEl = document.getElementById("asset-dialog-body");
      const assetCloseEl = document.getElementById("asset-close");
      const plotViewStorageKey = "plotGalleryViewMode";
      const sidebarWidthStorageKey = "plotGallerySidebarWidth";
      const selected = new Set([""]);
      const expanded = new Set([""]);
      let treeMode = "source";
      let plotViewMode = localStorage.getItem(plotViewStorageKey) === "static" ? "static" : "interactive";
      let folderMap = new Map();
      let root = null;
      let sidebarResizeActive = false;

      function rootLabel() {{
        return treeMode === "source" ? "Project source" : "docs/plots";
      }}

      function folderName(path) {{
        if (!path) return rootLabel();
        return path.split("/").at(-1).replaceAll("_", " ");
      }}

      function makeFolder(path, parent = null) {{
        if (folderMap.has(path)) return folderMap.get(path);
        const folder = {{ path, parent, children: new Map(), images: [], total: 0 }};
        folderMap.set(path, folder);
        return folder;
      }}

      function activeFolder(image) {{
        return treeMode === "source" ? image.source_folder : image.output_folder;
      }}

      function buildFolderMap() {{
        folderMap = new Map();
        root = makeFolder("");
        for (const image of images) {{
          const folder = activeFolder(image);
          const parts = folder ? folder.split("/") : [];
          let node = root;
          node.total += 1;
          let currentPath = "";
          for (const part of parts) {{
            currentPath = currentPath ? `${{currentPath}}/${{part}}` : part;
            const child = makeFolder(currentPath, node);
            node.children.set(part, child);
            node = child;
            node.total += 1;
          }}
          node.images.push(image);
        }}
      }}

      function sortedChildren(node) {{
        return Array.from(node.children.values()).sort((a, b) => folderName(a.path).localeCompare(folderName(b.path)));
      }}

      function renderTreeNode(node, depth = 0) {{
        const childNodes = sortedChildren(node);
        const li = document.createElement("li");
        const row = document.createElement("div");
        row.className = "folder-row";
        row.dataset.path = node.path;
        row.style.paddingLeft = `${{depth * 14}}px`;

        const toggle = document.createElement("button");
        toggle.type = "button";
        toggle.className = "toggle";
        toggle.textContent = expanded.has(node.path) ? "-" : "+";
        toggle.disabled = childNodes.length === 0;
        toggle.title = expanded.has(node.path) ? "Collapse folder" : "Expand folder";
        toggle.setAttribute("aria-label", toggle.title);
        toggle.addEventListener("click", () => {{
          if (expanded.has(node.path)) expanded.delete(node.path);
          else expanded.add(node.path);
          renderTree();
        }});

        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.className = "folder-check";
        checkbox.checked = selected.has(node.path);
        checkbox.setAttribute("aria-label", `Show ${{folderName(node.path)}}`);
        checkbox.addEventListener("change", () => {{
          if (checkbox.checked) {{
            if (!node.path) {{
              selected.clear();
              selected.add("");
            }} else {{
              selected.delete("");
              selected.add(node.path);
            }}
          }} else {{
            selected.delete(node.path);
          }}
          renderTree();
          renderGallery();
        }});

        const label = document.createElement("span");
        label.className = "folder-name";
        label.textContent = folderName(node.path);
        label.title = node.path || rootLabel();
        label.addEventListener("click", () => {{
          if (childNodes.length) {{
            if (expanded.has(node.path)) expanded.delete(node.path);
            else expanded.add(node.path);
            renderTree();
          }}
        }});

        const count = document.createElement("span");
        count.className = "folder-count";
        count.textContent = node.total;

        if (selected.has(node.path)) row.classList.add("is-selected");
        row.append(toggle, checkbox, label, count);
        li.append(row);

        if (expanded.has(node.path) && childNodes.length) {{
          const ul = document.createElement("ul");
          for (const child of childNodes) ul.append(renderTreeNode(child, depth + 1));
          li.append(ul);
        }}
        return li;
      }}

      function renderTree() {{
        buildFolderMap();
        treeEl.replaceChildren(renderTreeNode(root));
      }}

      function isInSelectedFolder(image) {{
        if (selected.size === 0) return false;
        const folderPath = activeFolder(image);
        for (const folder of selected) {{
          if (!folder || folderPath === folder || folderPath.startsWith(`${{folder}}/`)) return true;
        }}
        return false;
      }}

      function selectedLabels() {{
        return Array.from(selected)
          .sort((a, b) => a.localeCompare(b))
          .map((path) => path || rootLabel());
      }}

      function imageMatchesFilter(image, filterText) {{
        if (!filterText) return true;
        const haystack = `${{image.output_path}} ${{image.source_path}} ${{image.svg_path || ""}} ${{image.html_path || ""}} ${{image.data_path || ""}} ${{image.title}}`.toLowerCase();
        return haystack.includes(filterText);
      }}

      function renderChips(labels) {{
        chipsEl.replaceChildren();
        if (!labels.length) {{
          const chip = document.createElement("span");
          chip.className = "chip";
          chip.textContent = "No folders selected";
          chipsEl.append(chip);
          return;
        }}
        for (const label of labels) {{
          const chip = document.createElement("span");
          chip.className = "chip";
          chip.textContent = label;
          chipsEl.append(chip);
        }}
      }}

      function setPlotViewMode(mode) {{
        plotViewMode = mode === "static" ? "static" : "interactive";
        localStorage.setItem(plotViewStorageKey, plotViewMode);
        interactiveViewEl.classList.toggle("is-active", plotViewMode === "interactive");
        staticViewEl.classList.toggle("is-active", plotViewMode === "static");
        renderGallery();
      }}

      function sidebarWidthBounds() {{
        const minimum = 220;
        const maximum = Math.max(minimum, Math.min(620, Math.floor(window.innerWidth * 0.55)));
        return {{ minimum, maximum }};
      }}

      function setSidebarWidth(width, {{ persist = true }} = {{}}) {{
        const {{ minimum, maximum }} = sidebarWidthBounds();
        const clamped = Math.max(minimum, Math.min(maximum, Math.round(width)));
        document.documentElement.style.setProperty("--sidebar", `${{clamped}}px`);
        sidebarResizerEl.setAttribute("aria-valuenow", String(clamped));
        sidebarResizerEl.setAttribute("aria-valuemax", String(maximum));
        if (persist) localStorage.setItem(sidebarWidthStorageKey, String(clamped));
        return clamped;
      }}

      function restoreSidebarWidth() {{
        const stored = Number(localStorage.getItem(sidebarWidthStorageKey));
        if (Number.isFinite(stored) && stored > 0) setSidebarWidth(stored, {{ persist: false }});
        else sidebarResizerEl.setAttribute("aria-valuenow", String(Math.round(sidebarResizerEl.getBoundingClientRect().left)));
      }}

      function startSidebarResize(event) {{
        if (event.button !== undefined && event.button !== 0) return;
        sidebarResizeActive = true;
        document.body.classList.add("sidebar-resizing");
        sidebarResizerEl.setPointerCapture?.(event.pointerId);
        setSidebarWidth(event.clientX);
        event.preventDefault();
      }}

      function updateSidebarResize(event) {{
        if (!sidebarResizeActive) return;
        setSidebarWidth(event.clientX);
      }}

      function finishSidebarResize(event) {{
        if (!sidebarResizeActive) return;
        sidebarResizeActive = false;
        document.body.classList.remove("sidebar-resizing");
        if (event?.pointerId !== undefined) sidebarResizerEl.releasePointerCapture?.(event.pointerId);
      }}

      function makeAssetButton(label, path, image) {{
        const button = document.createElement("button");
        button.className = "resource-link";
        button.type = "button";
        button.textContent = label;
        button.addEventListener("click", () => showAssetPreview(label, path, image));
        return button;
      }}

      function makeDataButton(image) {{
        const button = document.createElement("button");
        button.className = "resource-link";
        button.type = "button";
        button.textContent = "Data";
        button.addEventListener("click", () => showDataTable(image));
        return button;
      }}

      function parseCsv(text) {{
        const rows = [];
        let row = [];
        let value = "";
        let quoted = false;
        for (let index = 0; index < text.length; index += 1) {{
          const char = text[index];
          const next = text[index + 1];
          if (quoted) {{
            if (char === '"' && next === '"') {{
              value += '"';
              index += 1;
            }} else if (char === '"') {{
              quoted = false;
            }} else {{
              value += char;
            }}
          }} else if (char === '"') {{
            quoted = true;
          }} else if (char === ",") {{
            row.push(value);
            value = "";
          }} else if (char === "\\n") {{
            row.push(value);
            rows.push(row);
            row = [];
            value = "";
          }} else if (char !== "\\r") {{
            value += char;
          }}
        }}
        if (value || row.length) {{
          row.push(value);
          rows.push(row);
        }}
        return rows;
      }}

      function renderDataRows(rows, path) {{
        dataBodyEl.replaceChildren();
        if (!rows.length) {{
          const note = document.createElement("p");
          note.className = "data-note";
          note.textContent = `No rows found in ${{path}}.`;
          dataBodyEl.append(note);
          return;
        }}
        const maxRows = 200;
        const [headers, ...dataRows] = rows;
        const shownRows = dataRows.slice(0, maxRows);
        const note = document.createElement("p");
        note.className = "data-note";
        note.textContent = `Showing ${{shownRows.length}} of ${{dataRows.length}} data row${{dataRows.length === 1 ? "" : "s"}} from ${{path}}.`;
        const table = document.createElement("table");
        table.className = "data-table";
        const thead = document.createElement("thead");
        const headRow = document.createElement("tr");
        for (const header of headers) {{
          const th = document.createElement("th");
          th.textContent = header;
          headRow.append(th);
        }}
        thead.append(headRow);
        const tbody = document.createElement("tbody");
        for (const sourceRow of shownRows) {{
          const tr = document.createElement("tr");
          for (let column = 0; column < headers.length; column += 1) {{
            const td = document.createElement("td");
            td.textContent = sourceRow[column] ?? "";
            tr.append(td);
          }}
          tbody.append(tr);
        }}
        table.append(thead, tbody);
        dataBodyEl.append(note, table);
      }}

      async function showDataTable(image) {{
        dataTitleEl.textContent = image.title;
        dataBodyEl.replaceChildren();
        const loading = document.createElement("p");
        loading.className = "data-note";
        loading.textContent = `Loading ${{image.data_path}}...`;
        dataBodyEl.append(loading);
        dataModalEl.classList.remove("hidden");
        try {{
          const response = await fetch(image.data_path);
          if (!response.ok) throw new Error(`HTTP ${{response.status}}`);
          renderDataRows(parseCsv(await response.text()), image.data_path);
        }} catch (error) {{
          dataBodyEl.replaceChildren();
          const note = document.createElement("p");
          note.className = "data-note";
          note.textContent = `Could not load ${{image.data_path}}: ${{error.message}}`;
          dataBodyEl.append(note);
        }}
      }}

      function closeDataTable() {{
        dataModalEl.classList.add("hidden");
        dataBodyEl.replaceChildren();
      }}

      function showAssetPreview(label, path, image) {{
        assetTitleEl.textContent = `${{label}} preview: ${{image.title}}`;
        assetBodyEl.replaceChildren();
        if (label === "SVG") {{
          const frame = document.createElement("iframe");
          frame.className = "asset-preview-frame";
          frame.title = `${{image.title}} SVG`;
          frame.src = path;
          assetBodyEl.append(frame);
        }} else {{
          const img = document.createElement("img");
          img.className = "asset-preview-image";
          img.src = path;
          img.alt = image.title;
          assetBodyEl.append(img);
        }}
        assetModalEl.classList.remove("hidden");
      }}

      function closeAssetPreview() {{
        assetModalEl.classList.add("hidden");
        assetBodyEl.replaceChildren();
      }}

      function renderPlotPreview(image) {{
        const preview = document.createElement("div");
        preview.className = "plot-preview";
        if (plotViewMode === "interactive" && image.html_path) {{
          preview.classList.add("is-interactive");
          const frame = document.createElement("iframe");
          frame.className = "interactive-frame";
          frame.loading = "lazy";
          frame.title = image.title;
          frame.src = image.html_path;
          preview.append(frame);
          return preview;
        }}

        const img = document.createElement("img");
        img.loading = "lazy";
        img.decoding = "async";
        img.src = image.output_path;
        img.alt = image.title;
        preview.append(img);
        return preview;
      }}

      function renderGallery() {{
        const filterText = searchEl.value.trim().toLowerCase();
        const visible = images
          .filter((image) => isInSelectedFolder(image) && imageMatchesFilter(image, filterText))
          .sort((a, b) => {{
            if (plotViewMode === "interactive" && Boolean(a.html_path) !== Boolean(b.html_path)) {{
              return a.html_path ? -1 : 1;
            }}
            return a.source_path.localeCompare(b.source_path);
          }});
        const labels = selectedLabels();

        gridEl.replaceChildren();
        renderChips(labels);
        titleEl.textContent = labels.length === 1 ? labels[0] : "Selected plot folders";
        metaEl.textContent = `${{visible.length}} visible image${{visible.length === 1 ? "" : "s"}} from ${{selected.size}} selected folder${{selected.size === 1 ? "" : "s"}}`;
        emptyEl.classList.toggle("hidden", visible.length !== 0);

        const fragment = document.createDocumentFragment();
        for (const image of visible) {{
          const card = document.createElement("article");
          card.className = "image-card";

          const head = document.createElement("div");
          head.className = "image-head";
          const imageTitle = document.createElement("div");
          imageTitle.className = "image-title";
          imageTitle.textContent = image.title;
          const imagePath = document.createElement("div");
          imagePath.className = "image-path";
          imagePath.textContent = image.source_path;
          imagePath.title = `Output: ${{image.output_path}}`;
          const actions = document.createElement("div");
          actions.className = "image-actions";
          if (image.output_path) actions.append(makeAssetButton("PNG", image.output_path, image));
          if (image.svg_path) actions.append(makeAssetButton("SVG", image.svg_path, image));
          if (image.data_path) actions.append(makeDataButton(image));
          const badge = document.createElement("span");
          badge.className = "asset-badge";
          if (image.interactive_source === "csv_backfill") {{
            badge.textContent = "CSV interactive";
            badge.title = "CSV-backed interactive reconstruction from plot data";
          }} else if (image.interactive_source === "native_plotly") {{
            badge.textContent = "Interactive";
          }} else {{
            badge.textContent = "Static only";
          }}
          actions.append(badge);
          head.append(imageTitle, imagePath, actions);

          card.append(head, renderPlotPreview(image));
          fragment.append(card);
        }}
        gridEl.append(fragment);
      }}

      document.getElementById("expand-all").addEventListener("click", () => {{
        for (const folder of folderMap.values()) expanded.add(folder.path);
        renderTree();
      }});
      document.getElementById("collapse-all").addEventListener("click", () => {{
        expanded.clear();
        expanded.add("");
        renderTree();
      }});
      document.getElementById("clear-selection").addEventListener("click", () => {{
        selected.clear();
        renderTree();
        renderGallery();
      }});
      function setTreeMode(mode) {{
        treeMode = mode;
        selected.clear();
        selected.add("");
        expanded.clear();
        expanded.add("");
        sourceModeEl.classList.toggle("is-active", treeMode === "source");
        outputModeEl.classList.toggle("is-active", treeMode === "output");
        renderTree();
        renderGallery();
      }}
      sourceModeEl.addEventListener("click", () => setTreeMode("source"));
      outputModeEl.addEventListener("click", () => setTreeMode("output"));
      interactiveViewEl.addEventListener("click", () => setPlotViewMode("interactive"));
      staticViewEl.addEventListener("click", () => setPlotViewMode("static"));
      sidebarResizerEl.addEventListener("pointerdown", startSidebarResize);
      sidebarResizerEl.addEventListener("pointermove", updateSidebarResize);
      sidebarResizerEl.addEventListener("pointerup", finishSidebarResize);
      sidebarResizerEl.addEventListener("pointercancel", finishSidebarResize);
      sidebarResizerEl.addEventListener("dblclick", () => {{
        localStorage.removeItem(sidebarWidthStorageKey);
        document.documentElement.style.removeProperty("--sidebar");
        restoreSidebarWidth();
      }});
      sidebarResizerEl.addEventListener("keydown", (event) => {{
        if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
        const current = Number(sidebarResizerEl.getAttribute("aria-valuenow")) || sidebarResizerEl.getBoundingClientRect().left;
        const delta = event.key === "ArrowLeft" ? -24 : 24;
        setSidebarWidth(current + delta);
        event.preventDefault();
      }});
      window.addEventListener("resize", () => {{
        const current = Number(sidebarResizerEl.getAttribute("aria-valuenow"));
        if (Number.isFinite(current) && current > 0) setSidebarWidth(current, {{ persist: false }});
      }});
      dataCloseEl.addEventListener("click", closeDataTable);
      assetCloseEl.addEventListener("click", closeAssetPreview);
      dataModalEl.addEventListener("click", (event) => {{
        if (event.target === dataModalEl) closeDataTable();
      }});
      assetModalEl.addEventListener("click", (event) => {{
        if (event.target === assetModalEl) closeAssetPreview();
      }});
      document.addEventListener("keydown", (event) => {{
        if (event.key === "Escape" && !dataModalEl.classList.contains("hidden")) closeDataTable();
        if (event.key === "Escape" && !assetModalEl.classList.contains("hidden")) closeAssetPreview();
      }});
      searchEl.addEventListener("input", renderGallery);
      tileSizeEl.addEventListener("input", () => {{
        gridEl.style.setProperty("--card-min", `${{tileSizeEl.value}}px`);
      }});

      interactiveViewEl.classList.toggle("is-active", plotViewMode === "interactive");
      staticViewEl.classList.toggle("is-active", plotViewMode === "static");
      restoreSidebarWidth();
      gridEl.style.setProperty("--card-min", `${{tileSizeEl.value}}px`);
      renderTree();
      renderGallery();
    }})();
  </script>
</body>
</html>
"""


def remove_stale_indexes() -> None:
    for index_path in PLOTS_ROOT.rglob("index.html"):
        if index_path != PLOTS_ROOT / "index.html":
            index_path.unlink()


def main() -> None:
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    pngs = collect_pngs(PLOTS_ROOT)
    remove_stale_indexes()
    (PLOTS_ROOT / "index.html").write_text(render_gallery_page(PLOTS_ROOT, pngs), encoding="utf-8")


if __name__ == "__main__":
    main()
