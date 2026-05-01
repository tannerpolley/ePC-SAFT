from __future__ import annotations

import html
import json
from pathlib import Path


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


def image_manifest(pngs: list[Path]) -> list[dict[str, str]]:
    manifest = []
    for png in pngs:
        output_path = png.relative_to(PLOTS_ROOT).as_posix()
        svg = png.with_suffix(".svg")
        data = png.parent / "data" / f"{png.stem}_plot_data.csv"
        source_path = _source_path_for_output(output_path)
        manifest.append(
            {
                "path": output_path,
                "folder": _folder_for(source_path),
                "output_path": output_path,
                "output_folder": _folder_for(output_path),
                "svg_path": svg.relative_to(PLOTS_ROOT).as_posix() if svg.exists() else "",
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
      --shadow: 0 8px 24px rgba(25, 43, 65, 0.08);
      --sidebar: clamp(280px, 27vw, 420px);
      --sidebar-resizer: 8px;
      --card-min: 340px;
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
    body.sidebar-resizing {{ cursor: col-resize; user-select: none; }}
    .app-shell {{
      display: grid;
      grid-template-columns: var(--sidebar) var(--sidebar-resizer) minmax(0, 1fr);
      height: 100vh;
    }}
    .sidebar {{ display: flex; min-width: 0; flex-direction: column; background: var(--panel); }}
    .sidebar-resizer {{
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
    .brand {{ padding: 18px 18px 14px; border-bottom: 1px solid var(--border); }}
    h1 {{ margin: 0; font-size: 1.18rem; line-height: 1.2; letter-spacing: 0; }}
    .summary {{ margin: 7px 0 0; color: var(--muted); font-size: 0.88rem; }}
    .toolbar {{ display: grid; gap: 10px; padding: 12px 14px; border-bottom: 1px solid var(--border); background: var(--panel-2); }}
    .search {{ width: 100%; min-height: 34px; border: 1px solid var(--border); border-radius: 6px; padding: 7px 10px; background: #fff; color: var(--text); }}
    .actions, .mode-buttons {{ display: flex; gap: 8px; align-items: center; flex-wrap: wrap; }}
    .mode-button, .icon-button {{
      min-height: 32px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #fff;
      color: var(--muted);
      cursor: pointer;
      font-size: 0.82rem;
    }}
    .mode-button {{ padding: 0 10px; }}
    .mode-button.is-active {{ border-color: var(--blue); background: var(--blue-soft); color: var(--blue); font-weight: 650; }}
    .icon-button {{ width: 34px; height: 34px; color: var(--text); }}
    .tree-wrap {{ min-height: 0; overflow: auto; padding: 10px 8px 18px; }}
    .tree, .tree ul {{ margin: 0; padding: 0; list-style: none; font-size: 0.91rem; }}
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
    .toggle {{ width: 22px; height: 22px; border: 0; border-radius: 4px; background: transparent; color: var(--muted); cursor: pointer; }}
    .toggle:disabled {{ cursor: default; opacity: 0.2; }}
    .folder-check {{ width: 15px; height: 15px; accent-color: var(--blue); }}
    .folder-name {{ overflow: hidden; text-overflow: ellipsis; white-space: nowrap; cursor: pointer; }}
    .folder-count {{ color: var(--muted); font-size: 0.78rem; font-variant-numeric: tabular-nums; }}
    .content {{ min-width: 0; min-height: 0; overflow: auto; padding: 18px; }}
    .content-head {{ display: flex; justify-content: space-between; align-items: center; gap: 16px; margin-bottom: 14px; }}
    .content-title {{ margin: 0; font-size: 1rem; }}
    .tile-control {{ display: flex; gap: 8px; align-items: center; color: var(--muted); font-size: 0.84rem; }}
    .gallery-grid {{
      --card-min-local: var(--card-min);
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(var(--card-min-local), 1fr));
      gap: 16px;
    }}
    .plot-card {{ min-width: 0; border: 1px solid var(--border); border-radius: var(--radius); background: var(--panel); box-shadow: var(--shadow); overflow: hidden; }}
    .plot-head {{ display: grid; gap: 5px; padding: 12px 12px 10px; border-bottom: 1px solid var(--border); }}
    .plot-title {{ margin: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-size: 0.95rem; }}
    .plot-path {{ color: var(--muted); font-size: 0.78rem; overflow-wrap: anywhere; }}
    .image-actions {{ display: flex; gap: 8px; align-items: center; flex-wrap: wrap; margin-top: 4px; }}
    .asset-button {{
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #fff;
      color: var(--blue);
      cursor: pointer;
      padding: 4px 8px;
      font-size: 0.8rem;
    }}
    .asset-button:disabled {{ color: var(--muted); cursor: default; background: #f4f6f8; }}
    .plot-preview {{ display: grid; place-items: center; min-height: 220px; padding: 10px; background: #f7f9fc; }}
    .plot-preview img {{ display: block; max-width: 100%; max-height: 520px; object-fit: contain; background: #fff; border: 1px solid var(--border); border-radius: 6px; }}
    .empty-state {{ padding: 28px; color: var(--muted); text-align: center; border: 1px dashed var(--border-strong); border-radius: var(--radius); background: #fff; }}
    .asset-modal, .data-modal {{
      position: fixed;
      inset: 0;
      z-index: 10;
      display: grid;
      place-items: center;
      padding: 20px;
      background: rgba(13, 25, 43, 0.55);
    }}
    .asset-modal.hidden, .data-modal.hidden {{ display: none; }}
    .dialog {{
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      width: min(1180px, 96vw);
      max-height: 92vh;
      background: #fff;
      border-radius: 8px;
      overflow: hidden;
      box-shadow: 0 18px 48px rgba(11, 22, 37, 0.28);
    }}
    .dialog-head {{ display: flex; justify-content: space-between; gap: 12px; align-items: center; padding: 12px 14px; border-bottom: 1px solid var(--border); }}
    .dialog-title {{ margin: 0; font-size: 0.98rem; }}
    .dialog-body {{ min-height: 0; overflow: auto; padding: 14px; }}
    .dialog-close {{ border: 1px solid var(--border); border-radius: 6px; background: #fff; cursor: pointer; padding: 4px 8px; }}
    .asset-preview {{ display: grid; place-items: center; min-height: 50vh; }}
    .asset-preview img {{ max-width: 100%; max-height: 78vh; object-fit: contain; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.82rem; }}
    th, td {{ border: 1px solid var(--border); padding: 4px 6px; text-align: left; vertical-align: top; }}
    th {{ position: sticky; top: 0; background: #eef4fb; }}
  </style>
</head>
<body>
  <div class="app-shell">
    <aside class="sidebar">
      <div class="brand">
        <h1>{html.escape(title)}</h1>
        <p class="summary">{html.escape(summary)}</p>
      </div>
      <div class="toolbar">
        <input id="search" class="search" type="search" placeholder="Search plots, folders, or source paths">
        <div class="mode-buttons" aria-label="Folder tree mode">
          <button id="source-mode" class="mode-button is-active" type="button">Source tree</button>
          <button id="output-mode" class="mode-button" type="button">Output tree</button>
        </div>
        <div class="actions">
          <button id="expand-all" class="icon-button" type="button" title="Expand all">+</button>
          <button id="collapse-all" class="icon-button" type="button" title="Collapse all">-</button>
          <button id="clear-selection" class="mode-button" type="button">Clear</button>
        </div>
      </div>
      <div class="tree-wrap"><ul id="folder-tree" class="tree" data-testid="folder-tree"></ul></div>
    </aside>
    <div id="sidebar-resizer" class="sidebar-resizer" role="separator" aria-label="Resize plot folder sidebar" tabindex="0"></div>
    <main class="content">
      <div class="content-head">
        <h2 id="content-title" class="content-title">Select folders on the left</h2>
        <label class="tile-control">Tile size <input id="tile-size" type="range" min="260" max="560" value="340"></label>
      </div>
      <div id="gallery-grid" class="gallery-grid" data-testid="gallery-grid"></div>
    </main>
  </div>
  <div id="data-modal" class="data-modal hidden" role="dialog" aria-modal="true" aria-labelledby="data-dialog-title">
    <div class="dialog">
      <div class="dialog-head"><h2 id="data-dialog-title" class="dialog-title">CSV data</h2><button id="data-close" class="dialog-close" type="button">Close</button></div>
      <div id="data-body" class="dialog-body"></div>
    </div>
  </div>
  <div id="asset-modal" class="asset-modal hidden" role="dialog" aria-modal="true" aria-labelledby="asset-dialog-title">
    <div class="dialog">
      <div class="dialog-head"><h2 id="asset-dialog-title" class="dialog-title">Asset preview</h2><button id="asset-close" class="dialog-close" type="button">Close</button></div>
      <div id="asset-body" class="dialog-body asset-preview"></div>
    </div>
  </div>
  <script>
    (() => {{
      const images = {_safe_json(manifest)};
      const selected = new Set();
      const expanded = new Set([""]);
      const folderTreeEl = document.getElementById("folder-tree");
      const galleryEl = document.getElementById("gallery-grid");
      const searchEl = document.getElementById("search");
      const titleEl = document.getElementById("content-title");
      const tileSizeEl = document.getElementById("tile-size");
      const sourceModeEl = document.getElementById("source-mode");
      const outputModeEl = document.getElementById("output-mode");
      const dataModalEl = document.getElementById("data-modal");
      const dataBodyEl = document.getElementById("data-body");
      const dataCloseEl = document.getElementById("data-close");
      const assetModalEl = document.getElementById("asset-modal");
      const assetBodyEl = document.getElementById("asset-body");
      const assetCloseEl = document.getElementById("asset-close");
      const sidebarResizerEl = document.getElementById("sidebar-resizer");
      const sidebarWidthStorageKey = "plotGallerySidebarWidth";
      const galleryStateStorageKey = "plotGalleryStaticStateV1";
      let treeMode = "source";
      let resizing = false;

      function folderPath(image) {{ return treeMode === "source" ? image.source_folder : image.output_folder; }}
      function makeFolderMap() {{
        const map = new Map([["", {{ path: "", name: "All plots", children: new Set(), count: images.length }}]]);
        for (const image of images) {{
          const folder = folderPath(image);
          const parts = folder ? folder.split("/") : [];
          let current = "";
          for (const part of parts) {{
            const parent = current;
            current = current ? `${{current}}/${{part}}` : part;
            if (!map.has(current)) map.set(current, {{ path: current, name: part, children: new Set(), count: 0 }});
            map.get(parent).children.add(current);
          }}
          const allFolders = [""];
          current = "";
          for (const part of parts) {{
            current = current ? `${{current}}/${{part}}` : part;
            allFolders.push(current);
          }}
          for (const path of allFolders) map.get(path).count += path === "" ? 0 : 1;
        }}
        return map;
      }}

      function saveGalleryState() {{
        const state = {{
          treeMode,
          selected: [...selected],
          expanded: [...expanded],
          search: searchEl.value,
          tileSize: tileSizeEl.value,
          contentScrollTop: document.querySelector(".content").scrollTop,
        }};
        localStorage.setItem(galleryStateStorageKey, JSON.stringify(state));
      }}

      function expandAncestors(path) {{
        const parts = path ? path.split("/") : [];
        let current = "";
        expanded.add("");
        for (const part of parts.slice(0, -1)) {{
          current = current ? `${{current}}/${{part}}` : part;
          expanded.add(current);
        }}
      }}

      function restoreGalleryState() {{
        let state = null;
        try {{
          state = JSON.parse(localStorage.getItem(galleryStateStorageKey) || "null");
        }} catch (_error) {{
          state = null;
        }}
        if (!state || typeof state !== "object") return;
        if (state.treeMode === "source" || state.treeMode === "output") treeMode = state.treeMode;
        if (typeof state.search === "string") searchEl.value = state.search;
        if (typeof state.tileSize === "string" && state.tileSize) tileSizeEl.value = state.tileSize;
        const map = makeFolderMap();
        if (Array.isArray(state.expanded)) {{
          expanded.clear();
          expanded.add("");
          for (const path of state.expanded) if (typeof path === "string" && map.has(path)) expanded.add(path);
        }}
        if (Array.isArray(state.selected)) {{
          selected.clear();
          for (const path of state.selected) {{
            if (typeof path === "string" && map.has(path)) {{
              selected.add(path);
              expandAncestors(path);
            }}
          }}
        }}
        sourceModeEl.classList.toggle("is-active", treeMode === "source");
        outputModeEl.classList.toggle("is-active", treeMode === "output");
        if (Number.isFinite(Number(state.contentScrollTop))) {{
          requestAnimationFrame(() => {{ document.querySelector(".content").scrollTop = Number(state.contentScrollTop); }});
        }}
      }}

      function descendants(folder) {{
        const map = makeFolderMap();
        const result = [];
        for (const path of map.keys()) {{
          if (path === folder || (folder && path.startsWith(`${{folder}}/`)) || folder === "") result.push(path);
        }}
        return result;
      }}

      function renderTree() {{
        const map = makeFolderMap();
        function renderNode(path) {{
          const node = map.get(path);
          const li = document.createElement("li");
          const row = document.createElement("div");
          row.className = "folder-row";
          row.style.paddingLeft = `${{path ? path.split("/").length * 12 : 0}}px`;
          if (selected.has(path)) row.classList.add("is-selected");
          const toggle = document.createElement("button");
          toggle.className = "toggle";
          toggle.type = "button";
          toggle.textContent = node.children.size ? (expanded.has(path) ? "▾" : "▸") : "";
          toggle.disabled = !node.children.size;
          toggle.addEventListener("click", () => {{
            if (expanded.has(path)) expanded.delete(path);
            else expanded.add(path);
            renderTree();
          }});
          const check = document.createElement("input");
          check.className = "folder-check";
          check.type = "checkbox";
          check.checked = selected.has(path);
          check.addEventListener("change", () => {{
            if (check.checked) selected.add(path);
            else selected.delete(path);
            expandAncestors(path);
            saveGalleryState();
            renderTree();
            renderGallery();
          }});
          const name = document.createElement("span");
          name.className = "folder-name";
          name.textContent = node.name;
          name.title = path || "All plots";
          name.addEventListener("click", () => {{
            if (selected.has(path)) selected.delete(path);
            else selected.add(path);
            expandAncestors(path);
            saveGalleryState();
            renderTree();
            renderGallery();
          }});
          const count = document.createElement("span");
          count.className = "folder-count";
          count.textContent = node.count;
          row.append(toggle, check, name, count);
          li.append(row);
          if (expanded.has(path) && node.children.size) {{
            const ul = document.createElement("ul");
            for (const child of [...node.children].sort()) ul.append(renderNode(child));
            li.append(ul);
          }}
          return li;
        }}
        folderTreeEl.replaceChildren(renderNode(""));
      }}

      function selectedFolders() {{ return selected.size ? [...selected] : []; }}
      function imageMatchesFolder(image, folder) {{
        const path = folderPath(image);
        return folder === "" || path === folder || path.startsWith(`${{folder}}/`);
      }}
      function visibleImages() {{
        const query = searchEl.value.trim().toLowerCase();
        const folders = selectedFolders();
        return images.filter((image) => {{
          const inFolder = !folders.length || folders.some((folder) => imageMatchesFolder(image, folder));
          const haystack = `${{image.output_path}} ${{image.source_path}} ${{image.svg_path}} ${{image.data_path}} ${{image.title}}`.toLowerCase();
          return inFolder && (!query || haystack.includes(query));
        }});
      }}

      function makeAssetButton(label, path, image) {{
        const button = document.createElement("button");
        button.className = "asset-button";
        button.type = "button";
        button.textContent = label;
        button.disabled = !path;
        button.addEventListener("click", () => {{
          if (!path) return;
          showAssetPreview(label, path, image);
        }});
        return button;
      }}
      function makeDataButton(image) {{
        const button = document.createElement("button");
        button.className = "asset-button";
        button.type = "button";
        button.textContent = "Data";
        button.disabled = !image.data_path;
        button.addEventListener("click", () => showDataTable(image));
        return button;
      }}
      function parseCsv(text) {{
        const rows = [];
        let row = [], field = "", quoted = false;
        for (let i = 0; i < text.length; i++) {{
          const char = text[i];
          if (quoted) {{
            if (char === '"' && text[i + 1] === '"') {{ field += '"'; i++; }}
            else if (char === '"') quoted = false;
            else field += char;
          }} else if (char === '"') quoted = true;
          else if (char === ",") {{ row.push(field); field = ""; }}
          else if (char === "\\n") {{ row.push(field); rows.push(row); row = []; field = ""; }}
          else if (char !== "\\r") field += char;
        }}
        if (field || row.length) {{ row.push(field); rows.push(row); }}
        return rows;
      }}
      async function showDataTable(image) {{
        if (!image.data_path) return;
        const response = await fetch(image.data_path);
        const text = await response.text();
        const rows = parseCsv(text);
        const table = document.createElement("table");
        for (const [index, cells] of rows.entries()) {{
          const tr = document.createElement("tr");
          for (const cell of cells) {{
            const el = document.createElement(index === 0 ? "th" : "td");
            el.textContent = cell;
            tr.append(el);
          }}
          table.append(tr);
        }}
        document.getElementById("data-dialog-title").textContent = image.data_path;
        dataBodyEl.replaceChildren(table);
        dataModalEl.classList.remove("hidden");
      }}
      function closeDataTable() {{ dataModalEl.classList.add("hidden"); dataBodyEl.replaceChildren(); }}
      function showAssetPreview(label, path, image) {{
        document.getElementById("asset-dialog-title").textContent = `${{label}}: ${{path}}`;
        assetBodyEl.replaceChildren();
        if (label === "PNG" || label === "SVG") {{
          const img = document.createElement("img");
          img.src = path;
          img.alt = image.title;
          assetBodyEl.append(img);
        }}
        assetModalEl.classList.remove("hidden");
      }}
      function closeAssetPreview() {{ assetModalEl.classList.add("hidden"); assetBodyEl.replaceChildren(); }}

      function renderGallery() {{
        const items = visibleImages();
        galleryEl.style.setProperty("--card-min-local", `${{tileSizeEl.value}}px`);
        titleEl.textContent = items.length ? `${{items.length}} plot${{items.length === 1 ? "" : "s"}}` : "No matching plots";
        if (!items.length) {{
          const empty = document.createElement("div");
          empty.className = "empty-state";
          empty.textContent = "Select folders on the left or adjust the search filter.";
          galleryEl.replaceChildren(empty);
          return;
        }}
        const fragment = document.createDocumentFragment();
        for (const image of items.sort((a, b) => a.output_path.localeCompare(b.output_path))) {{
          const card = document.createElement("article");
          card.className = "plot-card";
          const head = document.createElement("div");
          head.className = "plot-head";
          const imageTitle = document.createElement("h3");
          imageTitle.className = "plot-title";
          imageTitle.textContent = image.title;
          const imagePath = document.createElement("div");
          imagePath.className = "plot-path";
          imagePath.textContent = image.output_path;
          const actions = document.createElement("div");
          actions.className = "image-actions";
          actions.append(makeAssetButton("PNG", image.output_path, image));
          actions.append(makeAssetButton("SVG", image.svg_path, image));
          actions.append(makeDataButton(image));
          const preview = document.createElement("div");
          preview.className = "plot-preview";
          const img = document.createElement("img");
          img.src = image.output_path;
          img.alt = image.title;
          preview.append(img);
          head.append(imageTitle, imagePath, actions);
          card.append(head, preview);
          fragment.append(card);
        }}
        galleryEl.replaceChildren(fragment);
      }}

      function setTreeMode(mode) {{
        treeMode = mode;
        selected.clear();
        expanded.clear();
        expanded.add("");
        sourceModeEl.classList.toggle("is-active", treeMode === "source");
        outputModeEl.classList.toggle("is-active", treeMode === "output");
        saveGalleryState();
        renderTree();
        renderGallery();
      }}
      function restoreSidebarWidth() {{
        const saved = Number(localStorage.getItem(sidebarWidthStorageKey));
        if (Number.isFinite(saved) && saved > 180) document.documentElement.style.setProperty("--sidebar", `${{saved}}px`);
      }}
      function setSidebarWidth(clientX, persist = true) {{
        const shellLeft = document.querySelector(".app-shell").getBoundingClientRect().left;
        const width = Math.min(Math.max(clientX - shellLeft, 220), Math.min(window.innerWidth * 0.55, 620));
        document.documentElement.style.setProperty("--sidebar", `${{width}}px`);
        sidebarResizerEl.setAttribute("aria-valuenow", String(Math.round(width)));
        if (persist) localStorage.setItem(sidebarWidthStorageKey, String(Math.round(width)));
      }}
      function startSidebarResize(event) {{ resizing = true; document.body.classList.add("sidebar-resizing"); sidebarResizerEl.setPointerCapture(event.pointerId); }}
      function updateSidebarResize(event) {{ if (resizing) setSidebarWidth(event.clientX); }}
      function finishSidebarResize(event) {{ if (!resizing) return; resizing = false; document.body.classList.remove("sidebar-resizing"); sidebarResizerEl.releasePointerCapture(event.pointerId); }}

      document.getElementById("expand-all").addEventListener("click", () => {{ for (const folder of makeFolderMap().keys()) expanded.add(folder); saveGalleryState(); renderTree(); }});
      document.getElementById("collapse-all").addEventListener("click", () => {{ expanded.clear(); expanded.add(""); saveGalleryState(); renderTree(); }});
      document.getElementById("clear-selection").addEventListener("click", () => {{ selected.clear(); saveGalleryState(); renderTree(); renderGallery(); }});
      sourceModeEl.addEventListener("click", () => setTreeMode("source"));
      outputModeEl.addEventListener("click", () => setTreeMode("output"));
      sidebarResizerEl.addEventListener("pointerdown", startSidebarResize);
      sidebarResizerEl.addEventListener("pointermove", updateSidebarResize);
      sidebarResizerEl.addEventListener("pointerup", finishSidebarResize);
      sidebarResizerEl.addEventListener("pointercancel", finishSidebarResize);
      searchEl.addEventListener("input", () => {{ saveGalleryState(); renderGallery(); }});
      tileSizeEl.addEventListener("input", () => {{ saveGalleryState(); renderGallery(); }});
      document.querySelector(".content").addEventListener("scroll", () => {{
        clearTimeout(window.plotGalleryScrollSaveTimer);
        window.plotGalleryScrollSaveTimer = setTimeout(saveGalleryState, 120);
      }});
      dataCloseEl.addEventListener("click", closeDataTable);
      assetCloseEl.addEventListener("click", closeAssetPreview);
      dataModalEl.addEventListener("click", (event) => {{ if (event.target === dataModalEl) closeDataTable(); }});
      assetModalEl.addEventListener("click", (event) => {{ if (event.target === assetModalEl) closeAssetPreview(); }});
      document.addEventListener("keydown", (event) => {{
        if (event.key === "Escape") {{ closeDataTable(); closeAssetPreview(); }}
      }});
      restoreSidebarWidth();
      restoreGalleryState();
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
