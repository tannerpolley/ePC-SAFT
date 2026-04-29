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
        source_path = _source_path_for_output(output_path)
        manifest.append(
            {
                "path": output_path,
                "folder": _folder_for(source_path),
                "output_path": output_path,
                "output_folder": _folder_for(output_path),
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
    .app-shell {{
      display: grid;
      grid-template-columns: var(--sidebar) minmax(0, 1fr);
      height: 100vh;
    }}
    .sidebar {{
      display: flex;
      min-width: 0;
      flex-direction: column;
      border-right: 1px solid var(--border);
      background: var(--panel);
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
    .image-link {{
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
    .image-link:focus-visible {{
      outline: 3px solid rgba(25, 95, 184, 0.35);
      outline-offset: -3px;
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
      .tree-wrap {{ max-height: 42vh; }}
      .content {{ display: block; }}
      .gallery-scroll {{ overflow: visible; padding: 14px; }}
      .content-title-row {{ align-items: start; flex-direction: column; }}
      .view-tools {{ width: 100%; justify-content: space-between; }}
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
    <main class="content">
      <header class="content-header">
        <div class="content-title-row">
          <div>
            <h2 id="gallery-title">Selected plots</h2>
            <p id="gallery-meta" class="meta"></p>
          </div>
          <label class="view-tools">
            Tile width
            <input id="tile-size" class="range" type="range" min="220" max="520" step="20" value="320">
          </label>
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
      const selected = new Set([""]);
      const expanded = new Set([""]);
      let treeMode = "source";
      let folderMap = new Map();
      let root = null;

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
        const haystack = `${{image.output_path}} ${{image.source_path}} ${{image.title}}`.toLowerCase();
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

      function renderGallery() {{
        const filterText = searchEl.value.trim().toLowerCase();
        const visible = images.filter((image) => isInSelectedFolder(image) && imageMatchesFilter(image, filterText));
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
          head.append(imageTitle, imagePath);

          const link = document.createElement("a");
          link.className = "image-link";
          link.href = image.output_path;
          link.target = "_blank";
          link.rel = "noopener";

          const img = document.createElement("img");
          img.loading = "lazy";
          img.decoding = "async";
          img.src = image.output_path;
          img.alt = image.title;
          link.append(img);

          card.append(head, link);
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
      searchEl.addEventListener("input", renderGallery);
      tileSizeEl.addEventListener("input", () => {{
        gridEl.style.setProperty("--card-min", `${{tileSizeEl.value}}px`);
      }});

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
