from __future__ import annotations

import argparse
import socket
import sys
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
GALLERY_INDEX = PLOTS_ROOT / "index.html"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def find_available_port(host: str, preferred_port: int, *, attempts: int = 50) -> int:
    for port in range(preferred_port, preferred_port + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                probe.bind((host, port))
            except OSError:
                continue
            return port
    raise RuntimeError(f"No available localhost port found from {preferred_port} to {preferred_port + attempts - 1}.")


def build_gallery() -> None:
    from scripts.paper_validation.tools import build_analysis_galleries

    build_analysis_galleries.main()


class PlotGalleryRequestHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, auto_build: bool = True, **kwargs) -> None:
        self.auto_build = auto_build
        super().__init__(*args, **kwargs)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def send_head(self):
        request_path = self.path.split("?", 1)[0]
        if self.auto_build and request_path in {"/", "/index.html", "/docs/plots/", "/docs/plots/index.html"}:
            try:
                build_gallery()
            except Exception as exc:  # pragma: no cover - exercised through manual server use.
                self.send_error(HTTPStatus.INTERNAL_SERVER_ERROR, f"Could not rebuild plot gallery: {exc}")
                return None
            if request_path in {"/", "/index.html"}:
                self.path = "/docs/plots/index.html"
        return super().send_head()


def serve_gallery(host: str, port: int, *, auto_build: bool = True) -> None:
    handler = partial(PlotGalleryRequestHandler, directory=str(REPO_ROOT), auto_build=auto_build)
    server = ThreadingHTTPServer((host, port), handler)
    url = f"http://{host}:{port}/docs/plots/index.html"
    print(f"Serving ePC-SAFT plot gallery at {url}", flush=True)
    print(f"Root index: {GALLERY_INDEX}", flush=True)
    if auto_build:
        print("Auto-rebuilding root index on browser refresh.", flush=True)
    print("Press Ctrl+C to stop.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping plot gallery server.", flush=True)
    finally:
        server.server_close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and serve the generated docs/plots gallery on localhost.")
    parser.add_argument("--host", default="127.0.0.1", help="Host interface for the local server.")
    parser.add_argument("--port", type=int, default=8765, help="Preferred local server port.")
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Serve the existing gallery without rebuilding docs/plots/**/index.html first.",
    )
    parser.add_argument(
        "--no-auto-build",
        action="store_true",
        help="Do not rebuild docs/plots/index.html when the browser requests the root page.",
    )
    parser.add_argument(
        "--strict-port",
        action="store_true",
        help="Fail instead of choosing the next free port when the preferred port is busy.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.no_build:
        build_gallery()
    port = args.port if args.strict_port else find_available_port(args.host, args.port)
    serve_gallery(args.host, port, auto_build=not args.no_auto_build)


if __name__ == "__main__":
    main()
