from __future__ import annotations

import argparse
import socket
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PLOTS_ROOT = REPO_ROOT / "docs" / "plots"
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


def serve_gallery(host: str, port: int) -> None:
    handler = partial(SimpleHTTPRequestHandler, directory=str(PLOTS_ROOT))
    server = ThreadingHTTPServer((host, port), handler)
    url = f"http://{host}:{port}/"
    print(f"Serving ePC-SAFT plot gallery at {url}", flush=True)
    print(f"Root index: {PLOTS_ROOT / 'index.html'}", flush=True)
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
    serve_gallery(args.host, port)


if __name__ == "__main__":
    main()
