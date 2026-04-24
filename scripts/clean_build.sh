#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

rm -rf build dist .pytest_cache .ruff_cache .mypy_cache 2>/dev/null || true
find . \( -path "./build" -o -path "./dist" -o -path "./.git" -o -path "./.venv" \) -prune -o -type d -name "*.egg-info" -exec rm -rf {} +
find . \( -path "./build" -o -path "./dist" -o -path "./.git" -o -path "./.venv" \) -prune -o -type d -name "__pycache__" -exec rm -rf {} +
find src/epcsaft -maxdepth 1 -type f \( -name "_core*.so" -o -name "_core*.pyd" \) -delete
