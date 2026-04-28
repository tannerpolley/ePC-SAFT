$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent $ScriptDir
Set-Location $RepoRoot

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    throw "uv is not installed or is not on PATH. Install uv first: https://docs.astral.sh/uv/getting-started/installation/"
}

uv --version
uv python pin 3.12
uv sync --no-install-project
uv run python scripts\build_epcsaft.py
uv run python scripts\codex_doctor.py
uv run python run_pytest.py --confidence -q
