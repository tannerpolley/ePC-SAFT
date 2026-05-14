$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent (Split-Path -Parent $ScriptDir)
Set-Location $RepoRoot

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    throw "uv is not installed or is not on PATH. Install uv first: https://docs.astral.sh/uv/getting-started/installation/"
}

uv --version
uv python pin 3.13
uv sync --no-install-project
uv run python scripts\dev\build_epcsaft.py
uv run python scripts\dev\doctor.py
uv run python scripts\dev\validate_project.py quick
