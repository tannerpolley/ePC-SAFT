$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent $ScriptDir
Set-Location $RepoRoot

Remove-Item -Recurse -Force build, dist, .pytest_cache, .ruff_cache, .mypy_cache -ErrorAction SilentlyContinue
Get-ChildItem -Recurse -Directory -Filter "*.egg-info" | Remove-Item -Recurse -Force
Get-ChildItem -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force
Get-ChildItem -Path src\epcsaft -File -Filter "_core*.pyd" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem -Path src\epcsaft -File -Filter "_core*.so" -ErrorAction SilentlyContinue | Remove-Item -Force
