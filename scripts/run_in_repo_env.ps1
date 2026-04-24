param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$PythonArgs
)

$ErrorActionPreference = "Stop"

if (-not $PythonArgs -or $PythonArgs.Count -eq 0) {
    throw "Usage: scripts/run_in_repo_env.ps1 <python-args...>"
}

$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $repoRoot

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    throw "uv is not installed or is not on PATH. Install uv first, then run scripts/bootstrap_uv.ps1."
}

& uv run python @PythonArgs
exit $LASTEXITCODE
