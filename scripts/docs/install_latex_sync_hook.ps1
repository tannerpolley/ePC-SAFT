$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = (git rev-parse --show-toplevel).Trim()
if (-not $repoRoot) {
    throw "Unable to resolve repo root."
}

$hookPath = (git rev-parse --git-path hooks/post-commit).Trim()
if (-not $hookPath) {
    throw "Unable to resolve post-commit hook path."
}

$hookDir = Split-Path -Parent $hookPath
if (-not (Test-Path -LiteralPath $hookDir)) {
    New-Item -ItemType Directory -Path $hookDir | Out-Null
}

$syncScript = (Join-Path $repoRoot "scripts/docs/sync_latex_mirror.ps1").Replace("\", "/")
$hookContent = @"
#!/bin/sh
set -eu

if ! git diff-tree --no-commit-id --name-only -r HEAD | grep -q '^docs/latex/'; then
    exit 0
fi

if command -v pwsh >/dev/null 2>&1; then
    pwsh -NoProfile -ExecutionPolicy Bypass -File "$syncScript"
else
    powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$syncScript"
fi
"@

Set-Content -LiteralPath $hookPath -Value $hookContent -Encoding ASCII
Write-Output "Installed LaTeX mirror post-commit hook: $hookPath"
