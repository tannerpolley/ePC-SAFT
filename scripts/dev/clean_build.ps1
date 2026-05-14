$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent (Split-Path -Parent $ScriptDir)
Set-Location $RepoRoot

Write-Warning "REPAIR-ONLY: this removes build/cache/native artifacts. Use normal build and confidence tests for routine validation."

function Remove-PathIfPresent {
    param([Parameter(Mandatory = $true)][string]$Path)

    if (Test-Path -LiteralPath $Path) {
        try {
            Remove-Item -LiteralPath $Path -Recurse -Force -ErrorAction Stop
        }
        catch {
            Write-Warning "Could not remove ${Path}: $($_.Exception.Message)"
        }
    }
}

foreach ($path in @("build", "dist", ".pytest_cache", ".ruff_cache", ".mypy_cache")) {
    Remove-PathIfPresent $path
}

$excludedTopLevelDirs = @("build", "dist", ".git", ".venv")
$scanRoots = Get-ChildItem -LiteralPath . -Force -Directory -ErrorAction SilentlyContinue |
    Where-Object { $excludedTopLevelDirs -notcontains $_.Name }

foreach ($root in $scanRoots) {
    $dirs = @($root) + @(Get-ChildItem -LiteralPath $root.FullName -Recurse -Directory -Force -ErrorAction SilentlyContinue)
    $dirs | Where-Object { $_.Name -like "*.egg-info" -or $_.Name -eq "__pycache__" } |
        ForEach-Object { Remove-PathIfPresent $_.FullName }
}

Get-ChildItem -Path src\epcsaft -File -Filter "_core*.pyd" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem -Path src\epcsaft -File -Filter "_core*.so" -ErrorAction SilentlyContinue | Remove-Item -Force
