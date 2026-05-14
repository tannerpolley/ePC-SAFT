$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = (git rev-parse --show-toplevel).Trim()
if (-not $repoRoot) {
    throw "Unable to resolve repo root."
}

$latexDir = Join-Path $repoRoot "docs\latex"
$sourceTex = Join-Path $latexDir "equations.tex"
$builtPdf = Join-Path $latexDir "out\equations.pdf"
$trackedPdf = Join-Path $latexDir "equations.pdf"

if (-not (Test-Path $sourceTex)) {
    throw "Missing source file: $sourceTex"
}

Push-Location $latexDir
try {
    & latexmk -pdf "equations.tex"
    if (-not (Test-Path $builtPdf)) {
        throw "Expected built PDF not found: $builtPdf"
    }
    $shouldCopy = $true
    if (Test-Path $trackedPdf) {
        $builtHash = (Get-FileHash -LiteralPath $builtPdf -Algorithm SHA256).Hash
        $trackedHash = (Get-FileHash -LiteralPath $trackedPdf -Algorithm SHA256).Hash
        $shouldCopy = $builtHash -ne $trackedHash
    }
    if ($shouldCopy) {
        Copy-Item -Force $builtPdf $trackedPdf
        Write-Output "Updated tracked equations PDF: $trackedPdf"
    }
    else {
        Write-Output "Tracked equations PDF already up to date: $trackedPdf"
    }
}
finally {
    Pop-Location
}
