param(
    [string]$MirrorRoot = "C:\Users\Tanner\Documents\git\LaTeX-Projects\ePC-SAFT-LaTeX",
    [string]$RemoteUrl = "https://git@git.overleaf.com/686d837e6bb7fceea660fb5b",
    [string]$Branch = "master"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Invoke-Git {
    & git @args
    if ($LASTEXITCODE -ne 0) {
        throw "git $($args -join ' ') failed with exit code $LASTEXITCODE"
    }
}

$parent = Split-Path -Parent $MirrorRoot
if (-not (Test-Path -LiteralPath $parent)) {
    New-Item -ItemType Directory -Path $parent | Out-Null
}

if (-not (Test-Path -LiteralPath $MirrorRoot)) {
    Invoke-Git clone --branch $Branch $RemoteUrl $MirrorRoot
}
elseif (-not (Test-Path -LiteralPath (Join-Path $MirrorRoot ".git"))) {
    throw "Mirror path exists but is not a Git repository: $MirrorRoot"
}

$insideWorkTree = (& git -C $MirrorRoot rev-parse --is-inside-work-tree).Trim()
if ($LASTEXITCODE -ne 0 -or $insideWorkTree -ne "true") {
    throw "Mirror path is not a valid Git work tree: $MirrorRoot"
}

$actualRemote = (& git -C $MirrorRoot remote get-url origin).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Unable to read mirror origin remote at $MirrorRoot"
}
if ($actualRemote -ne $RemoteUrl) {
    throw "Mirror origin remote mismatch. Expected '$RemoteUrl' but found '$actualRemote'."
}

$actualBranch = (& git -C $MirrorRoot branch --show-current).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Unable to read mirror branch at $MirrorRoot"
}
if ($actualBranch -ne $Branch) {
    throw "Mirror branch mismatch. Expected '$Branch' but found '$actualBranch'."
}

Write-Output "LaTeX mirror is ready: $MirrorRoot"
