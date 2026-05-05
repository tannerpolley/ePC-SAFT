param(
    [string]$MirrorRoot = "C:\Users\Tanner\Documents\git\LaTeX-Projects\ePC-SAFT",
    [string]$RemoteUrl = "https://git@git.overleaf.com/686d837e6bb7fceea660fb5b",
    [string]$Branch = "master",
    [string]$Message = "",
    [switch]$NoCommit,
    [switch]$NoPush
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Invoke-Git {
    & git @args
    if ($LASTEXITCODE -ne 0) {
        throw "git $($args -join ' ') failed with exit code $LASTEXITCODE"
    }
}

$repoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$sourceRoot = Join-Path $repoRoot "docs\latex"
if (-not (Test-Path -LiteralPath $sourceRoot)) {
    throw "Missing LaTeX source directory: $sourceRoot"
}
if (-not (Test-Path -LiteralPath (Join-Path $sourceRoot "equations.tex"))) {
    throw "Missing LaTeX source file: $(Join-Path $sourceRoot 'equations.tex')"
}
if (-not (Test-Path -LiteralPath (Join-Path $MirrorRoot ".git"))) {
    throw "Missing LaTeX mirror Git repository: $MirrorRoot. Run scripts\setup_latex_mirror.ps1 first."
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

& robocopy $sourceRoot $MirrorRoot /MIR /XD .git /XF .git /R:2 /W:1 /NFL /NDL /NP
$robocopyExit = $LASTEXITCODE
if ($robocopyExit -gt 7) {
    throw "robocopy failed with exit code $robocopyExit"
}

Invoke-Git -C $MirrorRoot add -A
$status = & git -C $MirrorRoot status --porcelain
if ($LASTEXITCODE -ne 0) {
    throw "Unable to read mirror status at $MirrorRoot"
}
if (-not $status) {
    Write-Output "LaTeX mirror already matches docs\latex."
    exit 0
}

if ($NoCommit) {
    Write-Output "LaTeX mirror has staged changes. Commit skipped because -NoCommit was supplied."
    exit 0
}

if (-not $Message) {
    $shortSha = (& git -C $repoRoot rev-parse --short HEAD).Trim()
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to read ePC-SAFT commit SHA."
    }
    $Message = "Sync LaTeX mirror from ePC-SAFT $shortSha"
}

Invoke-Git -C $MirrorRoot commit -m $Message

if ($NoPush) {
    Write-Output "LaTeX mirror committed. Push skipped because -NoPush was supplied."
    exit 0
}

Invoke-Git -C $MirrorRoot push origin $Branch
Write-Output "LaTeX mirror synced, committed, and pushed: $MirrorRoot"
