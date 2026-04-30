# Create a project-local Git worktree and register it as safe for Codex Git use.

[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[^\\/:*?"<>|]+$')]
    [string] $Name,

    [Parameter(Mandatory = $true)]
    [string] $Branch,

    [string] $Base = "HEAD",

    [string] $RepoRoot,

    [switch] $SkipSafeDirectory
)

$ErrorActionPreference = "Stop"

function Convert-ToGitPath {
    param([Parameter(Mandatory = $true)][string] $Path)
    return [System.IO.Path]::GetFullPath($Path).Replace("\", "/")
}

function Invoke-Git {
    param([Parameter(Mandatory = $true)][string[]] $Arguments)
    & git @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git $($Arguments -join ' ') failed with exit code $LASTEXITCODE"
    }
}

if ([string]::IsNullOrWhiteSpace($RepoRoot)) {
    $RepoRoot = (& git rev-parse --show-toplevel 2>$null)
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($RepoRoot)) {
        throw "Could not determine the repository root. Run this from a trusted repo checkout or pass -RepoRoot."
    }
}

$repoRootPath = Convert-ToGitPath $RepoRoot
if (-not (Test-Path -LiteralPath $repoRootPath -PathType Container)) {
    throw "Repository root does not exist: $repoRootPath"
}

$worktreesRoot = Join-Path $repoRootPath ".worktrees"
$worktreePath = Join-Path $worktreesRoot $Name
$worktreeGitPath = Convert-ToGitPath $worktreePath

& git -c "safe.directory=$repoRootPath" -C $repoRootPath check-ignore -q .worktrees/
if ($LASTEXITCODE -ne 0) {
    throw ".worktrees/ is not ignored. Add '.worktrees/' to .gitignore and commit that safety rule before creating project-local worktrees."
}

if (Test-Path -LiteralPath $worktreePath) {
    throw "Worktree path already exists: $worktreeGitPath"
}

if ($PSCmdlet.ShouldProcess($worktreeGitPath, "create Git worktree for branch $Branch from $Base")) {
    New-Item -ItemType Directory -Force -Path $worktreesRoot | Out-Null
    Invoke-Git @(
        "-c", "safe.directory=$repoRootPath",
        "-C", $repoRootPath,
        "worktree", "add",
        $worktreeGitPath,
        "-b", $Branch,
        $Base
    )
}

if (-not $SkipSafeDirectory) {
    $existingSafeDirectories = @(& git config --global --get-all safe.directory 2>$null)
    if ($existingSafeDirectories -notcontains $worktreeGitPath) {
        if ($PSCmdlet.ShouldProcess($worktreeGitPath, "add Git safe.directory entry")) {
            Invoke-Git @("config", "--global", "--add", "safe.directory", $worktreeGitPath)
        }
    }
}

Write-Host "worktree_path: $worktreeGitPath"
Write-Host "branch: $Branch"
Write-Host "safe_directory: $(-not $SkipSafeDirectory)"
