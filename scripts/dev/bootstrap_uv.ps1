$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent (Split-Path -Parent $ScriptDir)
Set-Location $RepoRoot

function Invoke-CheckedNative {
    param(
        [Parameter(Mandatory = $true)]
        [string]$FilePath,
        [string[]]$ArgumentList = @()
    )

    & $FilePath @ArgumentList
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

function Resolve-UvCommand {
    $command = Get-Command uv -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }
    $localUv = Join-Path $env:USERPROFILE ".local\bin\uv.exe"
    if (Test-Path -LiteralPath $localUv) {
        return $localUv
    }
    throw "uv is not installed or is not on PATH. Install uv first: https://docs.astral.sh/uv/getting-started/installation/"
}

$uv = Resolve-UvCommand

Invoke-CheckedNative $uv @("--version")
Invoke-CheckedNative $uv @("python", "pin", "3.13")
Invoke-CheckedNative $uv @("sync", "--no-install-project")

$defaultIpoptRoot = Join-Path $env:USERPROFILE "Documents\deps\ipopt-msvc"
$ipoptRoot = $env:EPCSAFT_IPOPT_ROOT
if ([string]::IsNullOrWhiteSpace($ipoptRoot)) {
    $ipoptRoot = $env:EPCSAFT_PEP517_IPOPT_ROOT
}
if ([string]::IsNullOrWhiteSpace($ipoptRoot) -and (Test-Path -LiteralPath $defaultIpoptRoot)) {
    $ipoptRoot = $defaultIpoptRoot
}
if (-not ([string]::IsNullOrWhiteSpace($ipoptRoot))) {
    $env:EPCSAFT_IPOPT_ROOT = (Resolve-Path -LiteralPath $ipoptRoot).Path
    $ipoptBin = Join-Path $env:EPCSAFT_IPOPT_ROOT "bin"
    if (Test-Path -LiteralPath $ipoptBin) {
        $env:EPCSAFT_RUNTIME_DLL_DIRS = $ipoptBin
        $env:PATH = "$ipoptBin;$env:PATH"
    }
}

Invoke-CheckedNative $uv @("run", "python", "scripts\dev\build_epcsaft.py")
Invoke-CheckedNative $uv @("run", "python", "scripts\dev\doctor.py", "--require-ipopt")
Invoke-CheckedNative $uv @("run", "python", "scripts\dev\validate_project.py", "quick")
