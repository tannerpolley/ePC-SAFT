param(
    [ValidateSet("Setup", "Build", "Doctor")]
    [string]$Step = "Setup"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repoRoot = (git rev-parse --show-toplevel).Trim()
Set-Location $repoRoot

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    throw "uv is required for this repo workflow. Install uv, then rerun setup."
}

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

function Set-NativeIpoptEnvironment {
    $defaultIpoptRoot = "C:\ProgramData\miniconda3\envs\ePC-SAFT\Library"
    $candidate = $env:EPCSAFT_IPOPT_ROOT

    if ([string]::IsNullOrWhiteSpace($candidate)) {
        $candidate = $env:EPCSAFT_PEP517_IPOPT_ROOT
    }
    if ([string]::IsNullOrWhiteSpace($candidate) -and (Test-Path -LiteralPath $defaultIpoptRoot)) {
        $candidate = $defaultIpoptRoot
    }
    if ([string]::IsNullOrWhiteSpace($candidate)) {
        throw "Native Ipopt is required for this ePC-SAFT setup. Set EPCSAFT_IPOPT_ROOT or install Ipopt at $defaultIpoptRoot."
    }

    $ipoptRoot = (Resolve-Path -LiteralPath $candidate).Path
    $ipoptBin = Join-Path $ipoptRoot "bin"
    $ipoptLib = Join-Path $ipoptRoot "lib"

    if (-not (Test-Path -LiteralPath $ipoptBin)) {
        throw "Native Ipopt runtime DLL directory was not found: $ipoptBin"
    }
    if (-not (Test-Path -LiteralPath $ipoptLib)) {
        throw "Native Ipopt library directory was not found: $ipoptLib"
    }

    $env:EPCSAFT_IPOPT_ROOT = $ipoptRoot
    $env:EPCSAFT_RUNTIME_DLL_DIRS = $ipoptBin

    $pathEntries = $env:PATH -split ";"
    if ($pathEntries -notcontains $ipoptBin) {
        $env:PATH = "$ipoptBin;$env:PATH"
    }

    Write-Host "Using native Ipopt root: $ipoptRoot"
}

function Invoke-NativeBuild {
    Set-NativeIpoptEnvironment
    Invoke-CheckedNative uv @("run", "python", "scripts/dev/build_epcsaft.py")
}

function Invoke-NativeDoctor {
    Set-NativeIpoptEnvironment
    Invoke-CheckedNative uv @("run", "python", "scripts/dev/doctor.py", "--require-ipopt")
}

switch ($Step) {
    "Setup" {
        Invoke-CheckedNative uv @("sync", "--no-install-project")
        Invoke-NativeBuild
        Invoke-NativeDoctor
    }
    "Build" {
        Invoke-NativeBuild
    }
    "Doctor" {
        Invoke-NativeDoctor
    }
}
