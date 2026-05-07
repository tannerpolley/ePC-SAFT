param(
    [string]$IpoptPrefix = "",
    [switch]$SkipSync
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $repoRoot

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    throw "uv is not installed or is not on PATH. Install uv first, then rerun this script."
}

if (-not $IpoptPrefix) {
    if ($env:EPCSAFT_IPOPT_PREFIX) {
        $IpoptPrefix = $env:EPCSAFT_IPOPT_PREFIX
    } elseif ($env:IPOPTWINDIR) {
        $IpoptPrefix = $env:IPOPTWINDIR
    } elseif ($env:CONDA_PREFIX -and (Test-Path -LiteralPath (Join-Path $env:CONDA_PREFIX "Library\lib\ipopt.lib"))) {
        $IpoptPrefix = Join-Path $env:CONDA_PREFIX "Library"
    } elseif (Test-Path -LiteralPath "C:\ProgramData\miniconda3\envs\epcsaft-cyipopt-test\Library\lib\ipopt.lib") {
        $IpoptPrefix = "C:\ProgramData\miniconda3\envs\epcsaft-cyipopt-test\Library"
    }
}

if (-not $IpoptPrefix) {
    throw "Could not locate an IPOPT prefix. Pass -IpoptPrefix <conda-env>\Library or set EPCSAFT_IPOPT_PREFIX."
}

$IpoptPrefix = (Resolve-Path -LiteralPath $IpoptPrefix).Path
$includeDir = Join-Path $IpoptPrefix "include\coin-or"
$libDir = Join-Path $IpoptPrefix "lib"
$binDir = Join-Path $IpoptPrefix "bin"

foreach ($path in @($includeDir, $libDir, $binDir)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required IPOPT directory is missing: $path"
    }
}
if (-not (Test-Path -LiteralPath (Join-Path $libDir "ipopt.lib"))) {
    throw "Required IPOPT import library is missing: $(Join-Path $libDir 'ipopt.lib')"
}
if (-not (Get-Command pkg-config -ErrorAction SilentlyContinue)) {
    $pkgConfig = Join-Path $binDir "pkg-config.exe"
    if (-not (Test-Path -LiteralPath $pkgConfig)) {
        throw "pkg-config is required. For conda-forge IPOPT, run: conda install -n <env> -c conda-forge pkg-config"
    }
}

$shimDir = Join-Path $repoRoot "build\cyipopt-uv-shim"
New-Item -ItemType Directory -Force -Path $shimDir | Out-Null
$prefixForPc = $IpoptPrefix.Replace("\", "/")
@"
prefix=$prefixForPc
exec_prefix=`${prefix}
libdir=`${exec_prefix}/lib
includedir=`${prefix}/include/coin-or

Name: Ipopt
Description: Interior Point Optimizer
URL: https://github.com/coin-or/Ipopt
Version: 3.14
Cflags: -I`${includedir}
Libs: -L`${libdir} -lipopt
"@ | Set-Content -LiteralPath (Join-Path $shimDir "ipopt.pc") -Encoding ASCII

$env:PKG_CONFIG_PATH = $shimDir
$env:EPCSAFT_IPOPT_DLL_DIR = $binDir
$env:PATH = "$binDir;$libDir;$shimDir;$env:PATH"

if (-not $SkipSync) {
    & uv sync --group ipopt
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

& uv run python -c "import epcsaft; info=epcsaft.capabilities()['optimizers']['ipopt']; print(info)"
exit $LASTEXITCODE
