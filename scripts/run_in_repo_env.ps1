param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$PythonArgs
)

$ErrorActionPreference = 'Stop'

if (-not $PythonArgs -or $PythonArgs.Count -eq 0) {
    throw 'Usage: scripts/run_in_repo_env.ps1 <python-args...>'
}

$repoRoot = Split-Path -Parent $PSScriptRoot
$envName = Split-Path -Leaf $repoRoot

$pythonCandidates = @(
    "C:\ProgramData\Miniconda3\envs\$envName\python.exe",
    "C:\ProgramData\miniconda3\envs\$envName\python.exe",
    ($env:USERPROFILE + "\miniconda3\envs\$envName\python.exe"),
    ($env:USERPROFILE + "\anaconda3\envs\$envName\python.exe")
)

foreach ($candidate in $pythonCandidates) {
    if (Test-Path $candidate) {
        & $candidate @PythonArgs
        exit $LASTEXITCODE
    }
}

if (Get-Command conda -ErrorAction SilentlyContinue) {
    & conda run -n $envName python @PythonArgs
    exit $LASTEXITCODE
}

throw "Could not locate python.exe for conda environment '$envName'."
