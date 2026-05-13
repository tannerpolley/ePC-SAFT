param(
    [int]$PollIntervalSeconds = 120,
    [int]$MaxWaitMinutes = 480,
    [switch]$Once
)

$ErrorActionPreference = 'Stop'
$GoalRoot = Split-Path -Parent $PSScriptRoot
$GatePath = Join-Path $PSScriptRoot 'dependency_gate.yaml'
$StatusDir = Join-Path $PSScriptRoot '.goalbuddy-board'
$StatusPath = Join-Path $StatusDir 'watcher_status.json'
New-Item -ItemType Directory -Force -Path $StatusDir | Out-Null

function Write-WatcherStatus {
    param([string]$Status, [string]$Message)
    [pscustomobject]@{
        status = $Status
        message = $Message
        checked_at = (Get-Date).ToString('o')
        gate_path = $GatePath
        watcher_mode = 'bounded'
        auto_start_after_gate = $true
        poll_interval_seconds = $PollIntervalSeconds
        max_wait_minutes = $MaxWaitMinutes
    } | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

if (-not (Test-Path -LiteralPath $GatePath)) {
    Write-WatcherStatus -Status 'BLOCKED_DEPENDENCY_OR_REBASE' -Message 'Missing notes/dependency_gate.yaml.'
    exit 2
}

Write-WatcherStatus -Status 'PREPARED_WAITING' -Message 'Generic bounded watcher stub is installed; inspect dependency_gate.yaml and continue when the task-specific gate is satisfied.'
if ($Once) { exit 2 }
$Deadline = (Get-Date).AddMinutes($MaxWaitMinutes)
while ((Get-Date) -lt $Deadline) {
    Start-Sleep -Seconds $PollIntervalSeconds
    Write-WatcherStatus -Status 'PREPARED_WAITING' -Message 'Waiting for dependency gate to be satisfied.'
}
Write-WatcherStatus -Status 'PREPARED_WAITING' -Message "Timed out after $MaxWaitMinutes minutes."
exit 3
