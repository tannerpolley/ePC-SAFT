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

function Get-GateStatus {
    if (-not (Test-Path -LiteralPath $GatePath)) {
        return $null
    }

    $text = Get-Content -LiteralPath $GatePath -Raw
    $match = [regex]::Match($text, '(?m)^status:\s*([A-Z_]+)\s*$')
    if (-not $match.Success) {
        return $null
    }

    return $match.Groups[1].Value
}

function Write-WatcherStatus {
    param([string]$Status, [string]$Message, [string]$ObservedGateStatus)

    [pscustomobject]@{
        status = $Status
        message = $Message
        checked_at = (Get-Date).ToString('o')
        gate_path = $GatePath
        observed_gate_status = $ObservedGateStatus
        watcher_mode = 'bounded'
        auto_start_after_gate = $true
        poll_interval_seconds = $PollIntervalSeconds
        max_wait_minutes = $MaxWaitMinutes
        task_key = 'J'
    } | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

if (-not (Test-Path -LiteralPath $GatePath)) {
    Write-WatcherStatus -Status 'BLOCKED_DEPENDENCY_OR_REBASE' -Message 'Missing notes/dependency_gate.yaml.' -ObservedGateStatus 'missing'
    exit 2
}

$observedGateStatus = Get-GateStatus
if ($observedGateStatus -eq 'GATE_PASSED') {
    Write-WatcherStatus -Status 'GATE_PASSED' -Message 'Dependency gate is satisfied. Continue implementation without asking.' -ObservedGateStatus $observedGateStatus
    exit 0
}

if (-not $observedGateStatus) {
    Write-WatcherStatus -Status 'PREPARED_WAITING' -Message 'Dependency gate status is missing or unreadable.' -ObservedGateStatus 'unreadable'
} else {
    Write-WatcherStatus -Status 'PREPARED_WAITING' -Message "Dependency gate is $observedGateStatus." -ObservedGateStatus $observedGateStatus
}

if ($Once) { exit 2 }

$Deadline = (Get-Date).AddMinutes($MaxWaitMinutes)
while ((Get-Date) -lt $Deadline) {
    Start-Sleep -Seconds $PollIntervalSeconds
    $observedGateStatus = Get-GateStatus
    if ($observedGateStatus -eq 'GATE_PASSED') {
        Write-WatcherStatus -Status 'GATE_PASSED' -Message 'Dependency gate is satisfied. Continue implementation without asking.' -ObservedGateStatus $observedGateStatus
        exit 0
    }
    if (-not $observedGateStatus) {
        Write-WatcherStatus -Status 'PREPARED_WAITING' -Message 'Dependency gate status is missing or unreadable.' -ObservedGateStatus 'unreadable'
    } else {
        Write-WatcherStatus -Status 'PREPARED_WAITING' -Message "Waiting for dependency gate to reach GATE_PASSED; current status is $observedGateStatus." -ObservedGateStatus $observedGateStatus
    }
}

Write-WatcherStatus -Status 'PREPARED_WAITING' -Message "Timed out after $MaxWaitMinutes minutes." -ObservedGateStatus (Get-GateStatus)
exit 3
