param(
    [string]$Branch = "codex/generic-reactive-lle",
    [string[]]$IssueNumbers = @("92", "87", "89"),
    [string[]]$PrNumbers = @("110", "99", "105"),
    [int]$PollIntervalSeconds = 360,
    [int]$MaxWaitMinutes = 480,
    [ValidateSet("rebase", "ff-only")]
    [string]$UpdateMode = "ff-only",
    [switch]$Watch,
    [switch]$Once
)

$ErrorActionPreference = "Stop"
$StatusDir = Join-Path $PSScriptRoot ".goalbuddy-board"
$StatusPath = Join-Path $StatusDir "watcher_status.json"
New-Item -ItemType Directory -Force -Path $StatusDir | Out-Null

function Write-WatcherStatus {
    param([string]$Status, [string]$Message)

    [pscustomobject]@{
        status = $Status
        message = $Message
        checked_at = (Get-Date).ToString("o")
        branch = $Branch
        issue_numbers = $IssueNumbers
        pr_numbers = $PrNumbers
        watcher_mode = "bounded"
        auto_start_after_gate = $true
        poll_interval_seconds = $PollIntervalSeconds
        max_wait_minutes = $MaxWaitMinutes
        update_mode = $UpdateMode
        task_key = "J"
    } | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

function Stop-WithStatus {
    param([string]$Message)

    Write-WatcherStatus -Status ($Message.Split(" ")[0]) -Message $Message
    Write-Host $Message
    exit 1
}

function Require-Tool {
    param([string]$Name)

    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        Stop-WithStatus "BLOCKED_MISSING_TOOL $Name"
    }
}

function Test-Dependencies {
    if (($IssueNumbers.Count -gt 0) -or ($PrNumbers.Count -gt 0)) {
        Require-Tool "gh"
    }

    foreach ($issue in $IssueNumbers) {
        $state = gh issue view $issue --repo tannerpolley/ePC-SAFT --json state --jq ".state"
        if ($LASTEXITCODE -ne 0) {
            Stop-WithStatus "BLOCKED_DEPENDENCY issue=$issue state_lookup_failed"
        }
        if ($state -ne "CLOSED") {
            return $false
        }
    }

    foreach ($pr in $PrNumbers) {
        $mergedAt = gh pr view $pr --repo tannerpolley/ePC-SAFT --json mergedAt --jq ".mergedAt"
        if ($LASTEXITCODE -ne 0) {
            Stop-WithStatus "BLOCKED_DEPENDENCY pr=$pr state_lookup_failed"
        }
        if ([string]::IsNullOrWhiteSpace($mergedAt)) {
            return $false
        }
    }

    return $true
}

function Update-BranchFromMain {
    if ($UpdateMode -eq "ff-only") {
        git merge --ff-only origin/main
        if ($LASTEXITCODE -ne 0) {
            Stop-WithStatus "BLOCKED_REBASE_CONFLICT update_mode=ff-only"
        }
    } else {
        git rebase origin/main
        if ($LASTEXITCODE -ne 0) {
            git rebase --abort 2>$null
            Stop-WithStatus "BLOCKED_REBASE_CONFLICT update_mode=rebase"
        }
    }
}

function Invoke-GateCheck {
    git fetch origin --prune
    if ($LASTEXITCODE -ne 0) {
        Stop-WithStatus "BLOCKED_FETCH_ORIGIN"
    }

    $currentBranch = (git branch --show-current).Trim()
    if ($currentBranch -ne $Branch) {
        Stop-WithStatus "BLOCKED_BRANCH_MISMATCH expected=$Branch actual=$currentBranch"
    }

    if (Test-Dependencies) {
        Update-BranchFromMain
        Write-WatcherStatus -Status "GATE_PASS" -Message "Dependencies are closed and merged; branch updated from origin/main with $UpdateMode."
        Write-Host "GATE_PASS"
        exit 0
    }

    Write-WatcherStatus -Status "PREPARED_WAITING" -Message "Dependencies are not all closed and merged."
    Write-Host "PREPARED_WAITING dependencies not met."
    return $false
}

$deadline = (Get-Date).AddMinutes($MaxWaitMinutes)

do {
    $gatePassed = Invoke-GateCheck
    if ((-not $Watch) -or $Once) {
        exit 2
    }

    if ((Get-Date) -ge $deadline) {
        Write-WatcherStatus -Status "PREPARED_WAITING" -Message "Timeout reached without dependency gate pass."
        Write-Host "PREPARED_WAITING timeout reached without dependency gate pass."
        exit 2
    }

    Write-Host "PREPARED_WAITING dependencies not met. Sleeping $PollIntervalSeconds seconds."
    Start-Sleep -Seconds $PollIntervalSeconds
} while (-not $gatePassed)
