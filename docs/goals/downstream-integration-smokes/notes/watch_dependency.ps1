# Bounded dependency watcher template
# Per-agent copy for downstream-integration-smokes. Defaults are filled for this worktree branch.

param(
    [string]$Branch = "codex/downstream-integration-smokes",
    [string[]]$IssueNumbers = @("89", "90", "92", "93", "94"),
    [string[]]$PrNumbers = @(),
    [int]$PollIntervalSeconds = 120,
    [int]$MaxWaitMinutes = 480,
    [ValidateSet("rebase", "ff-only")]
    [string]$UpdateMode = "ff-only",
    [switch]$Continuous
)

$deadline = (Get-Date).AddMinutes($MaxWaitMinutes)

function Stop-WithStatus($Message) {
    Write-Host $Message
    exit 1
}

function Require-Tool($Name) {
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        Stop-WithStatus "BLOCKED_MISSING_TOOL $Name"
    }
}

function Test-Dependencies {
    if (($IssueNumbers.Count -gt 0) -or ($PrNumbers.Count -gt 0)) {
        Require-Tool "gh"
    }

    foreach ($issue in $IssueNumbers) {
        $state = gh issue view $issue --json state --jq ".state"
        if ($LASTEXITCODE -ne 0) {
            Stop-WithStatus "BLOCKED_DEPENDENCY issue=$issue state_lookup_failed"
        }
        if ($state -ne "CLOSED") {
            return $false
        }
    }

    foreach ($pr in $PrNumbers) {
        $mergedAt = gh pr view $pr --json mergedAt --jq ".mergedAt"
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
        Write-Host "GATE_PASS"
        exit 0
    }

    Write-Host "PREPARED_WAITING dependencies not met."
    return $false
}

if (-not $Continuous) {
    Invoke-GateCheck | Out-Null
    exit 2
}

while ((Get-Date) -lt $deadline) {
    Invoke-GateCheck | Out-Null
    Write-Host "PREPARED_WAITING sleeping $PollIntervalSeconds seconds."
    Start-Sleep -Seconds $PollIntervalSeconds
}

Write-Host "PREPARED_WAITING timeout reached without dependency gate pass."
exit 2
