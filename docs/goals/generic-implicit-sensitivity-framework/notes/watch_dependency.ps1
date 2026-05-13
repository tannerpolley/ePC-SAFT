param(
    [int]$PollIntervalSeconds = 120,
    [int]$MaxWaitMinutes = 480,
    [switch]$Once
)

$ErrorActionPreference = "Stop"

$ExpectedBranch = "codex/generic-implicit-sensitivity-framework"
$Repo = "tannerpolley/ePC-SAFT"
$TaskBBranch = "codex/cppad-explicit-parameter-derivatives"
$StartedAt = Get-Date
$Deadline = $StartedAt.AddMinutes($MaxWaitMinutes)

function Write-WatcherStatus {
    param(
        [string]$Status,
        [string]$Message,
        [object]$TaskBPr = $null
    )

    $payload = [ordered]@{
        status = $Status
        message = $Message
        checked_at = (Get-Date).ToString("o")
        expected_branch = $ExpectedBranch
        task_b_branch = $TaskBBranch
        poll_interval_seconds = $PollIntervalSeconds
        max_wait_minutes = $MaxWaitMinutes
        started_at = $StartedAt.ToString("o")
        deadline = $Deadline.ToString("o")
        task_b_pr = $TaskBPr
    }

    $statusPath = Join-Path $PSScriptRoot ".goalbuddy-board\watcher_status.json"
    if (-not (Test-Path -LiteralPath (Split-Path -Parent $statusPath))) {
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $statusPath) | Out-Null
    }
    $payload | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding UTF8
    Write-Output "$Status $Message"
}

function Get-TaskBMergedPr {
    $json = gh pr list --repo $Repo --state all --head $TaskBBranch --json number,title,state,mergedAt,mergeCommit,url,headRefName,baseRefName --limit 20
    if ($LASTEXITCODE -ne 0) {
        throw "gh pr list failed for $TaskBBranch"
    }

    $prs = @($json | ConvertFrom-Json)
    $merged = @($prs | Where-Object { $_.state -eq "MERGED" -and $_.mergedAt })
    if ($merged.Count -eq 0) {
        return $null
    }
    return $merged | Sort-Object mergedAt -Descending | Select-Object -First 1
}

function Confirm-And-UpdateBranch {
    git fetch origin --prune
    if ($LASTEXITCODE -ne 0) {
        throw "git fetch origin --prune failed"
    }

    $branch = (git branch --show-current).Trim()
    if ($branch -ne $ExpectedBranch) {
        throw "BRANCH_BOOTSTRAP_FAILED expected=$ExpectedBranch actual=$branch"
    }

    $status = git status --short
    $nonGoalChanges = @($status | Where-Object { $_ -and ($_ -notmatch "docs/goals/generic-implicit-sensitivity-framework/") })
    if ($nonGoalChanges.Count -gt 0) {
        throw "BLOCKED_DIRTY_WORKTREE non-goal changes present: $($nonGoalChanges -join '; ')"
    }

    git rebase origin/main
    if ($LASTEXITCODE -ne 0) {
        throw "BLOCKED_DEPENDENCY_OR_REBASE rebase_failed"
    }
}

do {
    $mergedPr = Get-TaskBMergedPr
    if ($null -ne $mergedPr) {
        Confirm-And-UpdateBranch
        Write-WatcherStatus -Status "GATE_PASSED" -Message "Task B merged; branch updated on origin/main. Continue GoalBuddy implementation without asking." -TaskBPr $mergedPr
        exit 0
    }

    Write-WatcherStatus -Status "PREPARED_WAITING" -Message "Task B is not merged yet."
    if ($Once) {
        exit 2
    }

    if ((Get-Date) -ge $Deadline) {
        Write-WatcherStatus -Status "PREPARED_WAITING" -Message "Timed out waiting for Task B merge after $MaxWaitMinutes minutes."
        exit 3
    }

    Start-Sleep -Seconds $PollIntervalSeconds
} while ($true)
