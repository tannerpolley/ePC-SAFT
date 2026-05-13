# Bounded dependency watcher template
# Save per-agent copy under docs/goals/<slug>/watch_dependency.ps1 if needed.

param(
    [string]$Branch,
    [string[]]$IssueNumbers = @(),
    [string[]]$PrNumbers = @(),
    [int]$PollIntervalSeconds = 120,
    [int]$MaxWaitMinutes = 480
)

$deadline = (Get-Date).AddMinutes($MaxWaitMinutes)

function Stop-WithStatus($Message) {
    Write-Host $Message
    exit 1
}

while ((Get-Date) -lt $deadline) {
    git fetch origin --prune

    $currentBranch = (git branch --show-current).Trim()
    if ($currentBranch -ne $Branch) {
        Stop-WithStatus "BLOCKED_BRANCH_MISMATCH expected=$Branch actual=$currentBranch"
    }

    $allIssuesClosed = $true
    foreach ($issue in $IssueNumbers) {
        $state = gh issue view $issue --json state --jq ".state"
        if ($state -ne "CLOSED") {
            $allIssuesClosed = $false
        }
    }

    $allPrsMerged = $true
    foreach ($pr in $PrNumbers) {
        $mergedAt = gh pr view $pr --json mergedAt --jq ".mergedAt"
        if ([string]::IsNullOrWhiteSpace($mergedAt)) {
            $allPrsMerged = $false
        }
    }

    if ($allIssuesClosed -and $allPrsMerged) {
        git rebase origin/main
        if ($LASTEXITCODE -ne 0) {
            git rebase --abort
            Stop-WithStatus "BLOCKED_REBASE_CONFLICT"
        }
        Write-Host "GATE_PASS"
        exit 0
    }

    Write-Host "PREPARED_WAITING dependencies not met. Sleeping $PollIntervalSeconds seconds."
    Start-Sleep -Seconds $PollIntervalSeconds
}

Write-Host "PREPARED_WAITING timeout reached without dependency gate pass."
exit 2
