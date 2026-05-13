param(
    [string]$Branch = 'codex/generic-electrolyte-lle',
    [string[]]$IssueNumbers = @('91', '86'),
    [string[]]$PrNumbers = @(),
    [int]$PollIntervalSeconds = 120,
    [int]$MaxWaitMinutes = 480,
    [ValidateSet('rebase', 'ff-only')]
    [string]$UpdateMode = 'ff-only',
    [switch]$Once
)

$ErrorActionPreference = 'Stop'
$Repo = 'tannerpolley/ePC-SAFT'
$StartedAt = Get-Date
$Deadline = $StartedAt.AddMinutes($MaxWaitMinutes)

function Write-WatcherStatus {
    param(
        [string]$Status,
        [string]$Message,
        [object]$DependencyState = $null
    )

    $payload = [ordered]@{
        status = $Status
        message = $Message
        checked_at = (Get-Date).ToString('o')
        branch = $Branch
        poll_interval_seconds = $PollIntervalSeconds
        max_wait_minutes = $MaxWaitMinutes
        update_mode = $UpdateMode
        started_at = $StartedAt.ToString('o')
        deadline = $Deadline.ToString('o')
        dependency_state = $DependencyState
        issue_numbers = $IssueNumbers
        pr_numbers = $PrNumbers
    }

    $statusPath = Join-Path $PSScriptRoot '.goalbuddy-board\watcher_status.json'
    $statusDir = Split-Path -Parent $statusPath
    if (-not (Test-Path -LiteralPath $statusDir)) {
        New-Item -ItemType Directory -Force -Path $statusDir | Out-Null
    }

    $payload | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding UTF8
    Write-Host "$Status $Message"
}

function Get-IssueState {
    param([string]$IssueNumber)

    $json = gh issue view $IssueNumber --repo $Repo --json state,url,title
    if ($LASTEXITCODE -ne 0) {
        throw "gh issue view failed for issue $IssueNumber"
    }

    return $json | ConvertFrom-Json
}

function Confirm-Branch {
    git fetch origin --prune
    if ($LASTEXITCODE -ne 0) {
        throw 'BLOCKED_FETCH_ORIGIN'
    }

    $currentBranch = (git branch --show-current).Trim()
    if ($currentBranch -ne $Branch) {
        throw "BLOCKED_BRANCH_MISMATCH expected=$Branch actual=$currentBranch"
    }
}

function Update-BranchFromMain {
    if ($UpdateMode -eq 'ff-only') {
        git merge --ff-only origin/main
        if ($LASTEXITCODE -ne 0) {
            throw 'BLOCKED_REBASE_CONFLICT update_mode=ff-only'
        }
        return
    }

    git rebase origin/main
    if ($LASTEXITCODE -ne 0) {
        git rebase --abort 2>$null
        throw 'BLOCKED_REBASE_CONFLICT update_mode=rebase'
    }
}

function Test-Dependencies {
    $dependencyState = [ordered]@{
        issues = @()
        prs = @()
    }

    foreach ($issueNumber in $IssueNumbers) {
        $issue = Get-IssueState -IssueNumber $issueNumber
        $issueState = [string]$issue.state
        if ([string]::IsNullOrWhiteSpace($issueState)) {
            $issueState = ''
        } else {
            $issueState = $issueState.ToUpperInvariant()
        }
        $dependencyState.issues += [ordered]@{
            number = [int]$issueNumber
            state = $issueState
            url = $issue.url
            title = $issue.title
        }

        if ($issueState -ne 'CLOSED') {
            return [pscustomobject]@{
                Passed = $false
                State = $dependencyState
            }
        }
    }

    foreach ($prNumber in $PrNumbers) {
        $mergedAt = gh pr view $prNumber --repo $Repo --json mergedAt --jq '.mergedAt'
        if ($LASTEXITCODE -ne 0) {
            throw "gh pr view failed for pr $prNumber"
        }

        $dependencyState.prs += [ordered]@{
            number = [int]$prNumber
            mergedAt = $mergedAt
        }

        if ([string]::IsNullOrWhiteSpace($mergedAt)) {
            return [pscustomobject]@{
                Passed = $false
                State = $dependencyState
            }
        }
    }

    return [pscustomobject]@{
        Passed = $true
        State = $dependencyState
    }
}

function Test-DependencyGate {
    Confirm-Branch
    $result = Test-Dependencies
    if ($result.Passed) {
        Update-BranchFromMain
        Write-WatcherStatus -Status 'GATE_PASS' -Message 'Dependencies are satisfied; branch updated from origin/main.' -DependencyState $result.State
        return $true
    }

    Write-WatcherStatus -Status 'PREPARED_WAITING' -Message 'Dependencies are not satisfied yet.' -DependencyState $result.State
    return $false
}

do {
    if (Test-DependencyGate) {
        exit 0
    }

    if ($Once) {
        exit 2
    }

    if ((Get-Date) -ge $Deadline) {
        Write-WatcherStatus -Status 'PREPARED_WAITING' -Message "Timed out after $MaxWaitMinutes minutes." -DependencyState $null
        exit 3
    }

    Start-Sleep -Seconds $PollIntervalSeconds
} while ($true)
