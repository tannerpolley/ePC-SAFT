param(
    [string]$Repo = "tannerpolley/ePC-SAFT",
    [int]$Issue = 67,
    [string]$RootCheckout = "C:\Users\Tanner\Documents\git\ePC-SAFT",
    [string]$Worktree = "C:\Users\Tanner\.codex\worktrees\7c19\ePC-SAFT",
    [string]$Branch = "codex/issue-68-ceres-cppad-ci-capabilities",
    [int]$IntervalSeconds = 300
)

$ErrorActionPreference = "Stop"
$Gh = "C:\Program Files\GitHub CLI\gh.exe"

function Invoke-GhJson {
    param([string[]]$Arguments)
    & $Gh @Arguments | ConvertFrom-Json
}

function Get-IssueState {
    Invoke-GhJson @("issue", "view", "$Issue", "--repo", $Repo, "--json", "number,title,state,stateReason,closedAt,url")
}

function Get-ClosingMergedPr {
    $owner, $name = $Repo -split "/", 2
    $query = @'
query($owner:String!,$name:String!,$number:Int!){
  repository(owner:$owner,name:$name){
    issue(number:$number){
      timelineItems(first:50,itemTypes:[CLOSED_EVENT]){
        nodes{
          __typename
          ... on ClosedEvent {
            createdAt
            closer {
              __typename
              ... on PullRequest {
                number
                title
                state
                mergedAt
                url
                baseRefName
                headRefName
              }
            }
          }
        }
      }
    }
  }
}
'@
    $result = & $Gh api graphql -f owner=$owner -f name=$name -F number=$Issue -f query=$query | ConvertFrom-Json
    $nodes = $result.data.repository.issue.timelineItems.nodes
    $nodes |
        ForEach-Object { $_.closer } |
        Where-Object { $_ -and $_.__typename -eq "PullRequest" -and $_.state -eq "MERGED" -and $_.baseRefName -eq "main" -and $_.mergedAt } |
        Select-Object -First 1
}

while ($true) {
    $now = Get-Date -Format "yyyy-MM-dd HH:mm:ssK"
    $issueState = Get-IssueState
    $closingPr = Get-ClosingMergedPr

    if ($issueState.state -eq "CLOSED" -and $closingPr) {
        Write-Host "[$now] Issue #$Issue is closed by merged PR #$($closingPr.number) at $($closingPr.mergedAt)."

        git -C $RootCheckout fetch origin main
        git -C $RootCheckout pull --ff-only origin main
        $localMain = git -C $RootCheckout rev-parse --short=12 main
        $originMain = git -C $RootCheckout rev-parse --short=12 origin/main

        if ($localMain -ne $originMain) {
            throw "Root checkout main ($localMain) does not match origin/main ($originMain)."
        }

        git -C $Worktree fetch origin main
        git -C $Worktree switch $Branch
        git -C $Worktree rebase origin/main

        Write-Host "[$now] Ready for issue #68. Root main and origin/main match at $localMain; $Branch rebased onto origin/main."
        exit 0
    }

    Write-Host "[$now] Issue #$Issue is not ready. State=$($issueState.state); closing merged PR found=$([bool]$closingPr). Checking again in $IntervalSeconds seconds."
    Start-Sleep -Seconds $IntervalSeconds
}
