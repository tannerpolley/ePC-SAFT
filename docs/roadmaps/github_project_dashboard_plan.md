# GitHub Project Dashboard Plan

## Decision

Use GitHub Projects as a visibility dashboard for the general reactive/electrolyte roadmap.

Do not use Project status alone as permission to start coding.

The actual start gate is:

```text
issue closed
PR merged
origin/main contains merge commit
branch rebases cleanly onto origin/main
dependency_gate.yaml passes
```

## Project name

```text
ePC-SAFT — General Reactive/Electrolyte Roadmap
```

## Status field values

```text
Draft
Goal Prep Ready
Blocked on Dependency
Ready to Start
In Progress
PR Open
Needs Review
Merged
Stopped / Needs Handoff
Deferred
```

## Recommended fields

```text
Status
Task Key
Phase
Branch
Depends On
PR
Agent
Priority
Risk
Last Checked
Blocker
```

## Setup using GitHub CLI

The setup agent should check what project commands are available:

```powershell
gh project --help
gh project create --help
gh project field-create --help
gh project item-add --help
```

If supported, create the project:

```powershell
gh project create --owner tannerpolley --title "ePC-SAFT — General Reactive/Electrolyte Roadmap"
```

Then create fields if the local GitHub CLI supports those commands.

If GitHub Project commands are unavailable, the setup agent should:

```text
create the issues
create docs/roadmaps files
create issue drafts
record that project creation must be done manually
```

## Manual fallback

If Projects cannot be created from the agent, create the GitHub Project manually and add all created issues A-M to it.

Use status values exactly as listed above.
