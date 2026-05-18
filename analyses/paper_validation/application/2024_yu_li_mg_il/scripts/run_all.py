from pathlib import Path


SOURCE_PATTERNS = (
    "docs/papers/md/Yu - 2024 - Highly efficient lithium extraction*.md",
    "docs/papers/md/Yu - 2024 - Supporting*.md",
)


def main() -> int:
    repo_root = Path(__file__).resolve().parents[5]
    missing = [
        pattern
        for pattern in SOURCE_PATTERNS
        if not list(repo_root.glob(pattern))
    ]
    if missing:
        print("Yu 2024 lane is scaffolded but blocked.")
        print("Missing source markdown or supplementary tables in this worktree:")
        for pattern in missing:
            print(f"  - {pattern}")
        print("Normalize source tables into data/input before enabling this lane.")
        return 1

    print("Yu 2024 source assets are present, but the executable model route is not implemented in this lane yet.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

