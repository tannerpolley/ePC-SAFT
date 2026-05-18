from pathlib import Path


SOURCE_PATTERNS = (
    "docs/papers/md/Hubach et al. - 2024 - Li+ Extraction*.md",
    "docs/papers/md/Hubach et al. - 2024 - Supporting Information*.md",
)


def main() -> int:
    repo_root = Path(__file__).resolve().parents[5]
    missing = [
        pattern
        for pattern in SOURCE_PATTERNS
        if not list(repo_root.glob(pattern))
    ]
    if missing:
        print("Hubach 2024 lane is scaffolded but blocked.")
        print("Missing source markdown in this worktree:")
        for pattern in missing:
            print(f"  - {pattern}")
        print("Normalize source tables into data/input before enabling this lane.")
        return 1

    print("Hubach 2024 source assets are present, but the executable model route is not implemented in this lane yet.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

