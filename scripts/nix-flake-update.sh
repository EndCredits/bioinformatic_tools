#!/usr/bin/env bash
# Recursively find all flake.nix and update their flake.lock to use the
# same nixpkgs revision.
#
# Usage:
#   ./scripts/nix-flake-update.sh           # update all
#   ./scripts/nix-flake-update.sh --check   # check if all locks are in sync

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

CHECK_ONLY=false
if [[ "${1:-}" == "--check" ]]; then
  CHECK_ONLY=true
fi

# Find all directories containing a flake.nix
FLAKE_DIRS=()
while IFS= read -r flake; do
  FLAKE_DIRS+=("$(dirname "$flake")")
done < <(find . -name flake.nix -not -path './.git/*' | sort)

if [[ ${#FLAKE_DIRS[@]} -eq 0 ]]; then
  echo "No flake.nix files found."
  exit 0
fi

echo "Found ${#FLAKE_DIRS[@]} flake(s):"
for dir in "${FLAKE_DIRS[@]}"; do
  echo "  $dir"
done
echo ""

# Use the root flake.lock as the source of truth for nixpkgs revision
ROOT_LOCK="./flake.lock"
if [[ ! -f "$ROOT_LOCK" ]]; then
  echo "No root flake.lock found. Generating from root flake.nix..."
  nix flake update
fi

# Extract nixpkgs rev from root lock
ROOT_REV=$(jq -r '.nodes.nixpkgs.locked.rev' "$ROOT_LOCK")
ROOT_HASH=$(jq -r '.nodes.nixpkgs.locked.narHash' "$ROOT_LOCK")
echo "Target nixpkgs: ${ROOT_REV:0:12}..."

OUT_OF_SYNC=0
for dir in "${FLAKE_DIRS[@]}"; do
  lock_file="$dir/flake.lock"
  if [[ ! -f "$lock_file" ]]; then
    echo "  [missing] $dir"
    OUT_OF_SYNC=1
    continue
  fi

  child_rev=$(jq -r '.nodes.nixpkgs.locked.rev // "none"' "$lock_file")
  if [[ "$child_rev" == "$ROOT_REV" ]]; then
    echo "  [  ok  ] $dir"
  else
    echo "  [drift ] $dir (rev: ${child_rev:0:12})"
    OUT_OF_SYNC=1
  fi
done

if $CHECK_ONLY; then
  exit $OUT_OF_SYNC
fi

# Update all flakes to pin to the same nixpkgs as root
if [[ $OUT_OF_SYNC -eq 0 ]]; then
  echo ""
  echo "All flakes in sync."
  exit 0
fi

echo ""
echo "Updating out-of-sync flakes..."

# Build override args: pin nixpkgs input to the root's revision
OVERRIDE_ARGS=()
if [[ -f "$ROOT_LOCK" ]]; then
  NIXPKGS_REF=$(jq -r '.nodes.nixpkgs.original.ref // empty' "$ROOT_LOCK")
  NIXPKGS_REV=$(jq -r '.nodes.nixpkgs.locked.rev' "$ROOT_LOCK")
  NIXPKGS_OWNER=$(jq -r '.nodes.nixpkgs.locked.owner' "$ROOT_LOCK")
  NIXPKGS_REPO=$(jq -r '.nodes.nixpkgs.locked.repo' "$ROOT_LOCK")
  OVERRIDE_ARGS=(--override-input nixpkgs "github:${NIXPKGS_OWNER}/${NIXPKGS_REPO}/${NIXPKGS_REV}")
fi

for dir in "${FLAKE_DIRS[@]}"; do
  lock_file="$dir/flake.lock"
  child_rev=$(jq -r '.nodes.nixpkgs.locked.rev // "none"' "$lock_file" 2>/dev/null)

  if [[ "$child_rev" != "$ROOT_REV" ]]; then
    echo "  Updating $dir ..."
    (cd "$dir" && nix flake update "${OVERRIDE_ARGS[@]}" --commit-lock-file 2>&1 | sed 's/^/    /')
  fi
done

echo ""
echo "Done. All flakes pinned to nixpkgs ${ROOT_REV:0:12}."
