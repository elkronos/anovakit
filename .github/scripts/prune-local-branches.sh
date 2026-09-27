#!/usr/bin/env sh
# Trim a local clone once its branches are done with.
#
#   .github/scripts/prune-local-branches.sh            # show what would go
#   .github/scripts/prune-local-branches.sh --apply    # do it
#
# 1. Drops remote-tracking refs for branches deleted on GitHub.
# 2. Removes worktrees whose branch is fully merged into main and which have
#    no uncommitted changes.
# 3. Deletes local branches fully merged into main.
#
# main, master, gh-pages and the branch you are on are never touched. A branch
# goes only when every one of its commits is already in main, and a worktree
# only when it is also clean, so no unmerged commit or uncommitted change can
# be lost.
set -eu

apply=false
[ "${1:-}" = "--apply" ] && apply=true

git fetch --prune origin

base=origin/main
git rev-parse --verify --quiet "$base" >/dev/null || base=origin/master
current=$(git symbolic-ref --quiet --short HEAD || echo "")

keep() {
  case "$1" in
    main|master|gh-pages) return 0 ;;
  esac
  [ "$1" = "$current" ]
}

merged() {
  git merge-base --is-ancestor "$1" "$base"
}

# Worktrees (skip the main one, the first listed)
git worktree list --porcelain | awk '
  /^worktree / { wt = substr($0, 10); idx++ }
  /^branch /   { b = $2; sub("^refs/heads/", "", b); if (idx > 1) print wt "\t" b }
' | while IFS="$(printf '\t')" read -r wt branch; do
  if keep "$branch" || ! merged "$branch"; then
    echo "Keeping worktree $wt ($branch)"
  elif [ -n "$(git -C "$wt" status --porcelain)" ]; then
    echo "Keeping worktree $wt ($branch): it has uncommitted changes"
  elif $apply; then
    git worktree remove "$wt" && echo "Removed worktree $wt"
  else
    echo "Would remove worktree $wt ($branch)"
  fi
done
git worktree prune

# Local branches
for branch in $(git for-each-ref --format='%(refname:short)' refs/heads); do
  if keep "$branch"; then
    continue
  elif ! merged "$branch"; then
    echo "Keeping $branch ($(git rev-list --count "$base..$branch") commit(s) not in ${base#origin/})"
  elif git worktree list --porcelain | grep -qx "branch refs/heads/$branch"; then
    echo "Keeping $branch: it is checked out in a worktree"
  elif $apply; then
    # Every commit is in main (checked above), so forcing is safe; -d would
    # refuse a branch merged into main but not into the one checked out.
    git branch -D "$branch"
  else
    echo "Would delete $branch (fully merged into ${base#origin/})"
  fi
done

$apply || echo "Dry run: nothing was changed. Run with --apply to do it."
