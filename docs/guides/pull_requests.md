# Reviewing and merging a pull request

A **pull request** (PR) is a proposal to add the changes on one branch to `main`. It
gives each change a description, a list of the files it touches, a test run and a place
for comments, before anything reaches `main`. In this project every backlog task arrives
as its own pull request, and only Chris merges.

*GitHub's page layout changes from time to time; the buttons may look slightly different
from this description.*

## 1. Open it

Go to the repository on GitHub, click **Pull requests**, then the one you want.

## 2. Check the tests

Near the bottom of the **Conversation** tab is a box listing the checks (lint, tests on
Python 3.10 to 3.13, docs). All should show a green tick. A red cross means something
failed: do not merge; ask Claude to look at it.

## 3. Review the changes

- The **Conversation** tab has the description: what changed, how it was checked, what
  to watch for.
- The **Files changed** tab shows every changed line: green added, red removed. You do
  not need to read every line of a large port. Read the description, and look at the
  parts you care about most (for example the physics modules).
- To comment on a line, hover over it and click the **+**. To ask for changes before
  merging, use **Review changes** (top right of the Files changed tab).

## 4. Merge

1. Back on the **Conversation** tab, find the green merge button at the bottom.
2. Click the small arrow beside it and choose **Squash and merge**. This turns all the
   commits on the branch into one commit on `main`, which keeps the history readable:
   one commit per task.
3. Click **Squash and merge**, leave the suggested message as it is (it keeps the
   co-author lines), and **Confirm**.
4. Click **Delete branch** when offered. The branch is finished; its history is kept in
   the merged commit and the pull request.

## 5. After merging

- The tests run again on `main` (the **Actions** tab). Check they pass.
- If you have a copy of the repository on another machine, update it with `git pull`.
- If you installed the package from GitHub, reinstall it (see
  [installing and reprocessing](install_and_reprocess.md)).

## Undoing a merge

A merged pull request can be undone with the **Revert** button on the merged pull
request page. That creates a new pull request that removes the changes, so nothing is
lost.

## Recommended one-off setting (Chris to decide)

GitHub can enforce the "work on a branch, merge by pull request, only when tests pass"
rule. In **Settings → Branches** (or **Settings → Rules → Rulesets**), add a rule for
`main` that requires a pull request before merging and requires the status checks to
pass. This stops anyone, including Claude sessions, from pushing straight to `main`.
