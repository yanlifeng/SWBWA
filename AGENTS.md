## Long-running task completion notification

For every user turn, as your first shell action, run:

TOKEN="$($HOME/.codex/notify_vscode.sh start)" && echo "NOTIFY_TOKEN=$TOKEN"

Remember the returned NOTIFY_TOKEN for this turn.

After all requested work, testing, verification, and other actions are fully complete,
and immediately before sending the final response, run exactly once:

$HOME/.codex/notify_vscode.sh done "<NOTIFY_TOKEN>"

Replace <NOTIFY_TOKEN> with the token returned at the beginning of this turn.

Do not run the completion command while work is still in progress, while waiting
for a command, while waiting for user input, or while waiting for approval.

Notification failure must not affect the task result and must not be retried.
