%tmux -CC
** tmux mode started **

Command Menu
----------------------------
esc    Detach cleanly.
  X    Force-quit tmux mode.
  L    Toggle logging.
  C    Run tmux command.
Detached

`tmux new -s session_name`
creates a new tmux session named session_name
`tmux attach -t session_name`
attaches to an existing tmux session named session_name
`tmux switch -t session_name`
switches to an existing session named session_name
`tmux list-sessions`
lists existing tmux sessions
`tmux detach (prefix + d)`
detach the currently attached session

tmux a (to attach to session)
ctrl+a d (to deattach from a session)
