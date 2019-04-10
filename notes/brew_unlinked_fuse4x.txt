
==> Summary
🍺  /usr/local/Cellar/sqlite/3.8.2: 9 files, 2.0M
==> Upgrading sshfs
==> Installing sshfs dependency: osxfuse
Error: Cannot install osxfuse because conflicting formulae are installed.

  fuse4x: because both install `fuse.pc`

Please `brew unlink fuse4x` before continuing.

Unlinking removes a formula's symlinks from /usr/local. You can
link the formula again after the install finishes. You can --force this
install, but the build may fail or cause obscure side-effects in the
resulting software.
                                                                                                                                                                                                       
20:34 glensk@mac /Users/glensk%brew unlink fuse4x
Unlinking /usr/local/Cellar/fuse4x/0.9.2... 6 links removed
