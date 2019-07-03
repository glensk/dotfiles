first: sudo vi /etc/shells and add /usr/local/Cellar/bash/4.3.18/bin/bash

On OSX, use dscl - the Directory Services CLI - to change the current user's shell.

To examine the current user's shell, use:

dscl . -read /Users/$USER UserShell


To change the current user's shell to, e.g. /usr/local/Cellar/bash/4.3.18/bin/bash 

sudo dscl . -change /Users/$USER UserShell /bin/bash /usr/local/Cellar/bash/4.3.18/bin/bash
