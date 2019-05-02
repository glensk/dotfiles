#!/bin/sh

cd /Users/glensk/.ssh

#If you get ""-bash: cd: .ssh: No such file or directory"
# You need to create it
#mkdir /Users/glensk/.ssh

# Set security on the folder so only you can read, write and execute
chmod 700 /Users/glensk/.ssh

# Change directory
#cd .ssh

# Create the config file
#touch config

# Edit the config file
#vi config

# When in vi hit i to enter insert mode and add the following lines to the config file
#ServerAliveCountMax 3
#ServerAliveInterval 10  (had it before at 120)

# Save and quit by hitting :wq

# Set correct permissions on the config file
chmod 644 /Users/glensk/.ssh/config

