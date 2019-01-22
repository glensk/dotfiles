#!/bin/sh
commands=`ls -1d $dotfiles/commands/* | sed 's|.*commands/||g'`
rm -f $dotfiles/commands/my_commands.txt
touch $dotfiles/commands/my_commands.txt

echo "$commands" > $dotfiles/commands/my_commands.txt

