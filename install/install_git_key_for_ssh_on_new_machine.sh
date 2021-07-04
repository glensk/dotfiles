#!/bin/sh
echo
echo in the next field use a passphrase
echo

ssh-keygen -t ed25519 -C "albert.glensk@gmail.com"
echo 
echo from https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#adding-your-ssh-key-to-the-ssh-agent
echo
echo "Note: The -K option is Apple's standard version of ssh-add, which stores the passphrase in your keychain for you when you add an ssh key to the ssh-agent. If you chose not to add a passphrase to your key, run the command without the -K option. If you don't have Apple's standard version installed, you may receive an error. For more information on resolving this error, see \"Error: ssh-add: illegal option -- K.\""
ssh-add -K ~/.ssh/id_ed25519

echo ""
echo "add the key to github (follow following link)"
echo ""
echo "https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account"
