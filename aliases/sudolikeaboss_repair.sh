#!/bin/sh
echo 'from https://support.1password.com/kb/201707/'

echo "close your brwoser (chrome)"
osascript -e 'quit app "Chrome"'

echo 'mkdir -p ~/Library/Application\ Support/Google/Chrome'
mkdir -p ~/Library/Application\ Support/Google/Chrome

echo "quit 1Password and 1Password mini (makt this script to do so by applescript)"
osascript -e 'quit app "1Password"'
osascript -e 'quit app "1Password mini"'

echo "Open 1Password again. (make this script do so)"
open -a 1Password 
echo "Restart your browser."
open /Applications/Google\ Chrome.app/

#cd ~/Library/Application Support/Google/Chrome/NativeMessagingHosts/
#cp 2bua8c4s2c.com.agilebits.1password.json ~/




#  What you should do
#  If you see this message right after updating 1Password, it doesn’t mean your browser lacks support for native messaging.
#  
#  First, try restarting your browser. If you’re still having trouble, follow these steps to allow 1Password to use native messaging in your browser:
#  
#  Open Terminal, which is in the Utilities folder of your Applications folder.
#  Copy and paste the following command and press Return:
#  
#  mkdir -p ~/Library/Application\ Support/Google/Chrome
#  Open 1Password and click the 1Password menu. Hold down the Control and Option keys on your keyboard. While holding the Control and Option keys, select “Quit 1Password and 1Password mini”.
#  
#  Open 1Password again.
#  
#  Restart your browser.
#  
#  Follow the additional steps below if you use Brave, Vivaldi, or Chrome Canary.
#




#
#    If you use Vivaldi or Chrome Canary
#    
#    If you use Chrome Canary or Vivaldi, follow these steps to allow 1Password to use native messaging in your browser:
#    
#    Make sure the 1Password app and your browser are both in the Applications folder.
#    In Finder, choose Go > Go to Folder. Copy and paste the following path and press Return:
#    
#    ~/Library/Application Support/Google/Chrome/NativeMessagingHosts/
#    Copy 2bua8c4s2c.com.agilebits.1password.json to your Desktop.
#    
#    Choose Go > Go to Folder. Copy and paste the correct path for your browser and press Return:
#    
#    Chrome Canary
#    
#    ~/Library/Application Support/Google/Chrome Canary/
#    Vivaldi
#    
#    ~/Library/Application Support/Vivaldi/
#    Open the folder NativeMessagingHosts or create it if it doesn’t exist.
#    
#    Drag 2bua8c4s2c.com.agilebits.1password.json from your Desktop into the NativeMessagingHosts folder.
#    
#    Restart your browser.
