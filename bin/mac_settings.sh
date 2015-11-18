#!/bin/sh

# Set finder highlight color to green
defaults write NSGlobalDomain AppleHighlightColor -string "0.764700 0.976500 0.568600"

# Disable the “Are you sure you want to open this application?” dialog
defaults write com.apple.LaunchServices LSQuarantine -bool false

# Disable Notification Center and remove the menu bar icon
launchctl unload -w /System/Library/LaunchAgents/com.apple.notificationcenterui.plist 2> /dev/null

# Enable full keyboard access for all controls
# (e.g. enable Tab in modal dialogs)
defaults write NSGlobalDomain AppleKeyboardUIMode -int 3

#defaults write NSGlobalDomain KeyRepeat -int 0    # very fast         then: logout and login again
#defaults write NSGlobalDomain KeyRepeat -int 0.5  # current           then: logout and login again
#defaults write NSGlobalDomain KeyRepeat -int 1    # fast              then: logout and login again
#defaults write NSGlobalDomain KeyRepeat -int 2    # normal            then: logout and login again
echo "set default screenshots folder to ~/Dropbox/Screenshots"
defaults write com.apple.screencapture location /Users/glensk/Dropbox/Screenshots/;killall SystemUIServer

echo "Disable press-and-hold for keys in favor of key repeat (if not good, disable) I think this also means no press-and-hold in xmgrace! better not!" 
#defaults write NSGlobalDomain ApplePressAndHoldEnabled -bool false 

echo "Set a blazingly fast keyboard repeat rate" 
#defaults write NSGlobalDomain KeyRepeat -int 0.1

echo "Set a shorter Delay until key repeat" 
#The amount of time it takes for key repetition to begin can also be set by adjusting the "Delay Until Repeat" option under System Preferences => Keyboard (Keyboard tab). Again, if this is still too slow for you (as it was for me), you can set an even faster speed by opening Terminal and typing:
#Where 4, again, can be adjusted (smaller is faster). I'd highly recommend you do not set this option under 4, though, because that would just be impossibly fast (touching a key for a mere split second would type about 10 repeating characters). I ended up setting mine to 7, which might still be too fast for me.
#defaults write NSGlobalDomain InitialKeyRepeat -int 12
#defaults write NSGlobalDomain InitialKeyRepeat -int 11 
defaults write NSGlobalDomain InitialKeyRepeat -int 14

echo "Force the Dock to only show running applications System"
defaults write com.apple.dock static-only -bool TRUE

echo "Disable the switching Spaces animation in the terminal by typing:"
defaults write com.apple.dock workspaces-swoosh-animation-off -bool YES && killall Dock


echo "Disable the sound effects on boot"
#sudo nvram SystemAudioVolume=" "  # this line did not work
sudo nvram SystemAudioVolume=%80


echo "###############################################################################"
echo "# Finder                                                                      #"
echo "###############################################################################"
echo ""
echo "Finder: allow quitting via ⌘ + Q; doing so will also hide desktop icons"
defaults write com.apple.finder QuitMenuItem -bool true

echo "to display full path in finder; Display full POSIX path as Finder window title"
defaults write com.apple.finder _FXShowPosixPathInTitle -bool true

# Enable spring loading for directories
defaults write NSGlobalDomain com.apple.springing.enabled -bool true

# Remove the spring loading delay for directories
defaults write NSGlobalDomain com.apple.springing.delay -float 0

# Avoid creating .DS_Store files on network volumes
defaults write com.apple.desktopservices DSDontWriteNetworkStores -bool true

# Finder: show status bar
defaults write com.apple.finder ShowStatusBar -bool true

# Finder: show path bar
defaults write com.apple.finder ShowPathbar -bool true

# Finder: allow text selection in Quick Look
defaults write com.apple.finder QLEnableTextSelection -bool true

# Automatically open a new Finder window when a volume is mounted
defaults write com.apple.frameworks.diskimages auto-open-ro-root -bool true
defaults write com.apple.frameworks.diskimages auto-open-rw-root -bool true
defaults write com.apple.finder OpenWindowForNewRemovableDisk -bool true

# Disable the warning before emptying the Trash
defaults write com.apple.finder WarnOnEmptyTrash -bool false

# Empty Trash securely by default
defaults write com.apple.finder EmptyTrashSecurely -bool true

# Enable AirDrop over Ethernet and on unsupported Macs running Lion
defaults write com.apple.NetworkBrowser BrowseAllInterfaces -bool true


# Show the ~/Library folder
chflags nohidden ~/Library

echo "####################################################################################"
echo " NOW YOU HAVE TO LOGOUT AND LOGIN AGAIN !!!!!!!!!!" 
echo "####################################################################################"

# Disable Dashboard
defaults write com.apple.dashboard mcx-disabled -bool true

# Don’t show Dashboard as a Space
defaults write com.apple.dock dashboard-in-overlay -bool true

exit
##########################################################################################
# Stuff which has to be run only once and is working now
##########################################################################################


