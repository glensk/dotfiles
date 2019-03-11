########################
Set Up Git
########################
yum install git                                                             #get git
git config --global user.name "glensk"                                      #set username
git config --global user.email "albert.glensk@gmail.com"                    #set email
git config --global credential.helper cache                                 #set Password caching you don't want to type your username and password every time
git config --global credential.helper 'cache --timeout=3600'                #By default git will cache your password for 15 minutes. You can change this if you like


git config --global user.name "gpleyson"
git config --global user.email "dagatkarimlan@gmail.com"

# set default so that all changes are always pushed to the repository
git config --global push.default "matching"

# set default so that you avoid unnecessary commits
git config --global branch.autosetuprebase always 

git config --global color.status auto
git config --global color.branch auto 

git config --global core.editor vim 

To query your Git settings, execute the following command:
git config --list 

# Configure Git to use this file
# as global .gitignor
git config --global core.excludesfile ~/.gitignore 

########################
########################
git commit -a -m "Initial commit"   # with -a add newstuff should not be necessary
git remote update
# to push use the command:
# git push [target]
# default for [target] is origin
# git push [remote-name] [branch-name=master]
git push origin 
git push origin master

# show the details of the remote repo called origin
git remote show origin 

# show the existing defined remote repositories
git remote

# pull in the latest changes of your remote repository
git pull

# To see remote repos(itories)
git remote -v

# schow branches
git branch -a

############### DONE #######################
git branch --set-upstream master origin/master
git branch --set-upstream master remotes/origin/master

#################################
# committing
#################################
git clone https://glensk@bitbucket.org/glensk/scripts.git
git clone https://glensk@bitbucket.org/glensk/scripts.git/Thermodynamics
git add -u  # so that it will only stage the modified files.
git commit -a   # commit only the modifications.

###################
git remote update
git pull to get latest version
#############
git add -A
git commit -m "txt"
git push
#############

der .git folder des repositories kann einfach in ein bestehendes svn repository kopiert werden.... wenn alles glatt laeuft kann weiterhin git und svn genutzt werden
mac_tools


##############
git branch
git reset --hard origin/master  ## when you see *(no branch) and mater underneeth
git checkout master             ## to get back to master branch when on *(no branch)

git push --force origin
git pull
git rebase --abort


### git status -u gives: 
# On branch master
# Your branch and 'origin/master' have diverged,
# and have 3 and 2 different commits each, respectively.
#
nothing to commit (working directory clean)

git submodule add https://github.com/jcf/vim-latex.git dotfiles/vim/bundle/vim-latex



###################################################################################
git co https://glensk@bitbucket.org/glensk/scriptsnew.git
###################################################################################
#
# git get remote url: 
git config --get remote.origin.url

###################################################################################
aiida-alloy
###################################################################################
https://gitlab.com/daniel.marchand/aiida-alloy
pw: mit zwei doppel s nix weiter
