
    ## unison.sh -ui text mpie_strato  !!!!!!! mann muss aber in einem der folder (articles/talks/scripts) sein!!!!!!!!!!!!!!!!!!!!!!!
    ## start from talks: works for: scripts


    # The algorithm that Unison uses to convert a (possibly) relative path
    # into an absolute path ("canonization") works by saving the current
    # working directory, changing the working directory to this path and
    # finally restoring the working directory. Apparently, this last step
    # fails in your case.
    
    ### --> ABHILFE: einfach in irgendeinen folder gehen der nicht /home/glenks/ist!
    
    # You should try to change the working directory to some safe place
    # before invoking unison:
    # 
    # #!/bin/bash
    # cd /cygdrive/c
    # /usr/bin/unison environment -debug all
    # --> THIS WORKS!!!
    
    # Roots of the synchronization
    root = /home/glensk
    root = /data/glensk/strato

   rshargs = -C
    mountpoint = talks

    # Paths to synchronize 
    path = articles
    path = talks
    path = scripts
    path = .unison/mpie_strato.prf
    path = .unison/unison.sh
    path = .tcshrc
    path = .vimrc
    path = .xmgracerc

    # Einstellungen
    auto = true                
    #batch = true            
    #silent = true           
    
    # auto = automatically accept default (nonconflicting) actions
    # batch = ask no questions at all
    # silent = print nothing except error messages
    
    # Versuche
    #fastcheck = true






#########################################
## set by some user
#########################################
#  ## Data not to be synchronized
#  ignore = Path home/sa/mm/audio/music
#  ignore = Path home/sa/.adobe/Acrobat/8.0/Synchronizer
#
#  ## Miscellaneous settings
#  rshargs = -C
#  auto =true
#  confirmbigdeletes = true
#  perms = -1
#  owner = true
#  group = true
#  times = true
#  #force = newer
#  sortbysize = true
#  sortnewfirst = true
#  maxthreads = 50
#  log = true
#  logfile = /home/sa/.unison/unison.log

#########################################
## set by some other user
#########################################
#mkdir ~/.unison
#touch ~/.unison/default.prf
#    INHALT:
#    cat <<END_ENTRY >> ~/.unison/default.prf
#    merge = diff3 -m CURRENT1 OLD CURRENT2 > NEW
#    backup = Name *
#    maxbackups = 10
#    log = true
#    logfile = /home/user01/.unison/unison.log
#    rshargs = -C
#    END_ENTRY
#

#unison -ui text /home/glensk/vortraege /data/glensk/strato/vortraege   ## -ui text = commandline   -ui graphic = GUI
#unison -ui text /home/glensk/talks/ /data/glensk/strato/talks/
#unison -ui text /home/glensk /data/glensk/strato -path talks -path article
#
#
#unison -ui text /home/glensk/talks/ /data/glensk/strato/talks/ -path blazej_abstracts/ -copymax 1
#

#rshargs = -C
#mountpoint = cv
