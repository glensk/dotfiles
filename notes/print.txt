lpr -Pcmcopy1166 dat.jpg  ## to print jpg file
lpq -Pcmcopy1166          ## to show queing list for this printer
lpq                       ## que for default printer
lprm -Pcmcopy1166 dat.jpg ## remove job from spooler

lpoptions -d cmcopy1166   ## set default printer 
setenv PRINTER cmcopy1166 ## in tcsh;      (export PRINTER=cmcopy1166 in bash)


display to show plot

172.16.6.246
cmcopy1166

--------------------------------------------------------------
CUPS
-------------------------------------------------------------
lpstat -p -d    #show printer
lpr -P cmcopy1166 ~/print.pdf   ## does not work
lp -d cmcopy1166 ~/print.pdf    ## 

vi ~/.cups/lpoptions ---> reinschreiben
Default cmcopy1166


--------------------------------------------------------------
mac
-------------------------------------------------------------
under system preferences got to printer
than add ip (since in defaults the cmcopy1166 did not appear)
ip reiter
add ip: 192.12.81.35
everything else is fine =-> printer works now
    Protocol: Was set : Line Printer Daemon -LPD
    Queue: EMPTY! (Leave blank for default queue)
    Use: Generic PostScript Printer

    cmfp34    Basement CM-Building:Buildingcmfp34.mpie.de     192.12.81.34

    cmfp35    2. Floor CM-Building:Buildingcmfp35.mpie.de     192.12.81.35

