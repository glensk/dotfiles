##########################################################
svn co http://cmcc1.mpie.de/data/pyCMD ~/svn
##########################################################

#####################################################################################
svn co http://cmcc1.mpie.de/data/Thermodynamics    #checkout == update/download # at cmpc # first time give username and pw (it is the password for the wiki/not the computer account cluster)
svn co --username glensk http://localhost:48540/data/Thermodynamics  # at cmmc
#####################################################################################

#####################################################################################
svn co --username guest svn://svn.ifm.liu.se/tdep/trunk .
cozo-ceci-do-vow
#####################################################################################

#### how did i add outcar folder
#### while being in ~/scripts/Thermodynamics
cp ~/scripts/outcar ~/scripts/Thermodynamics/outcar
svn add outcar
svn commit -m "message"
svn status
svnversion . ## make in scripts folder to get revisio number

#############################################################
svn co https://subversion.assembla.com/svn/scriptssvn/ scripts
#############################################################

svnignore
svn propset svn:ignore -F .gitignore .
You'll need to re-read the file whenever it's updated or changed (and the value will be overwritten each time).

svn folder kann beliebig umbenannt werden ... alle funktionalitaet bleibt erhalten
<<<<<<< HEAD
svn diff -r 30 dotfiles/subversion/config

if conflict: missng file/folder: add file/folder and make: svn revert file/folder    details at: http://stackoverflow.com/questions/4317973/svn-how-to-resolve-local-edit-incoming-delete-upon-update-message
=======
>>>>>>> d1e809b9529c5ce2a53db7a30dad93b231333b12


!     C              ti_dominique/manual_FORCES_SETS
      >   local delete, incoming delete upon update

--> svn revert ti_dominique/manual_FORCES_SETS


########################################################
https://publication.icams.rub.de/svn/RUW_defects
glensk
Vg1h0loW
###########################################################
svn co https://subversion.assembla.com/svn/scripts_svn/trunk scripts 


################################################################
svn revert --depth infinity *
svn co https://subversion.assembla.com/svn/thermodynamics_/ thermodynamics
in python folder do svn co https://subversion.assembla.com/svn/thermodynamics_/trunk thermodynamics
