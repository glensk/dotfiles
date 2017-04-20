find . -type f \( -name "lifetimesgood_*" -o -name "freqsgood_*" \)     # sucht nach 2 files
find . -name "file" -mtime +2                                           # sucht nach fiel weleche aelter als 2 tage ist
find . -name "RESTAR*"          					                    # sucht nach dateien mit dem Muster (-name) RESTAR  
find . -iname -name "RESTAR*" 						                    # [-iname waere ohne gros/kleinschreibung]
find . -name "*.jar" 							                        # finde datei im jetzigem ordner die .jar im titel hat
find -L -name "OUTCAR*" 						                        # durchsucht auch symbolische links


find /pfad/ -name "RESTAR*" -exec rm {} \;        			            # loescht alle gefundenen dateien in /pfad/...
find . -name "ANMERKUNG" -exec mv {} {}.txt \;                          # benennt alle ANMERKUNG files um zu ANMERKUNG.txt
find */KPOINTS -exec cp ~/pfad/KPOINTS {} \;				            # ersetzt alle gefundenen files mit der aus pfad

find -L . -name "*_*_*-shift" -type d		                            #findet auch alle ordner, * muss in "" stehen
find . -type d -perm -o=w 			                                    #finds only directories with permittion others: writable

find . -name DISP.phon | sed 's|\(.*\)DISP.phon|mv & \1DISP|' | sh


find . -maxdepth 1 -mindepth 1 -type d -name "ergebnis*" -print0 | sed 's|./ergebnis| |g'
		print0 prints erverything in one line
find -L 4.04Ang_250K 4.04Ang_500K 4.04Ang_700K 4.04Ang_850K 4.04Ang_934K -mindepth 2 -maxdepth 2 -name POSCAR -print -quit
ls -1d *Ang | xargs -n 1 cp INCAR

find . -type f \( -name "lifetimesgood_*" -o -name "freqsgood_*" \) -mtime +2-mtime +2

