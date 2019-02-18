##################################################
## snippets to insert  (copy not working so far)
##################################################
\f{Times-Italic}T\f{Times}\v{.4}\z{.52}\h{.15} melt\v{-0.94}\h{-2.2}exp
\f{Times-Italic}T\f{Times}\v{.4}\z{.52}\h{.15} melt\v{-0.94}\h{-2.2}exp\N / \f{Times-Italic}T
T\Smelt\N/T (1/K)

bei dem mg2si paper: erstmal pixel eingestellt -> dicke der striche ; die dpi machen dasss bild nur groesser
        in xmgrace aber haben ansonsten ueberhaupt keine auswirkung


grep "energy  w" OUTCAR | xmgrace -block - -bxy 0:4
grep POSITION OUTCAR -A38 | grep -v -e POS -e -- | xmgrace -block - -bxy 1:2


y=(s2.y-(2/3)*s0.y-(1/3)*s1.y)*0.096485341
Formation Enthalpy Mg\s2\NSi  [F(T)\sMg2Si\N - (2/3)F(T)\sMg\N - (1/3)F(T)\sSi\N]



The new way to add symbols in xmgrace is:
press Ctrl-E while you are in the text field: this will bring up the font dialog box -> select the Symbol font -> click the square root.

x -hardcopy dG_Mg2Si.xmg

\+ to increase Fontsize (textsize)
\- to decrease Fontsize (textsize)
\u begin unterlining, to stop underlining use \U
\f{Symbol}D\f{}G / kJ*mol\S-1
\f{Symbol}D\f{}           	# capital DELTA
\f{Symbol}G\f{}           	# capital Gamma  --> phonon dispersion
\f{Symbol}b			        # beta klein
\f{Symbol}l			        # lambda klein
\f{Symbol}q\f{}		        # Greek letters, example: theta
\cE\C				# Special symbols, example: Angstrom symbol
\S          writes Superscript
\N          writes Normal
\s          writes subscript
\u          writes underlined

#############################################
## einruecken Text
#############################################
\s\h{-1.55}   writes subscript further left
\s\v{-1.55}   writes subscript further down
\s\z{-1.55}   writes subscript fontsize

x -block thermo -bxy 1:3			print spalte1:spalte3


xmgrace -nxy
xmgrace ~/xmgraceTemplates/PhononDisp_hcp.agr -nxy dispersion.dat

- symbol char: 45


print directly to printer: device: postscript name:kprinter --> dann auf print



[928] 21:19 glensk@cmpc05 [~/v/Mg2Si--/flourite__]    
x 1 style2plot addthisadd 

x -settype xydy dUdL

xmgrace -settype xydy dUdL avg_dUdL -settype xy fit_*
#############################################
## axis rescale 
#############################################
unter axes: Reiter Tick labels: Extra: Tick labels: Axis transform kann man 1.1*$t eingeben, das skaliert nur die achse, nicht den graphen!

schnell kann man bilder zusammensetzen wenn man erst alles reinlaedt, dann die functionen die zusammen in ein graphen kommen exportiter und dann im andrem bild importiert

gibbs energy of formation vac:

y=exp(-s0.y/(0.08617332478*x))

#####################################################################
print: chooose PostScript
type: lpr -cmcopy1166

###################################################################
HARDCOPY DEVICE "PostScript"
DEVICE "JPEG" DPI 144
xmgrace -settype xydy file


# to make tikmarks at left (or right) side longer/make them show into Graph or line
# 	-> Viev, Axes properties
# 	-> Specia tab
# 	-> write the marks you want to have (Special ticks: Tick marks and labels)
# 	-> make an empty stirng in case of no text
#

# to send pic to printer:
select HARDCOPY DEVICE "PostScript" in ~/.xmgracerc
then: xmgrace -hardcopy file.agr

