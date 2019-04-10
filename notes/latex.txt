i> ./BetaT2.tex:271:  ==> Fatal error occurred, no output PDF file produced!
> Transcript written on BetaT2.log.
> execution error: Can’t make file ":BetaT2.pdf" into type alias. (-1700)
> ~
> Richard-Seguins-MacBookPro:~ richardseguin$ 


Howdy,

Did you create an updmap.cfg file for the MinionPro.map (or whatever it's called) and run

sudo -H mktexlsr
sudo -H updmap-sys

(as discussed in TeXLive2012Changes available at <https://dl.dropbox.com/u/10932738/index.html>) or run

sudo -H updmap-sys --enable Map MinionPro.map

to enable the map file?

Good Luck,

Herb Schulz


###########################################################
# description of packages
###########################################################
-->> \usepackage{ulem}

     folgende Befehle zu Verfügung:
        
        \uline{important}  % unterstreichen
        \uuline{urgent}    % doppelt unterstreichen
        \uwave{boat}       % unterschlängeln
        \sout{wrong}       % durchstreichen
        \xout{removed}     % ausstreichen mit //////.

    
-->> \usepackage{color}

    Die möglichen neuen Befehle rem und add um die Textabschnitte, die entfernt bzw. hinzugefügt werden müssen, deutlicher zu markieren, werden exemplarisch durch
        
        \newcommand{\rem}[1]{\textcolor{red}{\sout{#1}}}
        \newcommand{\add}[1]{\textcolor{blue}{\uline{#1}}}

        eingeführt. Diese können ab sofort in einer tex-Datei eingesetzt werden:
            
            This is a text line. \rem{This text was removed.} \add{This text added.} This is another text line


-->> \usepackage{multirow}
    To enable multirows in Tables


###########################################################
# Sonderzeichen
###########################################################
\AA             Angstrom

###########################################################
# citations 
###########################################################
to get correct citations instead of ? : run in command line: bibtex Paper    # where you have the Paper.txt file


###########################################################
# make citations clickable include
###########################################################
\usepackage{hyperref}
\hypersetup{
        linktoc=all,
        colorlinks=true,
        linkcolor=black,
        citecolor=black,
        filecolor=black,
        pagecolor=black,
        urlcolor=blue,
        bookmarks=true,
        pdfborder={0 0 0 0},
        pdftitle={\dissTitle},
        pdfauthor={Albert Glensk},
        pdfkeywords={Finite temperature DFT, Ab initio, Thermodynamics, Anharmonicity, Local Anharmonic Approximation, Point defects},
        pdfsubject={Finite temperature DFT, Ab initio, Thermodynamics, Anharmonicity, Local Anharmonic Approximation, Point defects},
        pdfdisplaydoctitle=true,
        pdftoolbar=true,                % Anzeigen der Acrobat toolbar oder nicht
        pdfmenubar=true,                % Anzeigen des Acrobat menu oder nicht
        bookmarksopen=true,
        bookmarksnumbered=true,
        % pdfstartview=X Y Z ! DONT USE DESTROYS TO OPEN WITH ACROBAT
        }
\definecolor{maroon}{cmyk}{0,0.87,0.68,0.32}
