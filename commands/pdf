cpdf infos on http://www.coherentpdf.com/usage-examples.html

pdftk infos on: http://www.pdflabs.com/docs/pdftk-cli-examples/
           and: http://linuxcommando.blogspot.de/2014/01/how-to-split-up-pdf-files-part-2.html

How to merge Pdfs:
pdftk 10011948_57.pdf 10011948_*.pdf cat output landold.pdf

How to cut margins from pdf
pdfcrop --margins '-170 0 0 -135' 20150925130717.pdf out2.pdf   # pixels for left, top, right and bottom 


how to split pdfs:
pdftk infile.pdf cat 51-100 output outfile.pdf
pdftk 1998_PR_Kraftmakher_Equilibrium_Vacancies_Thermal_prop_of_matter_Kraftmakher.pdf cat 5-5 output Kraftmakher_page83.pdf
or:
cpdf in.pdf 1-12 -o out.pdf
cpdf in.pdf 13-end -o out.pdf

how to split each page in the input file into a separate output file: (use the burst command)
pdftk myoldfile.pdf burst 
split in chunks of x pages
cpdf cpdf in.pdf -split -chunk 12 -o out%%%.pdf


pdftk kan split pdfs 
pdftk can rotate all or certain pages
pdftk can rotate only by n*90 degrees where n is an integer
pdftk cann fill out pdfs
Rotate by 19 degrees:
cpdf -rotate-contents 19 in.pdf -o out.pdf

acrobat professional can straighten misrotated pages (which are misrotated by a few degrees) but document looses quality.

remove certain page: (here remove page 10 to page 25)
pdftk myDocument.pdf cat 1-9 26-end output removedPages.pdf

um mit preview pdf's zusammenzufuegen muss man erst alle reinziehen, dann ein pdf auf das andere schieben so dass ein gruenes plus erscheint
