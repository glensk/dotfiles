added to /Users/glensk/anaconda/lib/python2.7/site-packages/IPython/core/magics/osm.py


    @skip_doctest
    @line_magic
    def mcd(self, folder):
        os.makedirs(folder)
        os.chdir(folder)


