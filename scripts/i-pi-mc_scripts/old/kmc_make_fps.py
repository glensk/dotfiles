#!/usr/bin/env python
import click

# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init

# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.option('-i','--infile',required=False,default="input.data.all",help="Path to inputfile in RuNNer format.")

def make_fps(infile):
    '''
    Given an infile ---which is supposed to have a runner input format---, e.g.  a fps is made.
    '''
    if not os.path.isfile(infile):
        sys.exit('missing file '+infile)

