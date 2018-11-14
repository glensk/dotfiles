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
@click.command(context_settings=CONTEXT_SETTINGS)


click.option('-nsi'  ,required=True, prompt=True, type=int, help="number of Si atoms")

#filename="input.data.all"
#rm -f $filename
#touch $filename
#files=`find . -maxdepth 4 -name simulation.pos_0.xyz`
#for i in $files;do
#    echo $i
#        xyz2runner.sh $i >> $filename
#        done

def gather_xyz(nsi):
    ''' This scipt looks for all simulation.pos_0.xyz files in the $searchfolder and converts
    those to one $outputfile
    '''
    pass

if __name__ == "__main__":
    gather_xyz()

