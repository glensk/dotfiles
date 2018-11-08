#!/usr/bin/python
import argparse
import numpy as np
import os


def ReadIndex(inputfile):
    """This function reads the input file to extract the index of the first line of each frame to use later"""
    frame_line = []
    print('inputfile',inputfile)
    with open(inputfile) as datafile:
        if inputfile.endswith('.data') or inputfile.endswith('.data.all'):
            for i, line in enumerate(datafile):
                if line.split()[0] == 'begin': # If the datafile is RuNNer-like it will look for the keyword "begin"
                    frame_line.append(i)
        elif inputfile.endswith('.xyz'):
            for i, line in enumerate(datafile):
                if line.split()[0].isdigit(): # If the datafile is xyz, it will look for a line starting with a number in it
                    frame_line.append(i)
        else: raise NameError('Unknown data type. Currently only .data (RuNNer) and .xyz are supported')
    frame_line.append(i+1)  # Appends the last line number + 1. Necessary for the way copy works in the Print_frames function
    return frame_line


def Precomputed(precomp_file, frame_line):
    """The file fed by the user is read and the selected frames are stored to be printed in a new file """
    sel_frames = []
    with open(precomp_file) as precomp: # Reads the provided file and appends the indices of the frames. It could probably be done with np.loadtxt
        for val in precomp.read().split():
            sel_frames.append(int(val))
    if len(sel_frames) > len(frame_line):
        raise ValueError('The number of chosen frames is higher than the total number of frames in the input file')
    sel_frames.sort() # Sorts the selected indices to speed up the writing process. This could be disabled at will but then several cycles of writing are needed in Print_frames
    return np.array(sel_frames)


def Random(nframes, frame_line):
    """A random selection among the frames is done"""
    if nframes > len(frame_line):
        raise ValueError('The number of randomly chosen frames is higher than the total number of frames in the input file')
    rand_frames = np.random.choice(range(len(frame_line)-1), nframes, replace = False)
    rand_frames.sort()
    np.savetxt('rndIndices.idx',rand_frames, fmt='%d', header='This file contains the indices of the randomly selected structures')
    return rand_frames


def Stride(nstride, frame_line): # !TODO Should add option to choose starting and ending frame
    """A selection of every n-th frame"""
    stride_frames = np.arange(0,frame_line[-1],nstride)
    return stride_frames


def Print_frames(inputfile, frames, frame_line, split, prefix):
    """Prints in a new file (or two new files, with the splitting flag active) the chosen structures. If no prefix is given, the standard output is DATAFILE_selected.EXTENSION"""
    if inputfile.endswith('.data') or inputfile.endswith('.data.all'):
        suffix = '.data'
        if prefix == '':
            prefix = inputfile[:-5]
    elif inputfile.endswith('.xyz'):
        suffix = '.xyz'
        if prefix == '':
            prefix = inputfile[:-4]
    sel_input = open(prefix+'_selected'+suffix,'w') # Maintains the same extension as the original datafile
    if split:
        rem_input = open(prefix+'_remaining'+suffix,'w')
    copy = 0 # Represents the lines to copy. If negative the lines will not be printed or will be printed in the split file. If positive will copy in the "_selected" file
    counter = 0 # Counter for the frames that are being copied
    frames = np.append(frames,len(frame_line)-1) # Appends the last frame to force the function to read the inputfile until the last line. Less efficient but necessary if one wants to split the dataset
    printed_frame = frames[counter]

    with open(inputfile) as datafile:
        for i, line in enumerate(datafile):
            if i == frame_line[printed_frame]:
                copy = frame_line[printed_frame+1] - frame_line[printed_frame]
                counter += 1
                printed_frame = frames[counter]
            if copy > 0:
                sel_input.write(line)
            if split and copy <= 0:
                rem_input.write(line)
            copy -= 1

    sel_input.close()
    if split:
        rem_input.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)

    parser.add_argument("inputfile", type=str, help="The input file out of which the structures are going to be selected")
    parser.add_argument("--split", action='store_true', help="If active, it splits the initial datafile in two different files")
    parser.add_argument("--prefix", type=str, default='', help="Prefix for output files (defaults to input file name)")
    subparsers = parser.add_subparsers(help='Choose either random or precomp selection of the datapoints', dest='mode')

    parser_random = subparsers.add_parser('random', help='Random selection of frames. Must be followed by the number of frames you want to choose')
    parser_random.add_argument('N_rand', type=int, help='Number of randomly chosen frames')
    parser_precomp = subparsers.add_parser('precomp', help='Precomputed selection of frames. A file containing the indices of the frames to be selected must be provided. It is 0-based')
    parser_precomp.add_argument('precomp_file', type=str, help='The path to the file containing the indices of the frames')
    parser_stride = subparsers.add_parser('stride', help='Select every N frames. N must be specified after calling this option')
    parser_stride.add_argument('N_stride', type=int, help='Every how many frames you want to print the structures')

    args = parser.parse_args()
    frame_line = ReadIndex(args.inputfile)
    if args.mode == 'precomp':
        frames = Precomputed(args.precomp_file, frame_line)
    elif args.mode == 'random':
        frames = Random(args.N_rand, frame_line)
    elif args.mode == 'stride':
        frames = Stride(args.N_stride, frame_line)
    Print_frames(args.inputfile, frames, frame_line, args.split, args.prefix)
