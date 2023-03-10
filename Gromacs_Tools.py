#!/usr/bin/env python
# coding: utf-8
__author__="Marcello Miceli"
__version__="0"
import sys
import argparse
import warnings
import os
import fileinput
import glob
import shutil
import subprocess
warnings.filterwarnings('ignore')




if shutil.which('gmx') is None:
    print('\n{0}\nThe program gmx is not present in the current environment\n{0}\n'.format('#'*80))
    sys.exit(10)


def parse_options():
    parser = argparse.ArgumentParser(description='''Gromacs mdrun and grompp wrapper''')
    parser.add_argument('--pdb',
                        metavar='complex.pdb',
                        dest='str',
                        type=str,
                        help='Input protein-ligand pdb file.',
                        required=True)
    parser.add_argument('--top',
                        metavar='topol.top',
                        dest='top',
                        type=str,
                        help='Input protein-ligand topology file.',
                        required=True)
    parser.add_argument('--trj',
                        metavar='trajectory.xtc',
                        dest='xtc',
                        type=str,
                        help='Input protein-ligand trajectory file, without any solvent molecule.',
                        required=True)
    parser.add_argument('--outDir',
                        metavar='outdir',
                        dest='outDir',
                        type=str,
                        help='Path to the output directory.',
                        required=True)
    parser.add_argument('--start',
                        metavar='ss',
                        dest='startFrame',
                        type=int,
                        help='Starting frame of the analysis window.',
                        required=True)
    parser.add_argument('--end',
                        metavar='ee',
                        dest='endFrame',
                        type=int,
                        help='Last frame of the analysis window.',
                        required=True)
    parser.add_argument('--windows',
                        metavar='nw',
                        dest='windows',
                        type=int,
                        help='Number of mm/gbsa calculations to perform in equally spaced sub-windows',
                        default=5)
    parser.add_argument('--protName',
                        metavar='Protein',
                        dest='protName',
                        type=str,
                        help='Name to identify all protein chains in input topology.'
                        ' E.g. "Chain" if in topology there is Chain_A, Chain_B, and Chain_C.'
                        ' Default is Protein',
                        default='Protein')
    parser.add_argument('--ligName',
                        metavar='LIG',
                        dest='ligName',
                        type=str,
                        help='Name of the ligand in the input tropology. Default is LIG',
                        default='LIG')
    parser.add_argument('--cpu',
                        metavar='cpus',
                        dest='cpus',
                        type=int,
                        help='Number of cpus to be used.'
                        ' Default is 8.',
                        default=8)
    parser.add_argument('--lwindow',
                        metavar='lw',
                        dest='lw',
                        type=int,
                        help='Width of each sub-window in ps. Default: 500 ps',
                        default=500)
    parser.add_argument('--dt',
                        metavar='dt',
                        dest='dt',
                        type=int,
                        help='Time step in ps of the input trajectory. Default: 5 ps',
                        default=5)
    parser.add_argument('--keep',
                        dest='keep',
                        help='Keep intermediate files. Default not kept.',
                        default=False,
                        action='store_true')
    parser.add_argument('--decomp',
                        dest='decomp',
                        help='Perform per-residue energy decomposition. Default not done.',
                        default=False,
                        action='store_true')


    args, unknown = parser.parse_known_args()
    expected = ["--pdb","--top","--trj","--outDir","--ligName",
                "--protName","--cpu","--keep","--start","--end",
                "--windows","--lwindow","--dt","--decomp"]

    for cmd in unknown:
        if cmd not in expected:
            print('\n\n{0}\nUnknown command found in your command line: "{1}".\n{0}\n\n'.format('*'*50,cmd))
            sys.exit(1)

    return args

def check_directories(dirs):
    for dir in dirs:
        if not os.path.isdir(dir):
            print("\n\n{0}\nError! Provided directory '{1}' does not exist or is not a directory.\n{0}\n\n".format('*'*50,dir))
            sys.exit(2)

def check_files(files):
    for file in files:
        if not os.path.isfile(file):
            print("\n\n{0}\nError! Provided file '{1}' does not exist or is not a file.\n{0}\n\n".format('*'*50,file))
            sys.exit(3)

def CheckExit(result,logpath="log.log"):
    """
    This function checks the ouput code of a subprocess call to an
    utility and warns and exits if it is non-zero
    """
    if result.returncode != 0:
        print("\n\n{0}\nFATAL: Error encountered in {1} Check {2} for details!\n\n{0}\n\n".format('*'*50,result.args[1],logpath))
        sys.exit(4)


def mdpcreator(mdpdict,mdpfold):
with open(mdpfold+"simfile.mdp", 'w') as f:
    for key, value in mdpdict.items():
        f.write('%s=%s\n' % (key, value))

def g_grompp(MdpFile,,TopolFile,TprFile):

    process = subprocess.run(['gmx grompp','-c'',CoordFile,
                              '-t',TopolFile,
                              '-f',MdpFile,
                              '-o',TprFile],
                              stdout=glog,
                              stderr=glog)



def g_mdrun(TprFile=None,nsteps=None,gpu=False,cpi_file):
    if TprFile is None:
        print("\n\n provide a Tpr file")
        return 0
    if gpu:
    process = subprocess.run(['gmx mdrun','-s ',tpr,
                                  '-o',outname,
                                  '-nt',cpu,
                                  '-nsteps',nsteps,
                                  '-cpi',cpi_file,],
                                  stdout=glog,
                                  stderr=glog)

    else:
    process = subprocess.run(['gmx mdrun','-s ',tpr,
                              '-o',outname,
                              '-nt',cpu,
                              '-nsteps',nsteps,
                              '-cpi',cpi_file,
                              '-nb gpu'],
                              stdout=glog,
                              stderr=glog)
