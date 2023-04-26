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


#controllo che sia presente gromacs sulla macchina

if shutil.which('gmx') is None:
    print('\n{0}\nThe program gmx is not present in the current environment\n{0}\n'.format('#'*80))
    sys.exit(10)

#definisco le opzioni per il parsing nel caso in cui il comando venga lanciato da terminale
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

#faccio il check che la directory esista
def check_directories(dirs):
    for dir in dirs:
        if not os.path.isdir(dir):
            print("\n\n{0}\nError! Provided directory '{1}' does not exist or is not a directory.\n{0}\n\n".format('*'*50,dir))
            sys.exit(2)
#faccio il check che il file esista
def check_files(files):
    for file in files:
        if not os.path.isfile(file):
            print("\n\n{0}\nError! Provided file '{1}' does not exist or is not a file.\n{0}\n\n".format('*'*50,file))
            sys.exit(3)

#controllo che il processo sia stato eseguito correttamente
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


#lancio grompp 
def g_grompp(MdpFile='grompp.mdp', TopolFile='topol.top', CoordFile='conf.gro', TprFile='topol.tpr', IndexFile=0, maxwarning=0):
    """
    Descrizione...
    """
    
    

    cmd=['gmx' ,'grompp','-c',CoordFile,'-r',CoordFile,'-p',TopolFile,'-f',MdpFile,'-n'*(IndexFile != 0), str(IndexFile)*(IndexFile!=0),'-o',TprFile,'-maxwarn', str(maxwarning)]
    cmd=[el for el in cmd if el != '']
    if not os.path.exists('log'):
        os.mkdir('log')
    with open('log/grompp.log', 'w') as glog:
        process = subprocess.run(cmd,stdout=glog,stderr=glog)

#lancio mdrun con o senza gpu
def g_mdrun(OutputFile, xtcFile,cpo_file, log_File, edr_File, TprFile=None,gpu=False,cpi_file=None):

    with open('md-sim/g_mdrun.log', 'w') as glog:
        if TprFile is None:
            print("\n\n provide a Tpr file")
            return 0

        if gpu:
            process = subprocess.run(['gmx', 'mdrun',
                                        '-s',TprFile,
                                        '-o',OutputFile,
                                        '-x',xtcFile,
                                        '-g', log_File,
                                        '-e', edr_File,
                                        '-cpo',cpo_file,
                                        '-nb', 'gpu'],
                                        stdout=glog,
                                        stderr=glog)


#lancio editconf
def g_editconf(f='conf.gro', o='out.gro',n=0,d='0',bt='cubic'):
    
    h="""
gmx editconf converts generic structure format to .gro, .g96 or .pdb.

Options to specify input files:

 -f      [<.gro/.g96/...>]  (conf.gro)
           Structure file: gro g96 pdb brk ent esp tpr
 -n      [<.ndx>]           (index.ndx)      (Opt.)
           Index file
 -bf     [<.dat>]           (bfact.dat)      (Opt.)
           Generic data file

Options to specify output files:

 -o      [<.gro/.g96/...>]  (out.gro)        (Opt.)
           Structure file: gro g96 pdb brk ent esp
 -mead   [<.pqr>]           (mead.pqr)       (Opt.)
           Coordinate file for MEAD

Other options:
 -d      <real>             (0)
           Distance between the solute and the box

    """
    if not os.path.exists('log'):
        os.mkdir('log')         
    with open('log/editconf.log', 'w') as glog:
        cmd=['gmx' ,'editconf','-f',f,'-o',o,'-d',d,'-n'*(n != 0),str(n)*(n != 0),'-bt',bt]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)
        #process = subprocess.run(cmd,stdout=glog,stderr=glog)


             
def g_insert_molecule(f='protein.gro',ci='insert.gro',o='out.gro',nmol='1',ip=0):
    if not os.path.exists('log'):
        os.mkdir('log')         
    with open('log/insert-molecules.log', 'w') as glog:
        cmd=['gmx' ,'insert-molecules','-f',f,'-ci',ci,'-o',o,'-nmol',nmol,'-ip'*(ip != 0),ip*(ip != 0)]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)    
        
def g_solvate(cp="protein.gro",p="topol.top",o="out.gro"):
    if not os.path.exists('log'):
            os.mkdir('log')         
    with open('log/solvate.log', 'w') as glog:
        cmd=['gmx' ,'solvate','-cp',cp,'-p',p,'-o',o]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)      
        
def g_genion(s="topol.tpr",p="topol.top",o="confout.gro",neutral="yes",conc=0,replace='System'):
    if not os.path.exists('log'):
            os.mkdir('log')         
    with open('log/genion.log', 'w') as glog:
        cmdecho=['echo','-e',f'{replace} \n']
        cmd=['gmx' ,'genion','-s',s,'-p',p,'-o',o,'-neutral',neutral,'-conc'*(conc !=0),str(conc)*(conc !=0)]
        cmd=[el for el in cmd if el != '']
        echo = subprocess.Popen(cmdecho, stdout=subprocess.PIPE)
        process = subprocess.run(cmd,stdout=glog,stderr=glog,stdin=echo.stdout)
            
def g_make_ndx(f="conf.gro",n=0,o="index.ndx",Sel="\n"):
    if not os.path.exists('log'):
            os.mkdir('log')         
    with open('log/make-ndx.log', 'w') as glog:
        SELE=' \n'.join(Sel)
        SELE=SELE+"\n q"
        cmdecho=['echo','-e',f'{SELE} \n']
        cmd=['gmx' ,'make_ndx','-f',f,'-o',o,'-n'*(n != 0),str(n)*(n != 0)]
        cmd=[el for el in cmd if el != '']
        echo = subprocess.Popen(cmdecho, stdout=subprocess.PIPE)
        process = subprocess.run(cmd,stdout=glog,stderr=glog,stdin=echo.stdout)  

def g_energy(f,o):
    if not os.path.exists('log'):
        os.mkdir('log')  

    with open('log/energy.log', 'w') as glog:
        cmd=['gmx' ,'energy','-f',f,'-o',o]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)
