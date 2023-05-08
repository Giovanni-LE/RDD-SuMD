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

    h="""
    This function checks if the directory exist or is not a directory
    """
    for dir in dirs:
        if not os.path.isdir(dir):
            print("\n\n{0}\nError! Provided directory '{1}' does not exist or is not a directory.\n{0}\n\n".format('*'*50,dir))
            sys.exit(2)
#faccio il check che il file esista
def check_files(files):

    h="""
    This function checks if the file exist or is not a file
    """
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
      
    h="""
GMX GROMPP
    Gmx Grompp  (the gromacs preprocessor) 
    reads a molecular topology file, checks the validity of the file,
    expands the topology 
    from a molecular description to an atomic description.


Options to specify input files:
            
-CoordFile [<.gro/.g96/…>] (conf.gro)
    Structure file: gro g96 pdb brk ent esp tpr

-CoordFile [<.gro/.g96/…>] (restraint.gro) (Optional)
    Structure file: gro g96 pdb brk ent esp tpr

-TopolFile [<.top>] (topol.top)
    Topology file

-MdpFile [<.mdp>] (grompp.mdp)
    grompp input file with MD parameters

-IndexFile [<.ndx>] (index.ndx) (Optional)
    Index file


Options to specify output files:

-TprFile [<.tpr>] (topol.tpr)
    Portable xdr run input file


Other options:

-maxwarn <int> (0)
    Number of allowed warnings during input processing. Not for normal use and may generate unstable systems 
"""
    
    cmd=['gmx' ,'grompp','-c',CoordFile,'-r',CoordFile,'-p',TopolFile,'-f',MdpFile,'-n'*(IndexFile != 0), str(IndexFile)*(IndexFile!=0),'-o',TprFile,'-maxwarn', str(maxwarning)]
    cmd=[el for el in cmd if el != '']
    if not os.path.exists('log'):
        os.mkdir('log')
    with open('log/grompp.log', 'w') as glog:
        process = subprocess.run(cmd,stdout=glog,stderr=glog)


       
def g_mdrun(OutputFile, xtcFile, log_File, edr_File,groFile, TprFile=None,gpu=False,append=False):

  h="""
GMX MDRUN
    GMX MDRUN gromacs command is used to perform molecular dynamics (MD) in the ) in gromacs simulations.
    MD is a molecular modelling simulation technique to study the behaviour of molecular systems under different conditions, such as such as temperature, pressure, density and composition

Options to specify input files:

-TprFile [<.tpr>] (topol.tpr)
    Portable xdr run input file

-cpi_file [<.cpt>] (state.cpt) (Optional)
     Checkpoint file


Options to specify output files:

-OutputFile [<.trr/.cpt/…>] (traj.trr)
     Full precision trajectory: trr cpt tng

-xtcFile [<.xtc/.tng>] (traj_comp.xtc) (Optional)
        Compressed trajectory (tng format or portable xdr format)

-cpo_file [<.cpt>] (state.cpt) (Optional)
        Checkpoint file

-edr_File [<.edr>] (ener.edr)
        Energy file

-log_File [<.log>] (md.log)
        Log file


Other options:

 -'gpu' <enum> (auto)
        Calculate non-bonded interactions on: auto, cpu, gpu

-[no]append (yes)
        Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names
 """
  if not os.path.exists('log'):
      os.mkdir('log')
  with open('log/mdrun.log', 'w') as glog:
    if TprFile is None:
      print("\n\n provide a Tpr file")
      return 0

    if append == True:
      append = '-append'

    if gpu:
      cmd = ['gmx', 'mdrun','-s',TprFile,'-o',OutputFile,'-c',groFile,'-x',xtcFile,'-g', log_File,'-e', edr_File,'-nb', 'gpu']
    else:
      cmd = ['gmx', 'mdrun','-s',TprFile,'-o',OutputFile,'-c',groFile,'-x',xtcFile,'-g', log_File,'-e', edr_File,'-nb', 'cpu']

    if append == '-append':
      cmd += [append]

    cmd=[el for el in cmd if el != '']
    process = subprocess.run(cmd,stdout=glog,stderr=glog)           


#lancio editconf
def g_editconf(f='conf.gro', o='out.gro',n=0,d='0',bt='cubic'):
    
    h="""
Gmx Editconf 
Gmx editconf converts generic structure format to .gro, .g96 or .pdb.

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

    h="""
GMX INSERT-MOLECULE
    Gmx insert molecule inserts -nmol copies of the system specified in the -ci input file.
        The insertions take place either into vacant space in the solute conformation given with -f, 
        or into an empty box given by -box. Specifying both -f and -box behaves like -f, 
        but places a new box around the solute before insertions. Any velocities present are discarded.

Options to specify input files:

-f [<.gro/.g96/…>] (protein.gro) (Optional)
    Existing configuration to insert into: gro g96 pdb brk ent esp tpr

-ci [<.gro/.g96/…>] (insert.gro)
    Configuration to insert: gro g96 pdb brk ent esp tpr

-ip [<.dat>] (positions.dat) (Optional)
    Predefined insertion trial positions


Options to specify output files:

-o [<.gro/.g96/…>] (out.gro)
    Output configuration after insertion: gro g96 pdb brk ent esp


Other options:

-nmol <int> (0)
    Number of extra molecules to insert
    """
    if not os.path.exists('log'):
        os.mkdir('log')         
    with open('log/insert-molecules.log', 'w') as glog:
        cmd=['gmx' ,'insert-molecules','-f',f,'-ci',ci,'-o',o,'-nmol',nmol,'-ip'*(ip != 0),ip*(ip != 0)]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)    

        
def g_solvate(cp="protein.gro",p="topol.top",o="out.gro"):

    h="""
GMX SOLVATE
    Gmx Solvate generate a box of solvent. Specify -cs and -box. Or specify -cs and -cp with a structure file with a box, but without atoms.


Options to specify input files:

-cp [<.gro/.g96/…>] (protein.gro) (Optional)
    Structure file: gro g96 pdb brk ent esp tpr


Options to specify input/output files:

-p [<.top>] (topol.top) (Optional)
    Topology file


Options to specify output files:

-o [<.gro/.g96/…>] (out.gro)
    Structure file: gro g96 pdb brk ent esp
    """

    if not os.path.exists('log'):
            os.mkdir('log')         
    with open('log/solvate.log', 'w') as glog:
        cmd=['gmx' ,'solvate','-cp',cp,'-p',p,'-o',o]
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)      


        
def g_genion(s="topol.tpr",p="topol.top",o="confout.gro",neutral="yes",conc=0,replace='System'):

    h="""
GMX GENION
    Gmx genion randomly replaces solvent molecules with monoatomic ions. 
    The group of solvent molecules should be continuous and all molecules should have the same number of atoms. 
    The user should add the ion molecules to the topology file or use the -p option to automatically modify the topology.

Options to specify input files:

-s [<.tpr>] (topol.tpr)
    Portable xdr run input file



Options to specify input/output files:

-p [<.top>] (topol.top) (Optional)
    Topology file


Options to specify output files:

-o [<.gro/.g96/…>] (out.gro)
    Structure file: gro g96 pdb brk ent esp



Other options:

-conc <real> (0)
    Specify salt concentration (mol/liter). 
    This will add sufficient ions to reach up to the specified concentration as computed from the volume of the cell in the input .tpr file.
    Overrides the -np and -nn options.

-[no]neutral (no)
    This option will add enough ions to neutralize the system. These ions are added on top of those specified with -np/-nn or -conc.
"""

    if not os.path.exists('log'):
            os.mkdir('log')         
    with open('log/genion.log', 'w') as glog:
        cmdecho=['echo','-e',f'{replace} \n']
        cmd=['gmx' ,'genion','-s',s,'-p',p,'-o',o,'-neutral',neutral,'-conc'*(conc !=0),str(conc)*(conc !=0)]
        cmd=[el for el in cmd if el != '']
        echo = subprocess.Popen(cmdecho, stdout=subprocess.PIPE)
        process = subprocess.run(cmd,stdout=glog,stderr=glog,stdin=echo.stdout)
            


def g_make_ndx(f="conf.gro",n=0,o="index.ndx",Sel="\n"):

    h="""
GMX MAKE NDX
    You ONLY have to use gmx make_ndx when you need SPECIAL index groups. 
    There is a default index group for the whole system, 9 default index groups for proteins, 
    and a default index group is generated for every other residue name.
    When no index file is supplied, also gmx make_ndx will generate the default groups. 
    With the index editor you can select on atom, residue and chain names and numbers.


Options to specify input files:

-f [<.gro/.g96/…>] (conf.gro) (Optional)
    Structure file: gro g96 pdb brk ent esp tpr

-n [<.ndx> […]] (index.ndx) (Optional)
    Index file


Options to specify output files:

-o [<.ndx>] (index.ndx)
    Index file
"""
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



def g_pdb2gmx(pdbFile,TopolFile, groFile,posre, ff,water):

    h="""
GMX pdb2gmx 
    reads a .pdb (or .gro) file, reads some database files, 
    adds hydrogens to the molecules and generates coordinates in GROMACS (GROMOS), 
    or optionally .pdb, format and a topology in GROMACS format. 
    These files can subsequently be processed to generate a run input file.
  


Options to specify input files:

-pdbFile [<.gro/.g96/…>] (conf.gro) (Optional)
    Structure file: gro g96 pdb brk ent esp tpr



Options to specify output files:

-groFile [<.gro/.g96/…>] (conf.gro)
    Structure file: gro g96 pdb brk ent esp

-TopolFile [<.top>] (topol.top)
    Topology file

-posre [<.itp>] (posre.itp)
    Include file for topology


Other options:

-ignh <string> (select)
    Force field, interactive by default. Use -h for information.

-water <enum> (select)
    Water model to use: select, none, spc, spce, tip3p, tip4p, tip5p, tips3p
"""
    if not os.path.exists('log'):
        os.mkdir('log')  

    with open('log/pdb2gmx.log', 'w') as glog:
        cmd=['gmx' ,'pdb2gmx','-f',pdbFile,'-p',TopolFile,'-o',groFile,'-i',posre ,'-water', water,'-ff', ff,'-ignh']
        cmd=[el for el in cmd if el != '']
        process = subprocess.run(cmd,stdout=glog,stderr=glog)