import sys
import argparse
import warnings
import os
import fileinput
import glob
import shutil
import subprocess
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import Gromacs_Tools #Import dello script Gromacs_Tools.py contenente le funzioni per eseguire i comandi di Gromacs
from MDAnalysis import transformations
import pyfiglet
import GmxMdp
import time

def parse_options():
  parser = argparse.ArgumentParser(description='''SuMD''')
  parser.add_argument('--pdb',
                        metavar='complex.pdb',
                        dest='str',
                        type=str,
                        help='Input protein-ligand pdb file.',
                        required=True)
  for cmd in unknown:
    if cmd not in expected:
      print('\n\n{0}\nUnknown command found in your command line: "{1}".\n{0}\n\n'.format('*'*50,cmd))
      sys.exit(1)
  return args    

class GmxMdp:
    def __init__(self) -> None:
        self.mdp_params=dict()
        self.mdpfile_out=False
        self.mdpfile_in=False
    def read_mdp(self,mdpfile_in=False):
    #This function read an mdp file
        if mdpfile_in:
            "Give a correct mdpfilename"
            sys.exit(-1)
        with open(mdpfile_in,'r') as f:
            #cicle in the mdp file with trivial check for the correctness of the file
            #check if a line starts with a space
            for lines in f:
                line_nospace=lines.replace.replace(" ","")
                if line_nospace[0]!=';' and line_nospace[0].isalpha():
                    #TODO Pensare una gestione più intelligente dei parametri
                    mdpparam=lines.strip().split("=")
                    mdpparam[0]=mdpparam[0].replace.replace(" ","")
                    self.mdp_params[mdpparam[0]]=mdpparam[1]
    def write_mdp(self,mdpfile_out=False):
        if mdpfile_out:
            "Give a correct mdpfilename"
            sys.exit(-1)
        self.mdpfile_out=mdpfile_out
        with open(mdpfile_out,'w') as f:
            for key in self.mdp_params.keys():
                f.write(f"{key} = {self.mdp_params[key]}\n")
    def list_mdp_params(self):
        for key in self.mdp_params:
            print(f"{key} = {self.mdp_params[key]}\n")


class SuMD:
  def __init__(self,PDBPROT=None,POCKET=None,PDBLIG=None,TopFOLD=None,distance_from_edge=0.5,distance_from_pocket=0.4):
    self.sumd_param={'m':0, # threshold of the angular coefficient values of the straight line that interpolates the saved points
      't_max':31, # maximum threshold of bankruptcy attempts granted in the preliminary run
      'counter1t':17, # maximum threshold of failed attempts granted in the SuMD run
      'ct02':19, # threshold relating to the times in which the distance between the com is between 0 and 2 Å during the final phase of unsupervised MD
      'ct25':19, # threshold relating to the times in which the distance between the com is between 2 and 5 Å during the final phase of unsupervised MD
      'ct59':19, # threshold relating to the times in which the distance between the com is between 5 and 9 Å during the final phase of unsupervised MD
      'distanza_suMD_FS':9, #Angstrom - rappresenta la distanza alla quale vogliamo che si entri nel final step
      'n_sample':5} 
    self.md_param={"coulombtype":"pme",
                   "coulomb-modifier":"Potential-shift-Verlet",
                   "rcoulomb-switch":"0.9",
                   "rcoulomb":"1",
                   "epsilon-r":"1",
                   "epsilon-rf":"1",
                   "vdwtype":"cutoff",
                   "vdw-modifier":"Potential-shift-Verlet",
                   "rvdw-switch":"0.9","rvdw":"1",
                   "constraints":"h-bonds",
                   "constraint_algorithm":"lincs"}
    if PDBPROT is None or POCKET is None or PDBLIG is None or TopFOLD is None:
      print("path to the protein pdb file:", PDBPROT)
      print("protein residues are:", POCKET)
      print("path to the lig pdb file:", PDBLIG)
      print("path to the topology folder:", TopFOLD)
      sys.exit(0,"Error one of the input file is not valid")
    
    print("path to the protein pdb file:", PDBPROT)
    print("protein residues are:", POCKET)
    print("path to the lig pdb file:", PDBLIG)
    print("path to the topology folder:", TopFOLD)    
    self.PDBPROT=PDBPROT
    self.POCKET=POCKET
    self.PDBLIG=PDBLIG
    self.TopFold=TopFOLD
    self.distance_from_edge=str(distance_from_edge)
    self.distance_from_pocket=distance_from_pocket*10
  
  def Linear_Fitting(self,PDB,XTC,n,POCKET=False,LIG="resname LIG"): #Funzione per il fitting lineare
    '''This function calculte the angular coefficent of the distance between POCKET and LIGAND over time
      PDB : str
      the path of the pdb file
      XTC :str 
      the path of the trajectory file
      n : int
      ?????
      POCKET : str
      the selection of the binding pocket
      LIG: str
      the selection of the ligand'''  
    #Check dei file
    if os.stat(PDB).st_size == 0: #check PDB file
      print('Il file PDB è vuoto')
      return -1
    
    if os.stat(XTC).st_size == 0: #check XTC file
      print('Il file XTC è vuoto')
      return -1

    u = mda.Universe(PDB,XTC) #universe creation
    workflow = [transformations.unwrap(u.atoms)]
    u.trajectory.add_transformations(*workflow)
    pocket = u.select_atoms(POCKET)
    LIG = u.select_atoms(LIG)
      
    n_frame = u.trajectory.n_frames
    #print(f"Numero di frame di u.trajectory: {n_frame}")

    step = n_frame/(int(n)-1) #calculate how many steps of the trajectory will be analyzed (used as integer)
    #print(step) #debug
    d = np.empty(len(u.trajectory[:n_frame:int(step)]), dtype=float) #distance vector

    p = 0 #aux variable for vector indexing

    for ts in u.trajectory[:n_frame:int(step)]: #trajectory analysis cycle
      d[p] = float(distances.distance_array(np.array(pocket.center_of_mass()),np.array(LIG.center_of_geometry())))
      p += 1

    x = np.linspace(0,n_frame-1,int(n), dtype=float) #x-axis time reference

      #print(f"Lunghezza di x: {len(x)} - x: {x}")
      #print(f"Lunghezza di y: {len(d.astype(np.float32))} - y: {d.astype(np.float32)}")

    m1 = np.polyfit(x,d.astype(np.float32),1) #linear fitting of distance vector
    return(m1[0],d)

  def g_mdrun(self,OutputFile, xtcFile, log_File, edr_File,groFile, TprFile=None,gpu=False,append=False):

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

#TODO tirare fuori la selezione del ligando dalla funzione

  def SuMdStep(self,MF,append):
      und_sim=self.und_sim
      Gromacs_Tools.g_grompp(MdpFile=MF, TopolFile= os.path.join(self.outdir,"topol","topol.top"), CoordFile= os.path.join(self.act_sim,"md.gro"), TprFile= os.path.join(self.und_sim,"md.tpr"), IndexFile= os.path.join(self.minim,"index.ndx"),maxwarning='2')   
      if not os.path.exists(os.path.join(self.und_sim,"md.tpr")):
          sys.exit('error during the creation of tpr file (grompp)')
          
      g_mdrun(OutputFile= os.path.join(und_sim,"md.trr"),xtcFile= os.path.join(und_sim,"md.xtc"), log_File= os.path.join(und_sim,"md.log"),edr_File=os.path.join(und_sim,"md.edr"),groFile=os.path.join(und_sim,"md.gro"), TprFile=os.path.join(und_sim,"md.tpr"),gpu=self.gpu,append=append)
      if not os.path.exists(os.path.join(und_sim,"md.xtc")):
          sys.exit('error during the creation of xtc file (mdrun)')
      
      Gromacs_Tools.g_editconf(f=os.path.join(und_sim,"md.tpr"), o=os.path.join(und_sim,"md.pdb")) # per avere il file pdb
      if not os.path.exists(os.path.join(und_sim,"md.pdb")):
          sys.exit('error during the creation of pdb file (editconf)')

      m_primo,dcm_vector = Linear_Fitting(PDB=os.path.join(und_sim,"md.tpr"),XTC = os.path.join(und_sim,"md.xtc"),n=n,POCKET='resid 267 or resid 252 or resid 264 or resid 256 or resid 174')
      return m_primo,dcm_vector


  def BuildGeometry(self,PDBPROT=None,PDBLIG=None,TopFOLD=None,newbox=True,d='5',mindist=100,solvate=False,Saltconc="0.15"):
      '''
      BuildGeometry
      Build the system geometry from a protein input file and a ligand input file, place randomly the ligand at a given
      position  
      '''
      
      try:
          from Gromacs_Tools import g_editconf
          from Gromacs_Tools import g_insert_molecule
          from Gromacs_Tools import g_solvate
          from Gromacs_Tools import g_make_ndx
          from Gromacs_Tools import g_genion
          from Gromacs_Tools import g_grompp
      except:
          sys.exit("Error in module loading")
      #Check dei file
      if os.stat(PDBPROT).st_size == 0: #check PDBPROT file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBLIG.split('/')[-1] )   
    
      if os.stat(PDBLIG).st_size == 0: #check PDBLIG file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBPROT.split('/')[-1])  
      
      if os.stat(TopFOLD).st_size == 0: #check PDBLIG file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBPROT.split('/'))  
      
          outfile=PDBPROT
      if newbox:
          outfile=os.path.join(outdir,'tmp.box.pdb')#define the path of the file
          g_editconf(f=PDBPROT,d=d,o=outfile)#create the new configuration file
      prot=mda.Universe(outfile)#load the new configuration in order to place the ligand in a random position above mindist
      dist=0 #flag to go out
      insert_try=0 #flag to avoid infinite loop
      
      while dist<mindist:#check to insert the ligand
          cord=np.random.rand(3)*np.min(prot.dimensions)#create a random position
          dist=distances.distance_array(cord,prot.atoms.positions).min()#chek if the COM of the new position is above mindist
          insert_try+=1 #increase try variable
          if insert_try>1000:#check if i Tried to inser the molecule more than 1000 times to avoid infinte loop
              sys.exit(f"Try to inser the molecule more than {insert_try} times, something wrong in the parameters") #exit in case of many insertion with an error
      print(f"moleculer will be inserted in {cord/10}, after {insert_try} attemps")#inform the user about the position of the placed molecule
      posfile=os.path.join(outdir,"pos.dat")# define the path the position file
      
      with open(posfile,'w') as f:
          ele=[f"{numb/10:>.2f}" for numb in cord.tolist()]#traslate the cordinate to a string /10 to convert Amstrong to nm, 
          f.write(' '.join(ele))
          f.write('\n')
      outsys=os.path.join(outdir,"sys.pdb")#define system path
      g_insert_molecule(f=outfile,o=outsys,ci=PDBLIG,ip=posfile)#insert the molecule in the random position defined before  
      TopFold=os.path.join(outdir,"topol")#define the position of the topology file

      if os.path.exists(TopFold): #chekc if the file exist
          shutil.rmtree(TopFold)

      shutil.copytree(TopFOLD,TopFold) #copy the topology folder in order to don't alterate the input
      
      if solvate:
          print(f"I will solvate the system, and neutralyze it with {Saltconc} mol/L")
          solvfile=os.path.join(outdir,"solv.pdb")# define the solvent file output
          topfile=os.path.join(TopFold,"topol.top")#define the topology outputpaht
          g_solvate(cp=outsys,p=topfile,o=solvfile)#add water to the system and modify the topology file
          
          with open("tmp.mdp","w") as f: #create a dummy mdp file
              f.write("\n")
          ionstprout=os.path.join(outdir,"ions.tpr")#define the output of the tpr file to add ions
          g_grompp(MdpFile="tmp.mdp",TopolFile=topfile,CoordFile=solvfile,TprFile=ionstprout)#create the tpr file to add ions
          os.remove("tmp.mdp")#remove the mdp file
          os.remove(outsys)# remove the old system file
          g_genion(s=ionstprout,neutral='yes',conc=Saltconc,replace='SOL',o=outsys,p=topfile) #add ions
          print(f"The system has been solvated, placed in a cubic box with size {prot.dimensions[0:3]},\nsystem path {outsys}")

      else:
          print(f"The system has not been solvated, placed in a cubic box with size {str(prot.dimensions[0])},\nsystem path {outsys}")
      return 0

  def RandSphere(r=1,cent=[0,0,0]):
      ''' RandSphere return a random point on a sphere of radius r, centered in cent
      r: float the radius of the sphere in nm
      cent: the center of the sphere'''
      tetha,phy=np.random.rand(2)
      tetha=tetha*2*np.pi
      phy=phy*np.pi
      zi=r*np.sin(phy)+cent[2]
      xy=r*np.cos(phy)
      xi=xy*np.cos(tetha)+cent[0]
      yi=xy*np.sin(tetha)+cent[1]
      return np.array([xi,yi,zi])

  def BuildGeometry(PDBPROT=None,PDBLIG=None,TopFOLD=None,newbox=True,d='5',dist=100,solvate=False,Saltconc="0.15",pocket="all"):
      '''
      BuildGeometry
      Build the system geometry from a protein input file and a ligand input file, place randomly the ligand at a given  
      '''
      
      try:
          from Gromacs_Tools import g_editconf
          from Gromacs_Tools import g_insert_molecule
          from Gromacs_Tools import g_solvate
          from Gromacs_Tools import g_make_ndx
          from Gromacs_Tools import g_genion
          from Gromacs_Tools import g_grompp
      except:
          sys.exit("Error in module loading")
      #Check dei file
      if os.stat(PDBPROT).st_size == 0: #check PDBPROT file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBLIG.split('/')[-1] )   
    
      if os.stat(PDBLIG).st_size == 0: #check PDBLIG file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBPROT.split('/')[-1])  
      
      if os.stat(TopFOLD).st_size == 0: #check PDBLIG file
          raise ValueError('The file %s if empty, file doesn\'t exixst'%PDBPROT.split('/'))  
      outfile=PDBPROT
      if newbox:
          outfile=os.path.join(outdir,'tmp.box.pdb')#define the path of the file
          g_editconf(f=PDBPROT,d=d,o=outfile)#create the new configuration file
      #Load protein and Ligand
      
      prot=mda.Universe(outfile)#load the new configuration in order to place the ligand in a random position above mindist
      LIG=mda.Universe(PDBLIG)
      rg=LIG.atoms.radius_of_gyration()*1.1#10% security coefficent
      insert_try=0 #flag to avoid infinite loop
      mindist=0 #flag to go out from simulation
      while rg>=mindist:
          cent=prot.select_atoms(POCKET).center_of_mass()
          cord=RandSphere(r=dist,cent=cent) 
          mindist=distances.distance_array(cord,prot.atoms.positions).min()#chek if the COM of the new position is above mindist
          insert_try+=1 #increase try variable
          print(f"{cord}")
          if insert_try>1000:#check if i Tried to inser the molecule more than 1000 times to avoid infinte loop
              sys.exit(f"Try to inser the molecule more than {insert_try} times, something wrong in the parameters") #exit in case of many insertion with an error
      '''
      while dist<mindist:#check to insert the ligand
          cord=np.random.rand(3)*np.min(prot.dimensions)#create a random position
          dist=distances.distance_array(cord,prot.atoms.positions).min()#chek if the COM of the new position is above mindist
          insert_try+=1 #increase try variable
          if insert_try>1000:#check if i Tried to inser the molecule more than 1000 times to avoid infinte loop
              sys.exit(f"Try to inser the molecule more than {insert_try} times, something wrong in the parameters") #exit in case of many insertion with an error
      '''
      print(f"moleculer will be inserted in {cord/10}, after {insert_try} attemps")#inform the user about the position of the placed molecule
      posfile=os.path.join(outdir,"pos.dat")# define the path the position file
      
      with open(posfile,'w') as f:
          ele=[f"{numb/10:>.2f}" for numb in cord.tolist()]#traslate the cordinate to a string /10 to convert Amstrong to nm, 
          f.write(' '.join(ele))
          f.write('\n')
      outsys=os.path.join(outdir,"sys.pdb")#define system path
      g_insert_molecule(f=outfile,o=outsys,ci=PDBLIG,ip=posfile)#insert the molecule in the random position defined before  
      TopFold=os.path.join(outdir,"topol")#define the position of the topology file

      if os.path.exists(TopFold): #chekc if the file exist
          shutil.rmtree(TopFold)

      shutil.copytree(TopFOLD,TopFold) #copy the topology folder in order to don't alterate the input
      
      if solvate:
          print(f"I will solvate the system, and neutralyze it with {Saltconc} mol/L")
          solvfile=os.path.join(outdir,"solv.pdb")# define the solvent file output
          topfile=os.path.join(TopFold,"topol.top")#define the topology outputpaht
          g_solvate(cp=outsys,p=topfile,o=solvfile)#add water to the system and modify the topology file
          
          with open("tmp.mdp","w") as f: #create a dummy mdp file
              f.write("\n")
          ionstprout=os.path.join(outdir,"ions.tpr")#define the output of the tpr file to add ions
          g_grompp(MdpFile="tmp.mdp",TopolFile=topfile,CoordFile=solvfile,TprFile=ionstprout)#create the tpr file to add ions
          os.remove("tmp.mdp")#remove the mdp file
          os.remove(outsys)# remove the old system file
          g_genion(s=ionstprout,neutral='yes',conc=Saltconc,replace='SOL',o=outsys,p=topfile) #add ions
          print(f"The system has been solvated, placed in a cubic box with size {prot.dimensions[0:2]},\nsystem path {outsys}")

      else:
          print(f"The system has not been solvated, placed in a cubic box with size {str(prot.dimensions[0])},\nsystem path {outsys}")
      return 0

  def minimization_equilibration(MF_m,MF_e):
      '''
      this function make the minimization and equilibration process
      MF_m mdp file for the minimization
      MF_e mdp file for the equilibration phase
      '''
      print('\nMinimization ...', end='')
      Gromacs_Tools.g_grompp(MdpFile=MF_m, TopolFile=os.path.join(outdir,"topol","topol.top"), CoordFile= os.path.join(outdir,"sys.pdb"), TprFile= os.path.join(minim,"min.tpr"),maxwarning='2')
      if not os.path.exists(os.path.join(minim,"min.tpr")):
          sys.exit("\nError during the minimization step - grompp doesn't work")
      
      g_mdrun(OutputFile= os.path.join(minim,"min.trr"),xtcFile='', log_File= os.path.join(minim,"min.log") ,edr_File= os.path.join(minim,"min.edr") ,groFile= os.path.join(minim,"min.gro"), TprFile= os.path.join(minim,"min.tpr"),gpu=gpu)
      if os.path.exists(os.path.join(minim,"min.trr")):
          print(' minimization completed')
      else:
          sys.exit("Error during the minimization step - mdrun doesn't work")

      Gromacs_Tools.g_make_ndx(f=os.path.join(minim,"min.gro") , o= os.path.join(minim,"index.ndx"),Sel=["keep 0","rSOL|rLIG|rNA|rCL","name 1 solvent","!1","name 2 solute","\n"])
    # Gromacs_Tools.g_make_ndx(f=os.path.join(minim,"min.gro") , o= os.path.join(minim,"index.ndx"))

      if not os.path.exists(os.path.join(minim,"index.ndx")):
          sys.exit('Error during the creation of index file')


      print('Equilibration (NVT) ...', end='')
      Gromacs_Tools.g_grompp(MdpFile=MF_e, TopolFile=os.path.join(outdir,"topol","topol.top"), CoordFile= os.path.join(minim,"min.gro"), IndexFile= os.path.join(minim,"index.ndx") , TprFile= os.path.join(equil,"md.tpr"),maxwarning='2')
      if not os.path.exists(os.path.join(equil,"md.tpr")):
          sys.exit("\nError during the equilibration step - grompp doesn't work")

      g_mdrun(OutputFile= os.path.join(equil,"md.trr"),xtcFile='', log_File= os.path.join(equil,"md.log"), edr_File= os.path.join(equil,"md.edr"), groFile=os.path.join(equil,"md.gro"), TprFile=os.path.join(equil,"md.tpr"),gpu=gpu)
      if os.path.exists(os.path.join(equil,"md.trr")):
          print(' Equilibration completed')
      else:
          sys.exit("Error during the equilibration step - mdrun doesn't work")
      return 0

  def CreateWorkingFold(self):
    self.workfold="sumd"#define workfold folder --> folder che conterrà ogni operazione dell'algoritmo
    if not os.path.exists(self.workfold): #chekc if the workfold exists
      os.mkdir(self.workfold) #in case doesn't exist create the directory

    self.fal_sim=os.path.join(self.workfold,"failure_simulation")# define the failure_simulation folder --> folder che conterrà i file fallimentari
    if not os.path.exists(self.fal_sim): #chekc if the fal_sim exists
      os.mkdir(self.fal_sim) #in case doesn't exist create the directory

    self.act_sim=os.path.join(self.workfold,"actual_simulation")# define the actual_simulation folder --> folder che conterrà i file in input alla simulazione md i-esima
    if not os.path.exists(self.act_sim): #chekc if the act_sim exists
      os.mkdir(self.act_sim) #in case doesn't exist create the directory

    self.und_sim=os.path.join(self.workfold,"undefined_simulation") #define und_sim folder --> folder che conterrà i file in output alla simulazione md i-esima
    if not os.path.exists(self.und_sim): #chekc if the und_sim exists
      os.mkdir(self.und_sim) #in case doesn't exist create the directory

    self.corr_sim=os.path.join(self.workfold,"correct_simulation") #define und_sim folder --> folder che conterrà i file in output alla simulazione md i-esima che soddisfano la condizione m>0
    if not os.path.exists(self.corr_sim): #chekc if the corr_sim exists
      os.mkdir(self.corr_sim) #in case doesn't exist create the directory

    self.outdir= os.path.join(self.workfold,"sysprep") #define output folder --> folder che conterrà i file in output alla minimizzazione ed equilibrazione dopo il cambio della geometria
    if not os.path.exists(self.outdir): #chekc if the outputfolder exists
      os.mkdir(self.outdir) #in case doesn't exist create the directory

    self.minim= os.path.join(self.outdir,"minimization") #define minimization folder --> folder che conterrà l'output della minimizzazione 
    if not os.path.exists(self.minim): #chekc folder exists
      os.mkdir(self.minim) #in case doesn't exist create the directory

    self.equil= os.path.join(self.outdir,"equilibration") #define equilibration folder --> folder che conterrà l'output dell'equalizzazione NVT 
    if not os.path.exists(self.equil): #chekc folder exists
      os.mkdir(self.equil) #in case doesn't exist create the directory
    self.mdp_fold= os.path.join(self.outdir,"mdp") #define mdp folder --> folder che conterrà gli mdp
    if not os.path.exists(self.mdp_fold): #chekc folder exists
      os.mkdir(self.mdp_fold) #in case doesn't exist create the directory
  
  def EmMdpFile(self,emMdpinput=False):
    '''Questa funzione crea il file di input della minimizzazione
    '''
    self.EmMdp=GmxMdp()
    defaultEmparams={"integrator":"steep",
                     "nsteps":"10000",
                     "emtol":"1000",
                     "emstep":"0.01"}
    self.EmMdp.mdp_params.update(self.md_param)
    self.EmMdp.mdp_params.update(defaultEmparams)
    if emMdpinput:
      self.emMdp.read_mdp(emMdpinput)
    self.EmMdp.write_mdp(self.mdp_fold+"01-em.mdp")

  def NvtMdpFile(self,nvtMdpInput=False):
    '''Questa funzione crea il file di input della minimizzazione
    '''
    self.NvtMdp=GmxMdp()
    defaultNvtparams={"define":"-DPOSRES",
                      "integrator":"md",
                      'tinit':'0','dt':'0.002',
                      'nsteps':'300000',
                      'tcoupl':'V-rescale',
                      'nsttcouple':'-1',
                      'nh-chain-length':'10',
                      'tc-grps':'System',
                      'tau_t':'0.1',
                      'ref_t':'300',
                      'gen_vel':'yes',
                      'gen_temp':'300',
                      'gen_seed':'-1'}
    self.NvtMdp.mdp_params.update(self.md_param)
    self.NvtMdp.mdp_params.update(defaultNvtparams)
    if nvtMdpInput:
      self.NvtMdp.read_mdp(nvtMdpInput)
    self.NvtMdp.write_mdp(self.mdp_fold+"02-nvt.mdp")
  
  
  def NptMdpFile(self,nptMdpinput=False):
    '''Questa funzione crea il file di input della Npt di equilibrazione
    '''
    self.NptMdp=GmxMdp()
    defaultNptparams={"define":"-DPOSRES",
                      "integrator":"md",
                      'tinit':'0',
                      'dt':'0.002',
                      'nsteps':'300000',
                      'tcoupl':'V-rescale',
                      'nsttcouple':'-1',
                      'nh-chain-length':'10',
                      'tc-grps':'System',
                      'tau_t':'0.1',
                      'ref_t':'300',
                      'gen_vel':'no',
                      'pcoupl':'Parrinello-Rahman',
                      'pcoupltype':'isotropic',
                      'nstpcouple':'-1',
                      'tau_p':'5.0',
                      'compressibility':'4.5e-5',
                      'ref_p':'1.0'}
    self.NptMdp.mdp_params.update(self.md_param)
    self.NptMdp.mdp_params.update(defaultNptparams)
    if nptMdpInput:
      self.NptMdp.read_mdp(nptMdpInput)
    self.NptMdp.write_mdp(self.mdp_fold+"03-npt.mdp")

  def MdMdpFile(self,MdMdpInput=False):
    '''Questa funzione crea il file di input della minimizzazione
    '''
    self.MdMdp=GmxMdp()
    defaultMdparams={"integrator":"md",
                      'tinit':'0',
                      'dt':'0.002',
                      'nsteps':'300000',
                      'tcoupl':'V-rescale',
                      'nsttcouple':'-1',
                      'nh-chain-length':'10',
                      'tc-grps':'System',
                      'tau_t':'0.1',
                      'ref_t':'300',
                      'gen_vel':'no',
                      'pcoupl':'Parrinello-Rahman',
                      'pcoupltype':'isotropic',
                      'nstpcouple':'-1',
                      'tau_p':'5.0',
                      'compressibility':'4.5e-5',
                      'ref_p':'1.0',
                      'nstxout':'50000',
                      'nstvout':'50000',
                      'nstfout':'50000',
                      'nstlog':'500',
                      'nstenergy':'500'}
    self.MdMdp.mdp_params.update(self.md_param)
    self.MdMdp.mdp_params.update(defaultMdparams)
    if MdMdpInput:
      self.NptMdp.read_mdp(MdMdpInput)
    self.MdMdp.write_mdp(self.mdp_fold+"04-md.mdp")


  def run(self):
    #SuMD params inizialization
    print("The SuMD simulation will start with the following Parameters")
    for param_key in self.sumd_param.keys():
        print(f"|{param_key}\t{self.sumd_param[param_key]}|")
    
    '''    self.sumd_param={'m':0, # threshold of the angular coefficient values of the straight line that interpolates the saved points
      't_max':31, # maximum threshold of bankruptcy attempts granted in the preliminary run
      'counter1t':17, # maximum threshold of failed attempts granted in the SuMD run
      'ct02':19, # threshold relating to the times in which the distance between the com is between 0 and 2 Å during the final phase of unsupervised MD
      'ct25':19, # threshold relating to the times in which the distance between the com is between 2 and 5 Å during the final phase of unsupervised MD
      'ct59':19,
      'distanza_suMD_FS':9, #Angstrom - rappresenta la distanza alla quale vogliamo che si entri nel final step
      'n':5} # threshold relating to the times in which the distance between the com is between 5 and 9 Å during the final phase of unsupervised MD
    '''
    #Inizializzazione

    m_primo=0 # value of the angular coefficient of the line which is updated and compared with m each short MD simulation
    t_step=1 # counter that is increased by 1 and compared with t_max every time there is a failed step in the preliminary run
    counter1=0 # counter that is incremented by 1 and compared with counter1t every time there is a failed step in the SuMD run
    c_step=0 # counters that allow us to determine when a SuMD run step is productive
    k_step=0 # counters that allow us to determine when a SuMD run step is productive
    d_cm_vector=np.zeros(self.sumd_param['n'])# vector in which the distances of the com taken in the 5 steps are entered within each short MD simulation during the SuMD run
    d_cm_out=0 # final distance value dcm stored at the end of the short MD simulations that occur once supervision is stopped
    counter02=0 # counter that is updated and compared with ct02 every time that  dcm is between 0 and 2 Å
    counter25=0 # counter that is updated and compared with ct25 every time that dcm is between 2 and 5 Å
    counter59=0 # counter that is updated and compared with ct59 every time that dcm is between 5 and 9 Å
    prima_volta = 1
    flag = 1
    corrette = 1
    fallimento = 1 #tiene conto delle volte in cui è stato necessario ricominciare la simulazione
    

    #Definition folders
    act_sim=self.act_sim
    und_sim=self.und_sim
    corr_sim=self.corr_sim
    fal_sim=self.fal_sim
    mdp_fold=self.mdp_fold
    
    
    
    
    
    
    
    
    
    
    
    
    while flag == 1:
      if prima_volta == 1: #Se è la prima volta che si esegue lo script, allora è necessario generare un ambiente 
        print("Creation of the box containing protein and ligand\n")
        BuildGeometry(self.PDBPROT,self.PDBLIG,self.TopFOLD,solvate=True,d=self.distance_from_edge,dist=self.distance_from_pocket,pocket=self.POCKET)
        
        minimization_equilibration(MF_m = self.EmMdp.mdpfile_out, MF_e = self.NvtMdp.mdpfile_out)
        shutil.rmtree(act_sim), shutil.copytree(equil, act_sim) #Copio l'output dell'equalizzazione nella cartella actual_simulation
        prima_volta = 0
        cambiare_velocità = 0 # la velocità degli atomi della prima simulazione non deve essere cambiata 
      
      avviso_preliminary_run = True
      while m_primo >= 0: #PRELIMINARY RUN
        if avviso_preliminary_run:
          print("\nsuMD: The algorithm has ENTERED the 'PRELIMINARY RUN' step\n")
          avviso_preliminary_run = False

        if cambiare_velocità == 0: #Eseguire MD senza cambiare le velocità
          print(f"MD simulation #{corrette} ...", end='')
          self.MdMdp.mdp_params['gen-vel']='no'
          self.MdMdp.mdp_params['continuation']='yes'
          self.MdMdp.write_mdp(mdp_fold+'04-md.mdp')
          m_primo,dcm_vector = SuMdStep(self.sumd_param['n'],MF=self.MdMdp.mdpfile_out,append=True)
        else: #Eseguire MD cambiando le velocità
          self.MdMdp.mdp_params['gen-vel']='yes'
          self.MdMdp.mdp_params['continuation']='no'
          self.MdMdp.write_mdp(mdp_fold+'04-md.mdp')
          print(f"MD simulation #{corrette} - change of velocity ...", end='')
          m_primo,dcm_vector = SuMdStep(self.sumd_param['n'],MF=self.MdMdp.mdpfile_out,append=True)
        if m_primo > 0: #Non produttivo
          cambiare_velocità = 1
          print(f" The step was not productive (Attempt {t_step})")
          t_step += 1
          shutil.rmtree(und_sim), os.mkdir(und_sim) #Cancello la cartella dei risultati
        else: #m_primo < 0, il preliminary run è terminato!
          cambiare_velocità = 0
          print(f" The step was productive (Distance of the last step: {dcm_vector[-1]:.2f} Anstrong)")
      
          shutil.copytree(und_sim,os.path.join(corr_sim,f"sim{corrette}")) #In modo da raccogliere tutte le simulazioni corrette
          shutil.rmtree(act_sim), shutil.copytree(und_sim,act_sim) #In modo che alla prossima simulazione md si riparta dall'ultimo file .gro (le ultime velocità)
          shutil.rmtree(und_sim), os.mkdir(und_sim) #elimino cartella dei risultati
          corrette += 1
          print('\nsuMD: The algorithm has EXITED the "preliminary run" step')

        if t_step > self.sumd_param['t_max']: #Riassegnare velocità e geometria (partire dall'inizio)
          print(f"\nThe configuration under consideration has failed for {t_step} times. The system will be updated with new geometry and new speeds")
          t_step = 0
          corrette = 1 #bisogna ricominciare dall'inizio
          cambiare_velocità = 0
          BuildGeometry(self.PDBPROT,self.PDBLIG,self.TopFOLD,solvate=True,d=self.distance_from_edge,dist=self.distance_from_pocket,pocket=self.POCKET)
          minimization_equilibration(MF_m = mdp_em, MF_e = mdp_nvt)
          shutil.rmtree(act_sim), shutil.copytree(equil, act_sim) #Copio l'output dell'equalizzazione nella cartella actual_simulation
          shutil.copytree(corr_sim,os.path.join(fal_sim,f"suMD{fallimento}")) #tengo conto del fallimenti e salvo tutto in fal_sim
          shutil.rmtree(corr_sim), os.mkdir(corr_sim) #La cartella delle vecchie simulazioni corrette deve essere eliminata dato che si ricomincia dall'inizio
          fallimento +=1

      avviso_suMD = True
      flag2 = 1
      while flag2 == 1: #suMD
        
        if avviso_suMD:
          print("\nsuMD: The algorithm has ENTERED the 'suMD' step\n")
          avviso_suMD = False
      
        if cambiare_velocità == 0: #Eseguire MD senza cambiare le velocità
          print(f"MD simulation #{corrette} ...", end='')
          self.MdMdp.mdp_params['gen-vel']='no'
          self.MdMdp.mdp_params['continuation']='yes'
          self.MdMdp.write_mdp(mdp_fold+'04-md.mdp')
          m_primo,dcm_vector = SuMdStep(self.sumd_param['n_sample'],MF=self.MdMdp.mdpfile_out,append=True)
        else: #Eseguire MD cambiando le velocità
          print(f"MD simulation #{corrette} - change of velocity ...", end='')
          self.MdMdp.mdp_params['gen-vel']='yes'
          self.MdMdp.mdp_params['continuation']='no'
          self.MdMdp.write_mdp(mdp_fold+'04-md.mdp')
          print(f"MD simulation #{corrette} - change of velocity ...", end='')
          m_primo,dcm_vector = SuMdStep(n,MF=self.MdMdp.mdpfile_out,append=True)
        if m_primo >= 0: #Non produttivo
          cambiare_velocità = 1 #Cambiare le velocità dopo essere rientrato nel while   
          print(f" The step was not productive (Consecutive attempt {counter1})")
          counter1 += 1
          shutil.rmtree(und_sim), os.mkdir(und_sim) #elimino cartella dei risultati
        else:
          print(f" The step was productive (Distance of the last step: {dcm_vector[-1]:.2f} Anstrong)")
          
          shutil.copytree(und_sim,os.path.join(corr_sim,f"sim{corrette}")) #In modo da raccogliere tutte le simulazioni corrette
          shutil.rmtree(act_sim), shutil.copytree(und_sim,act_sim) #In modo che alla prossima simulazione md si riparta dall'ultimo file .gro (le ultime velocità)
          shutil.rmtree(und_sim), os.mkdir(und_sim) #elimino cartella dei risultati
          corrette += 1

          if dcm_vector[-1] > self.sumd_param['distanza_suMD_FS']: #distanza ligando-proteina maggiore di distanza_suMD_FS
            k_step = c_step + 1
            counter1 = 0
            c_step += 1
            cambiare_velocità = 0 # avviare l'md run senza cambiare le velocità
          else:
            flag2 = 0
            print(f"\nProtein and ligand are less than {self.sumd_param['distanza_suMD_FS']} Angstrom apart ({dcm_vector[-1]:.2f} Angstrom)")
            print('\nsuMD: The algorithm has EXITED the "suMD" step')

        if counter1 >= self.sumd_param['counter1t']: #Riassegnare velocità e geometria (partire dall'inizio)
          print(f"\nThe configuration under consideration has failed for {counter1} times. The system will be updated with new geometry and new speeds")
          corrette = 1 #bisogna ricominciare dall'inizio
          t_step = 1 # counter that is increased by 1 and compared with t_max every time there is a failed step in the preliminary run
          counter1 = 0 # counter that is incremented by 1 and compared with counter1t every time there is a failed step in the SuMD run
          c_step = 0 # counters that allow us to determine when a SuMD run step is productive
          k_step = 0
          cambiare_velocità = 1 #bisogna cambiare geometria e velocità

          BuildGeometry(self.PDBPROT,self.PDBLIG,self.TopFOLD,solvate=True,d=self.distance_from_edge,dist=self.distance_from_pocket,pocket=self.POCKET)
          minimization_equilibration( MF_m = self.EmMdp.mdpfile_out, MF_e = self.NvtMdp.mdpfile_out)
          shutil.rmtree(act_sim), shutil.copytree(equil, act_sim) #Copio l'output dell'equalizzazione nella cartella actual_simulation
          shutil.copytree(corr_sim,os.path.join(fal_sim,f"suMD{fallimento}")) #tengo conto del fallimenti e salvo tutto in fal_sim
          shutil.rmtree(corr_sim), os.mkdir(corr_sim) #La cartella delle vecchie simulazioni corrette deve essere eliminata dato che si ricomincia dall'inizio
          fallimento +=1
        
      avviso_final_step = True
      flag3 = 0
      while flag3 == 0: #FINAL STEP
        if avviso_final_step:
            print("\nsuMD: The algorithm has ENTERED the final step\n")
            avviso_final_step = False

        #MD finale quando la distanza ligando-proteina è minore di distanza_suMD_FS
        print(f"MD simulation #{corrette} ...", end='')
        self.MdMdp.mdp_params['gen-vel']='no'
        self.MdMdp.mdp_params['continuation']='yes'
        self.MdMdp.write_mdp(mdp_fold+'04-md.mdp')
        m_primo,dcm_vector = mdrun(self.sumd_param['n'],MF=mdp_classic,append=True)
        d_cm_out = dcm_vector[-1]

        if d_cm_out > self.sumd_param['distanza_suMD_FS']:
          print(f"\nThe ligand has broken away from the protein (distance > {self.sumd_param['distanza_suMD_FS']} Angstrom)")
          flag3 = 1 # in questo modo si esce dal final step per poi rientrare nel suMD step
          cambiare_velocità = 1 #nella prossima simulazione in suMD si deve fare in modo da cambiare le velocità degli atomi
          m_primo = -1 # Imposto arbitrariamente m_primo < 0 in modo che non entrerà nel preliminary run
          counter1 = 0 #inizializzo il contatore dei fallimenti del suMD step
          shutil.rmtree(und_sim), os.mkdir(und_sim) #elimino cartella dei risultati
        else:
          shutil.copytree(und_sim,os.path.join(corr_sim,f"sim{corrette}")) #In modo da raccogliere tutte le simulazioni corrette
          shutil.rmtree(act_sim), shutil.copytree(und_sim,act_sim) #In modo che alla prossima simulazione md si riparta dall'ultimo file .gro (le ultime velocità)
          shutil.rmtree(und_sim), os.mkdir(und_sim)
          corrette += 1
            
          if d_cm_out>=0 and d_cm_out<=2:   
            counter02 += 1
            if counter02 > self.sumd_param['ct02']:
              print(" Binding site reached")
              flag3 = 1 
              flag = 0
          
          if d_cm_out>=2 and d_cm_out<=5: 
            counter25 += 1  
            if counter25> self.sumd_param['ct25']:     
              print(" Neighbor binding site reached")
              flag3 = 1
              flag = 0

          if d_cm_out>=5 and d_cm_out<=self.sumd_param['distanza_suMD_FS']: 
            counter59 += 1   
            if counter59> self.sumd_param['ct59']:  
              print(" Meta binding site reached")  
              flag3 = 1
              flag = 0
          allert=f'''\nBindig site step {counter02}/{self.sumd_param['ct02']}
    Neighbour site step {counter25}/{self.sumd_param['ct25']}
    Meta binding site {counter59}/{self.sumd_param['ct59']}'''
          print(allert)


if __name__=="__main__":
  start_time=time.localtime()

  # Creare la scritta stilizzata
  text_sumd = pyfiglet.figlet_format("suMD", font="slant", width=80, justify="center")
  # Definire la larghezza della cornice
  border_width = 80

  # Creare la cornice intorno alla scritta
  border = '+' + '-' * border_width + '+\n'

  # Centrare la scritta all'interno della larghezza della cornice
  centered_text = text_sumd.center(border_width)

  # Stampare la scritta con la cornice centrata
  print(border + centered_text + '\n' + border)
  print(time.asctime(start_time))
