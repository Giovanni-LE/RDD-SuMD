# # IMPORT DELLE LIBRERIE
# %%
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


# Creare la scritta stilizzata
text = pyfiglet.figlet_format("suMD", font="slant", width=80, justify="center")
# Definire la larghezza della cornice
border_width = 80

# Creare la cornice intorno alla scritta
border = '+' + '-' * border_width + '+\n'

# Centrare la scritta all'interno della larghezza della cornice
centered_text = text.center(border_width)

# Stampare la scritta con la cornice centrata
print(border + centered_text + '\n' + border)



sys.exit(0)
def Linear_Fitting(PDB,XTC,n,POCKET=False,LIG="resname LIG"): #Funzione per il fitting lineare
    '''
    This function calculte the angular coefficent of the distance between POCKET and LIGAND over time
    PDB : str
    the path of the pdb file
    XTC :str 
    the path of the trajectory file
    n : int
    ?????
    POCKET : str
    the selection of the binding pocket
    LIG: str
    the selection of the ligand
    '''  
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

    """ print('\nVorresti visualizzare anche la curva non fittata? Se si premi 1, altrimenti premi un tasto qualsiasi')
      ch = input()"""

    """#if ch == str(1):
    plt.plot(x,d.astype(np.float32), label = "RAW Distance") 
    plt.xlabel('Time(ps)')
    plt.ylabel('Distance (Angstrom)')
    plt.show

    plt.plot(x, m1[0]*x + m1[1], marker='o', label = "Linear Fitting") #plotting
    plt.xlabel('Time(ps)')
    plt.ylabel('Distance (Angstrom)')
    plt.show 
    plt.legend()"""
    return(m1[0],d)

def g_mdrun(OutputFile, xtcFile, log_File, edr_File,groFile, TprFile=None,gpu=False,append=False):

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
    

# %% [markdown]
# Funzione mdrun senza cambiare velocità degli atomi

# %%

#TODO tirare fuori la selezione del ligando dalla funzione

def mdrun(n,MF,append):
    
    Gromacs_Tools.g_grompp(MdpFile=MF, TopolFile= os.path.join(outdir,"topol","topol.top"), CoordFile= os.path.join(act_sim,"md.gro"), TprFile= os.path.join(und_sim,"md.tpr"), IndexFile= os.path.join(minim,"index.ndx"),maxwarning='2')   
    if not os.path.exists(os.path.join(und_sim,"md.tpr")):
        sys.exit('error during the creation of tpr file (grompp)')
        
    g_mdrun(OutputFile= os.path.join(und_sim,"md.trr"),xtcFile= os.path.join(und_sim,"md.xtc"), log_File= os.path.join(und_sim,"md.log"),edr_File=os.path.join(und_sim,"md.edr"),groFile=os.path.join(und_sim,"md.gro"), TprFile=os.path.join(und_sim,"md.tpr"),gpu=gpu,append=append)
    if not os.path.exists(os.path.join(und_sim,"md.xtc")):
        sys.exit('error during the creation of xtc file (mdrun)')
    
    Gromacs_Tools.g_editconf(f=os.path.join(und_sim,"md.tpr"), o=os.path.join(und_sim,"md.pdb")) # per avere il file pdb
    if not os.path.exists(os.path.join(und_sim,"md.pdb")):
        sys.exit('error during the creation of pdb file (editconf)')

    m_primo,dcm_vector = Linear_Fitting(PDB=os.path.join(und_sim,"md.tpr"),XTC = os.path.join(und_sim,"md.xtc"),n=n,POCKET='resid 267 or resid 252 or resid 264 or resid 256 or resid 174')
    return m_primo,dcm_vector


def BuildGeometry(PDBPROT=None,PDBLIG=None,TopFOLD=None,newbox=True,d='5',mindist=100,solvate=False,Saltconc="0.15"):
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
    cent: the '''
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

# %% [markdown]
# Funzione minimization per minimizzare il sistema

# %%
def minimization_equilibration(MF_m,MF_e):

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



#----------------------------------------------------
# Richiedi all'utente se una GPU è disponibile
while True:
    print('\nIs a GPU available? (Y|N)')
    gpu = input().strip().upper()
    if gpu == 'Y':
        gpu = True
        break
    elif gpu == 'N':
        gpu = False
        break
    else:
        print('Invalid input. Please enter Y or N.')

#----------------------------------------------------
# Richiedi all'utente il pdb della proteina e del ligando e la cartella da cui trarre la topologia
print("\nSelect the path to the PROTEIN pdb file")
#PDBPROT = askopenfilename()
PDBPROT="PDB/3RFM/0_sys/3RFM.gro"

# stampa il percorso del file selezionato
print("path to the protein pdb file:", PDBPROT)

print("\nSelect protein residues (ex. resid ... or resid ...)")
#POCKET = input()
POCKET='resid 267 or resid 252 or resid 264 or resid 256 or resid 174'
print("protein residues are:", POCKET)

print("\nSelect the path to the LIG pdb file")
#PDBLIG = askopenfilename()
PDBLIG="PDB/LIG/LIG.pdb"

# stampa il percorso del file selezionato
print("path to the lig pdb file:", PDBLIG)

print("\nSelect the path to the topology folder")
#PDBLIG = askopenfilename()
TopFOLD="PDB/3RFM_LIG_NOWATER/topol/"

# stampa il percorso del file selezionato
print("path to the topology folder:", TopFOLD)

#----------------------------------------------------
#Richiedi all'utente il numero di step
while True:
    try:
        n = int(input("Insert an integer number between 3 and 61: "))
        if n < 3 or n > 61:
            print("The number is not available. Please enter a number between 3 and 61.")
        else:
            break
    except ValueError:
        print("Invalid input. Please enter an integer number between 3 and 61.")

print("\nThe input number is:", n)  

# %%
m = 0 # threshold of the angular coefficient values of the straight line that interpolates the saved points
t_max = 31 # maximum threshold of bankruptcy attempts granted in the preliminary run

counter1t = 17 # maximum threshold of failed attempts granted in the SuMD run
ct02 = 19 # threshold relating to the times in which the distance between the com is between 0 and 2 Å during the final phase of unsupervised MD
ct25 = 19 # threshold relating to the times in which the distance between the com is between 2 and 5 Å during the final phase of unsupervised MD
ct59 = 19 # threshold relating to the times in which the distance between the com is between 5 and 9 Å during the final phase of unsupervised MD

#Inizializzazione
m_primo = 0 # value of the angular coefficient of the line which is updated and compared with m each short MD simulation
t_step = 1 # counter that is increased by 1 and compared with t_max every time there is a failed step in the preliminary run
counter1 = 0 # counter that is incremented by 1 and compared with counter1t every time there is a failed step in the SuMD run
c_step = 0 # counters that allow us to determine when a SuMD run step is productive
k_step = 0 # counters that allow us to determine when a SuMD run step is productive
d_cm_vector = np.zeros(5) # vector in which the distances of the com taken in the 5 steps are entered within each short MD simulation during the SuMD run
d_cm_out = 0 # final distance value dcm stored at the end of the short MD simulations that occur once supervision is stopped
counter02 = 0 # counter that is updated and compared with ct02 every time that  dcm is between 0 and 2 Å
counter25 = 0 # counter that is updated and compared with ct25 every time that dcm is between 2 and 5 Å
counter59 = 0 # counter that is updated and compared with ct59 every time that dcm is between 5 and 9 Å

# %%
#shutil.rmtree("sysprep/")
#shutil.rmtree("PDB/3RFM_LIG_NOWATER/05-md/")

workfold="sumd"#define workfold folder --> folder che conterrà ogni operazione dell'algoritmo
if not os.path.exists(workfold): #chekc if the workfold exists
  os.mkdir(workfold) #in case doesn't exist create the directory

fal_sim=os.path.join(workfold,"failure_simulation")# define the failure_simulation folder --> folder che conterrà i file fallimentari
if not os.path.exists(fal_sim): #chekc if the fal_sim exists
  os.mkdir(fal_sim) #in case doesn't exist create the directory

act_sim=os.path.join(workfold,"actual_simulation")# define the actual_simulation folder --> folder che conterrà i file in input alla simulazione md i-esima
if not os.path.exists(act_sim): #chekc if the act_sim exists
  os.mkdir(act_sim) #in case doesn't exist create the directory

und_sim=os.path.join(workfold,"undefined_simulation") #define und_sim folder --> folder che conterrà i file in output alla simulazione md i-esima
if not os.path.exists(und_sim): #chekc if the und_sim exists
  os.mkdir(und_sim) #in case doesn't exist create the directory

corr_sim=os.path.join(workfold,"correct_simulation") #define und_sim folder --> folder che conterrà i file in output alla simulazione md i-esima che soddisfano la condizione m>0
if not os.path.exists(corr_sim): #chekc if the corr_sim exists
  os.mkdir(corr_sim) #in case doesn't exist create the directory

outdir= os.path.join(workfold,"sysprep") #define output folder --> folder che conterrà i file in output alla minimizzazione ed equilibrazione dopo il cambio della geometria
if not os.path.exists(outdir): #chekc if the outputfolder exists
  os.mkdir(outdir) #in case doesn't exist create the directory

minim= os.path.join(outdir,"minimization") #define minimization folder --> folder che conterrà l'output della minimizzazione 
if not os.path.exists(minim): #chekc folder exists
  os.mkdir(minim) #in case doesn't exist create the directory

equil= os.path.join(outdir,"equilibration") #define equilibration folder --> folder che conterrà l'output dell'equalizzazione NVT 
if not os.path.exists(equil): #chekc folder exists
  os.mkdir(equil) #in case doesn't exist create the directory

'''mdp_em = 'mdp_nowater/02-em.mdp'
#mdp_nvt = 'mdp_nowater/03-nvt.mdp'
#mdp_velocity = 'mdp_nowater/06-md.mdp'
#mdp_classic = 'mdp_nowater/07-md.mdp'
mdp_nvt = 'mdp_nowater_acorciata/03-nvt.mdp'
mdp_velocity = 'mdp_nowater_acorciata/06-md.mdp'
mdp_classic = 'mdp_nowater_acorciata/07-md.mdp'

'''
mdp_em = 'mdp/02-em.mdp'
#mdp_nvt = 'mdp_nowater/03-nvt.mdp'
#mdp_velocity = 'mdp_nowater/06-md.mdp'
#mdp_classic = 'mdp_nowater/07-md.mdp'
mdp_nvt = 'mdp/03-nvt.mdp'
mdp_velocity = 'mdp/06-md.mdp'
mdp_classic = 'mdp/07-md.mdp'
#------------------------------------------------------------------------
prima_volta = 1
flag = 1
corrette = 1
fallimento = 1 #tiene conto delle volte in cui è stato necessario ricominciare la simulazione
distanza_suMD_FS = 9 #Angstrom - rappresenta la distanza alla quale vogliamo che si entri nel final step


def run(self):
  while flag == 1:
    if prima_volta == 1: #Se è la prima volta che si esegue lo script, allora è necessario generare un ambiente 
      print("Creation of the box containing protein and ligand\n")
      BuildGeometry(PDBPROT,PDBLIG,TopFOLD,solvate=True,d='0.5',dist=40,pocket=POCKET)
      minimization_equilibration(MF_m = mdp_em, MF_e = mdp_nvt)
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
        m_primo,dcm_vector = mdrun(n,MF=mdp_classic,append=True)
      else: #Eseguire MD cambiando le velocità
        print(f"MD simulation #{corrette} - change of velocity ...", end='')
        m_primo,dcm_vector = mdrun(n,MF=mdp_velocity,append=True)

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

      if t_step > t_max: #Riassegnare velocità e geometria (partire dall'inizio)
        print(f"\nThe configuration under consideration has failed for {t_step} times. The system will be updated with new geometry and new speeds")
        t_step = 0
        corrette = 1 #bisogna ricominciare dall'inizio
        cambiare_velocità = 0
        BuildGeometry(PDBPROT,PDBLIG,TopFOLD,solvate=True,d='0.5',dist=40,pocket=POCKET)
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
        m_primo,dcm_vector = mdrun(n,MF=mdp_classic,append=True)
      else: #Eseguire MD cambiando le velocità
        print(f"MD simulation #{corrette} - change of velocity ...", end='')
        m_primo,dcm_vector = mdrun(n,MF=mdp_velocity,append=True)

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

        if dcm_vector[-1] > distanza_suMD_FS: #distanza ligando-proteina maggiore di distanza_suMD_FS
          k_step = c_step + 1
          counter1 = 0
          c_step += 1
          cambiare_velocità = 0 # avviare l'md run senza cambiare le velocità
        else:
          flag2 = 0
          print(f'\nProtein and ligand are less than {distanza_suMD_FS} Angstrom apart ({dcm_vector[-1]:.2f} Angstrom)')
          print('\nsuMD: The algorithm has EXITED the "suMD" step')

      if counter1 >= counter1t: #Riassegnare velocità e geometria (partire dall'inizio)
        print(f"\nThe configuration under consideration has failed for {counter1} times. The system will be updated with new geometry and new speeds")
        corrette = 1 #bisogna ricominciare dall'inizio
        t_step = 1 # counter that is increased by 1 and compared with t_max every time there is a failed step in the preliminary run
        counter1 = 0 # counter that is incremented by 1 and compared with counter1t every time there is a failed step in the SuMD run
        c_step = 0 # counters that allow us to determine when a SuMD run step is productive
        k_step = 0
        cambiare_velocità = 1 #bisogna cambiare geometria e velocità

        BuildGeometry(PDBPROT,PDBLIG,TopFOLD,solvate=True,d='0.5',dist=40,pocket=POCKET)
        minimization_equilibration(MF_m = mdp_em, MF_e = mdp_nvt)
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
      m_primo,dcm_vector = mdrun(n,MF=mdp_classic,append=True)
      d_cm_out = dcm_vector[-1]

      if d_cm_out > distanza_suMD_FS:
        print(f"\nThe ligand has broken away from the protein (distance > {distanza_suMD_FS} Angstrom)")
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
          
          if counter02 > ct02:
            print(" Binding site reached")
            flag3 = 1 
            flag = 0
        
        if d_cm_out>=2 and d_cm_out<=5: 
          counter25 += 1  
          if counter25> ct25:     
            print(" Neighbor binding site reached")
            flag3 = 1
            flag = 0

        if d_cm_out>=5 and d_cm_out<=distanza_suMD_FS: 
          counter59 += 1   
          if counter59> ct59:  
            print(" Meta binding site reached")  
            flag3 = 1
            flag = 0
        allert=f'''\nBindig site step {counter02}/{ct02}
  Neighbour site step {counter25}/{ct25}
  Meta binding site {counter59}/{ct59}'''
        print(allert)


