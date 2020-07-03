# REFINE MEMBRANE PROTEIN STRUCTURE WITH OS-SSNMR CSA AND DC RESTRAINTS
# Written by D. K. Weber (last revised May 14 2020)
#
# DESCRIPTION
# -----------
# The XPLOR-NIH script implements the hybrid method (Shi...Veglia, 2009, DOI: 10.1007/s10858-009-9328-9)
# with several updates to utilize new functions and force fields in XPLOR-NIH. The script requires a
# well folded structure(s) to start. The topology of the structure is refined by geometrically 
#
# KEY UPDATES
# -----------
# - Choice of utilizing EEFx/IMMx potentials introduced by Marassi Lab
# - Choice of using Ez potentials for depth of insertion (not implemented yet)
# - Use of torsion and distance restraints optional (terms not used if restraints table not specified)
#
# PROTOCOL
# --------
# This follows the refinement protocol from Tian...Marassi, 2015, DOI: 10.1016/j.bpj.2015.06.047
# The template in the XPLOR-NIH tutorials (eefx-membrane) was used. Modified
# to use EEFx or E-Repel, appropriate DC/CSA weighting factors, optional restraints.
#
# 0. Get input PDB(s) from previous simulated annealing procedure
#    Use structure pre-inserted into membrane if using EEFx. The protocol wont forcibly translate along Z
# 1. Initial Powell Cartesian minimization (500 steps)
# 2. High temperature torsion dynamics (3000 K for 10 ps, 10000 steps)
# 3. Torsion angle dynamics with simulated annealing (3000 K to 25 K in 12.5 K steps, 0.2 ps/200 steps per increment)
# 4. Powell torsion angle minimization (500 steps)
# 5. Powell Cartesian minimization (500 steps)

xplor.requireVersion('3.0')
import protocol
import argparse
import numpy as np

# Argument Parser
def parse_args():
    parser = argparse.ArgumentParser(description='Initial folding for membrane protein structure.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--structure_in', type=str,
        help='Starting structure (PDB file)'
    )
    parser.add_argument(
        '--DC_NH', type=str,
        help='DC restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1', type=str,
        help='CSA restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_gly', type=str,
        help='CSA restraints table for glycines. Default: None', default=''
    )
    parser.add_argument(
        '--DIHE', type=str,
        help='Dihedral restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--DC_NH_max', type=float,
        help='Maximum dipolar coupling. Default: 10.735', default=10.735
    )
    parser.add_argument(
        '--nstructures', type=int,
        help='Dihedral restraints table. Default: 100', default=100
    )
    parser.add_argument(
        '--CSA_N1_tensor', type=float, nargs='+',
        help='Non-glycine N15 principal axis system components in ppm (57.3, 81.2, 228.1).',
        default=(57.3,81.2,228.1)
    )
    parser.add_argument(
        '--CSA_N1_tensor_gly', type=float, nargs='+',
        help='Glycine N15 principal axis system components in ppm (45.6, 66.3, 211.6).',
        default=(45.6,66.3,211.6)
    )
    parser.add_argument(
        '--CSA_N1_beta', type=float,
        help='Non-glycine beta Euler angle for N15 tensor. Default: -17.0', default=-17.0
    )
    parser.add_argument(
        '--CSA_N1_beta_gly', type=float,
        help='Glycine beta Euler angle for N15 tensor. Default: -21.6', default=-21.6
    )
    parser.add_argument(
        '--tm_domain', type=int, nargs='+',
        help='Start/Stop Residue IDs of helical segments. Default: None',
        default=''
    )
    parser.add_argument(
        '--immx_thickness', type=float,
        help='IMMx Membrane thickness. Default: 25.8 (DMPC/POPC bicelle)', default=25.8
    )
    parser.add_argument(
        '--immx_nparameter', type=int,
        help='IMMx n parameter. Default: 10', default=10
    )
    parser.add_argument(
        '--w_slf', type=float,
        help='SLF/CDIH force scaling. Default: 5.0', default=5.0
    )
    parser.add_argument(
        '--w_r', type=float,
        help='DC/CSA force scaling. Default: 3.0', default=3.0
    )
    parser.add_argument(
        '--HBDA', type=str,
        help='H-bond restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--NOE', type=str,
        help='NOE distance restraints table. Default: None', default=''
    )
    args = parser.parse_args()
    return args
    

# Read arguments
args = parse_args()
tm_domain = 'resid %s:%s and (name CA)' % (args.tm_domain[0],args.tm_domain[1])
nstructures = args.nstructures         # Number of structures to calculate
outname = "SCRIPT_STRUCTURE.sa"        # File prefix same as python script
seed = 3421                            # Random seed


# Input files
structure_in = args.structure_in
DC_NH  = args.DC_NH
CSA_N1 = args.CSA_N1
CSA_N1_gly = args.CSA_N1_gly
DIHE = args.DIHE
NOE = args.NOE
HBDA = args.HBDA


# 1H-15N DC
DC_NH_max = args.DC_NH_max  # Maximum dipolar coupling in kHz


# 15N tensor (correct to reduced sigma notation and reorder for XPLOR-NIH input)
sxx, syy, szz = args.CSA_N1_tensor; iso = np.mean(args.CSA_N1_tensor)
CSA_N1_tensor = (iso-sxx, iso-szz, iso-syy)
sxx, syy, szz = args.CSA_N1_tensor_gly; iso = np.mean(args.CSA_N1_tensor_gly)
CSA_N1_tensor_gly = (iso-sxx, iso-szz, iso-syy)
CSA_N1_beta = args.CSA_N1_beta
CSA_N1_beta_gly = args.CSA_N1_beta_gly


# EEFx Membrane parameters
# Membrane hydrophobic thickness: DMPC (25.4A), DPPC (28.6A), POPC (27.0A), DOPC (29.6A)
# See Marsh, Handbook of lipid bilayers, 2nd Ed. p. 379. (FROM XPLOR EXAMPLES)
immx_com = 'resid %s:%s and (name CA)' % (args.tm_domain[0],args.tm_domain[1])
immx_thickness = args.immx_thickness        # 25.4*0.75 + 27.0*0.25 = 25.8 A for DMPC/POPC bicelle ()
immx_nparameter = args.immx_nparameter      # IMMx n parameter of membrane profile


# Set force constants in kcal/mol/(err^2) - initial and final values
ini_tordb = 0.002   ;	fin_tordb = 2.0   # Statistical torsion angles - Marassi setting
ini_bond  = 1.0	    ;	fin_bond  = 1.0   # Bond lengths - Marassi settings
ini_angl  = 0.4	    ;	fin_angl  = 1.0   # Bond angles - Marassi settings
ini_impr  = 0.1	    ;	fin_impr  = 1.0   # Improper torsion angles - Marassi settings
ini_dihd  = 10.0    ;	fin_dihd  = 200.0 # Experimental torsions (CDIH)
ini_noe   = 2.0	    ;	fin_noe	  = 40.0  # NOE distances - Veglia setting
ini_vrep  = 0.90    ;	fin_vrep  = 0.75  # VDW Repulsion - Veglia seeting (ignored with EEFX)
ini_vrcn  = 0.001   ;	fin_vrcn  = 10.0  # 


# Set intial and final CSA and DC force constants - based on weighting factors to CDIH terms
# See for details: Shi et al (2009) J Biol NMR, DOI: 10.1007/s10858-009-9328-9 
w_slf = args.w_slf    # SLF/CDIH force scaling 
w_r = args.w_r        # DC/CSA force scaling
ini_dc    = 0.00010 * w_slf *  w_r/(w_r+1)   # Dipolar coupling force initial
fin_dc    = fin_dihd * w_slf * w_r/(w_r+1)   # Dipolar coupling force final
ini_cs    = 0.00010 * w_slf * 1/(w_r+1)      # CSA force inital
fin_cs    = fin_dihd * w_slf * 1/(w_r+1)     # CSA force final


# Initial structure
ini_model=structure_in
import psfGen
psfGen.pdbToPSF(ini_model)
protocol.genExtendedStructure()


# Load EEFX2 topology
protocol.parameters['protein']="eefx/protein_eefx2.par"
protocol.topology['protein']  ="eefx/protein_eefx2.top"


# Initialize potentials and parameters for high temperature dynamics simulated annealing
pots = PotList()
from simulationTools import StaticRamp, MultRamp
from simulationTools import InitialParams
highTempParams1 = []                  # Settings for high T torsion dynamics - Stage 1 EEFX
highTempParams2 = []                  # Settings for high T torsion dynamics - Stage 2 EEFX
highTempParams = []                   # Settings for high T torsion dynamics - REPEL and Stage 3 EEFX
rampedParams = []                     # Settings for annealing stage
lowTempParams = []                    # Settingd for low T dynamics (i.e., Powell minimization)


### 1H-15N DIPOLAR COUPLING RESTRAINTS ###

# Initialize DC tensor
dc_tensor = {}
from varTensorTools import create_VarTensor
oTensor = create_VarTensor('amide_NH')
oTensor.setDaMax(DC_NH_max/2)			      			
oTensor.setDa(DC_NH_max/2)
oTensor.setRh(0)
dc_tensor['amide_NH'] = oTensor

# Load working DC restraints if specified
if DC_NH:
    dc = PotList('DIPL')
    from rdcPotTools import create_RDCPot
    for (name, file, tensor, scale) in [('amide_NH', DC_NH, dc_tensor['amide_NH'], 1)]:
        DIPL=create_RDCPot(name=name, file=file, oTensor=tensor)
        DIPL.setScale(scale)
        DIPL.setShowAllRestraints(True)
        DIPL.setUseSign(False)
        DIPL.setThreshold(0.01)
        DIPL.setPotType('square')
        dc.append(DIPL)
        pass
    pots.append(dc)
    #highTempParams.append(StaticRamp("dc.setScale(ini_dc)"))
    rampedParams.append(MultRamp(ini_dc, fin_dc, "dc.setScale(VALUE)"))
    lowTempParams.append(StaticRamp("dc.setScale(fin_dc)"))
    
    for t in dc_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass


### CSA RESTRAINTS ###

# Initialize CSA tensor
csa_tensor = {}
from varTensorTools import create_VarTensor
oTensor = create_VarTensor('amide_N')
oTensor.setDaMax(DC_NH_max)			      			
oTensor.setDa(DC_NH_max/2)
oTensor.setRh(0)
csa_tensor['amide_N'] = oTensor

# Load working CSA restraints if specified

### DEVELOPMENT SECTION
# Add more CSA tensors here if need. Ensure AtomOrder is corrent in tbl file.
CSA_all = []
if CSA_N1:
    CSA_all.append(('amide_N1', CSA_N1, CSA_N1_tensor, CSA_N1_beta, 1))
if CSA_N1_gly:
    CSA_all.append(('amide_N1_gly', CSA_N1_gly, CSA_N1_tensor_gly, CSA_N1_beta_gly, 1))
### END DEVELOPMENT SECTION

if CSA_all:
    csa = PotList('CS')
    from csaPotTools import create_CSAPot
    for (name, file, atom_sigma, ang_b, scale) in CSA_all:
        CS=create_CSAPot(name=name, file=file, oTensor=csa_tensor['amide_N'])
        CS.setScale(scale)
        CS.setTensorClass("bond")
        CS.setThreshold(0.01)
        CS.setPotType('square')
        CS.setSigma(atom_sigma)
        CS.setBeta(ang_b)
        CS.setAtomOrder("123")
        CS.setShowAllRestraints(True)
        CS.setVerbose(True)
        CS.setDaScale(DC_NH_max*-1)
        csa.append(CS)
    pots.append(csa)
    #highTempParams.append(StaticRamp("csa.setScale(ini_cs)"))
    rampedParams.append(MultRamp(ini_cs, fin_cs, "csa.setScale(VALUE)"))
    lowTempParams.append(StaticRamp("csa.setScale(fin_cs)"))

    for t in csa_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass


### DIHEDRAL RESTRAINTS ###
    
# Load dihedral restraints (i.e., TALOS) if specified
if DIHE:
    protocol.initDihedrals(DIHE)
    from xplorPot import XplorPot
    pots.append(XplorPot('CDIH'))
    highTempParams1.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    highTempParams2.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    highTempParams.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    rampedParams.append(StaticRamp("pots['CDIH'].setScale(fin_dihd)"))
    lowTempParams.append(StaticRamp("pots['CDIH'].setScale(fin_dihd)"))
    pots['CDIH'].setThreshold(5.00)     # dflt [2.0]

    
### CHEMISTRY KNOWLEDGE ###
    
# Initialize chemisty knowledge
from torsionDBPotTools import create_TorsionDBPot
torsionDB = create_TorsionDBPot(name='torsionDB')
pots.append(torsionDB)
rampedParams.append(MultRamp(ini_tordb, fin_tordb, "torsionDB.setScale(VALUE)"))
lowTempParams.append(StaticRamp("pots['torsionDB'].setScale(fin_tordb)"))

from xplorPot import XplorPot
pots.append(XplorPot('BOND'))          # Dflt scale [1]
pots.append(XplorPot('ANGL'))          # Dflt scale [1]
pots.append(XplorPot('IMPR'))          # Dflt scale [1]
rampedParams.append(MultRamp(ini_bond, fin_bond, "pots['BOND'].setScale(VALUE)"))
rampedParams.append(MultRamp(ini_angl, fin_angl, "pots['ANGL'].setScale(VALUE)"))
rampedParams.append(MultRamp(ini_impr, fin_impr, "pots['IMPR'].setScale(VALUE)"))
lowTempParams.append(StaticRamp("pots['BOND'].setScale(fin_bond)"))
lowTempParams.append(StaticRamp("pots['ANGL'].setScale(fin_angl)"))
lowTempParams.append(StaticRamp("pots['IMPR'].setScale(fin_impr)"))

# Set threshold for violations of PotList terms
pots['BOND'].setThreshold(0.05)     # dflt [0.05]
pots['ANGL'].setThreshold(5.00)     # dflt [2.0]
pots['IMPR'].setThreshold(5.00)     # dflt [2.0]


### H-BOND RESTRAINTS ###

if HBDA:
    # hbda - distance/angle bb hbond term
    protocol.initHBDA(HBDA)
    pots.append( XplorPot('HBDA') )

    # hbdb - knowledge-based backbone hydrogen bond term
    protocol.initHBDB()
    pots.append( XplorPot('HBDB') )


### NOE DISTANCE RESTRAINTS

if NOE:
    dsts = PotList('NOE')
    from noePotTools import create_NOEPot
    for (exp, file, scale) in [('noe', NOE, 1)]:
        NOE = create_NOEPot(exp,file,)
        NOE.setScale(scale)
        NOE.setPotType('soft')
        dsts.append(NOE)
    pots.append(dsts)
    rampedParams.append(MultRamp(ini_noe, fin_noe, "dsts.setScale(VALUE)"))
    lowTempParams.append( StaticRamp("dsts.setScale(fin_noe )"))


# Initialize IMMx membrane
Zpos = 0
from eefxPotTools import create_EEFxPot, param_LK
from eefxPotTools import setCenter, setCenterXY
eefx=create_EEFxPot("EEFX","ALL",paramSet=param_LK,verbose=False)
eefx.setScale(1)
#eefx.setVerbose(1)
eefx.setIMMx(1)
eefx.useGROUp(1)
eefx.setMoveTol(0.5)
print(eefx.showParam())
eefx.setThickness(immx_thickness)
eefx.setProfileN(immx_nparameter)
eefx.setA(0.85)
pots.append(eefx)
setCenter(immx_com, Zpos)		# Translate selected center of mass to IMMx Zpos.
highTempParams1.append(StaticRamp("eefx.setScale(0)"))
highTempParams.append(StaticRamp("eefx.setScale(0.004)"))
rampedParams.append(MultRamp(0.004,1.0,"eefx.setScale(VALUE)"))
lowTempParams.append(StaticRamp("eefx.setScale(1.0)"))
    
    
# Initialize REPEL if EEFX not used - XPLOR-NIH Tutorial Settings - eefx-membrane
pots.append(XplorPot('VDW'))          # dflt scale [1]
highTempParams1.append(StaticRamp("pots['VDW'].setScale(0.004)"))
highTempParams1.append(StaticRamp("""protocol.initNBond(cutnb=100,
                                                         repel=1.2,
                                                         nbxmod=4,
                                                         tolerance=45,
                                                         onlyCA=1)"""))
# standard all-atom repel
highTempParams2.append(StaticRamp("protocol.initNBond(nbxmod=4)"))
highTempParams2.append(StaticRamp("pots['VDW'].setScale(0.1)"))
highTempParams.append(StaticRamp("pots['VDW'].setScale(0)"))
rampedParams.append(StaticRamp("pots['VDW'].setScale(0)"))
lowTempParams.append(StaticRamp("pots['VDW'].setScale(0)"))
    

# SET UP IVM OBJECTS (dyn, minc) THAT PERFORM DYNAMICS AND MINIMIZATION.
# IVM (internal variable module) is used to perform dynamics and minimization in
# both torsion-angle and Cartesian space. Bonds, angles and many impropers cannot
# change with internal torsion-angle dynamics.
#===============================================================================
from ivm import IVM
dyn = IVM()                     # IVM object for torsion-angle dynamics.
dyn.reset()                     # reset ivm topology for torsion-angle dynamics.
protocol.torsionTopology(dyn)

minc = IVM()                    # IVM object for Cartesian minimization.
protocol.cartesianTopology(minc)


# DYNAMICS SETTINGS FOR HIGH T AND ANNEALING STAGES.
#===============================================================================
ini_temp = 3500.0   ; fin_temp = 25.0   # Initial and final temperatures.
protocol.massSetup()                    # Give atoms uniform weights except for axes.


# Calculate structure module for folding - Standard protocol
def calcOneStructure(loopInfo):
    """
    This function calculates a single structure, performs analysis on the
    structure and then writes out a pdb file with remarks.
    """        
    # Generate initial structure and minimize.
    #===========================================================================
    # Generate initial structure by randomizing torsion angles.
    # Then set torsion angles from restraints (this shortens high T dynamics).
    # Then set selected center of mass to membrane center (IMMx z=0).
    import monteCarlo
    monteCarlo.randomizeTorsions(dyn)
    protocol.fixupCovalentGeom(maxIters=100, useVDW=1)
    import torsionTools
    if DIHE:
        import torsionTools
        torsionTools.setTorsionsFromTable(DIHE)
    
    InitialParams(rampedParams)     # parameters for SA.
    InitialParams(highTempParams1)  # reset some rampedParams.
    setCenter(immx_com, Zpos)        # Translate selected center of mass to IMMx Zpos.

    # Torsion angle minimization.
    protocol.initMinimize(dyn,
                          potList=pots,
                          numSteps=100,
                          printInterval=50)
    dyn.run()

    
    # High temperature dynamics.
    #===========================================================================
    # High temperature dynamics stage 1.
    protocol.initDynamics(dyn,
                          potList=pots,         # potential terms to use.
                          bathTemp=ini_temp,    # set bath temperature.
                          initVelocities=1,     # uniform initial velocities.
                          finalTime=3,          # run for finalTime or
                          numSteps=3001,        # numSteps * 0.001, whichever is less.
                          printInterval=100)    # printing rate in steps.
    dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
    dyn.run()
    
    # High temperature dynamics stage 2.
    InitialParams(highTempParams2)    
    setCenter(immx_com, Zpos)            # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,         # potential terms to use.
                          bathTemp=ini_temp,    # set bath temperature.
                          initVelocities=1,     # uniform initial velocities.
                          finalTime=3,          # run for finalTime or
                          numSteps=3001,        # numSteps * 0.001, whichever is less.
                          printInterval=100)    # printing rate in steps.
    dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
    dyn.run()

    # High temperature dynamics stage 3.
    InitialParams(highTempParams)
    setCenter(immx_com, Zpos)            # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,         # potential terms to use.
                          bathTemp=ini_temp,    # set bath temperature.
                          initVelocities=1,     # uniform initial velocities.
                          finalTime=26,			# run for finalTime or
                          numSteps=26001,		# numSteps * 0.001, whichever is less.
                          printInterval=100)    # printing rate in steps.
    dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
    dyn.run()

    
    # Initialize integrator and loop for simulated annealing and run.
    #===========================================================================
    # Dynamics for annealing.
    setCenter(immx_com, Zpos)            # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,
                          finalTime=0.4,        # run for finalTime or
                          numSteps=201,         # numSteps*0.001, whichever is less.
                          printInterval=100)    

    # Set up cooling loop and run.
    from simulationTools import AnnealIVM
    AnnealIVM(initTemp=ini_temp,
                       finalTemp=fin_temp,
                       tempStep=12.5,
                       ivm=dyn,
                       rampedParams=rampedParams
                       ).run()

    
    # Final minimization.
    #===========================================================================
    # Torsion angle minimization.
    protocol.initMinimize(dyn,
                          numSteps=500,         # dflt [500 steps]
                          potList=pots,
                          printInterval=50)
    dyn.run()

    # Final Cartesian minimization.
    protocol.initMinimize(minc,
                          numSteps=500,         # dflt [500 steps]
                          potList=pots,
                          dEPred=10)
    minc.run()

    # Recenter coordinates in XY plane.
    setCenterXY()                       # translate protein coordinates to XY center.
    
    # Do analysis and write structure when this routine is finished.
    pass
    

# Loop control for folding stage
from simulationTools import StructureLoop, FinalParams
StructureLoop(structLoopAction=calcOneStructure,
              numStructures=nstructures,
              pdbTemplate=outname,
              doWriteStructures=1,              # analyze and write coords after calc
              genViolationStats=1,              # print stats file
              averageContext=FinalParams(lowTempParams),
              averageFitSel="",                 # selection for bkbn rmsd [CA]
              averageCompSel="",                # selection for heavy atom rmsd
              averageTopFraction=0.1,           # Report stats on top 10%
              averagePotList=pots,              # Terms for stats and avg
              averageSortPots=[                 # Terms used to sort models
                  term for term in pots \
                  if term.instanceName() not in ("torsionDB","VDW")],
).run()


