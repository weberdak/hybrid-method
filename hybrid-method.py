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

# Argument Parser
def parse_args():
    parser = argparse.ArgumentParser(description='Initial folding for membrane protein structure.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--structures_in', type=str, nargs='+',
        help='Starting structure (PDB file)'
    )
    parser.add_argument(
        '--DC_NH_work_in', type=str,
        help='Working DC restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_work_in', type=str,
        help='Working CSA restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--DIHE_in', type=str,
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
        help='Non-glycine N15 principal axis system components in ppm (64.9,-105.9,41.0).',
        default=(64.9,-105.9,41.0)
    )
    parser.add_argument(
        '--CSA_N1_beta', type=float,
        help='Non-glycine beta Euler angle for N15 tensor. Default: -17.0', default=-17.0
    )
    parser.add_argument(
        '--eefx', type=bool,
        help='Use EEFx force field (True or False). Default: True', default=True, choices=(True, False)
    )
    parser.add_argument(
        '--helices', type=int, nargs='+',
        help='Residue IDs of helical segments. Default: None',
        default=''
    )
    parser.add_argument(
        '--immx_thickness', type=float,
        help='IMMx Membrane thickness. Default: 25.8', default=25.8
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
        '--HBDA_in', type=str,
        help='H-bond restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--NOE_in', type=str,
        help='NOE distance restraints table. Default: None', default=''
    )
    parser.add_argument(
        '--method', type=str,
        help='NOE distance restraints table. Default: None', default='fold', choices=('fold','refine')
    )
    parser.add_argument(
        '--EzPot', type=bool,
        help='Use Ez potential for depth of insrtion. Default: False', default=False, choices=(True, False)
    )
    args = parser.parse_args()
    return args
    

# Read arguments
args = parse_args()
method = args.method

# Refinement parameters
structures_in = args.structures_in     # Best structure from folding step
orderedRegion = 'resid %s:%s and (name CA)' % (args.helices[0],args.helices[1])


nstructures = args.nstructures         # Number of structures to calculate
outname = "SCRIPT_STRUCTURE.sa"        # File prefix same as python script
seed = 3421                            # Random seed

# Initial and final temperatures for High T and simulated annealin
ini_temp = 3000
fin_temp = 25
ste_temp = 12.5        # Initial and final temperatures for High T and simulated annealing 


# Input files
structures_in = args.structures_in
DC_NH_work_in  = args.DC_NH_work_in
DC_NH_free_in = ''
CSA_N1_work_in = args.CSA_N1_work_in
CSA_N1_free_in = ''
DIHE_in = args.DIHE_in
NOE_in = args.NOE_in
HBDA_in = args.HBDA_in


# OS-ssNMR parameters
DC_NH_max = args.DC_NH_max  # Maximum dipolar coupling in kHz
CSA_N1_tensor = args.CSA_N1_tensor
CSA_N1_beta = args.CSA_N1_beta

# EzPot Potential
EzPot = args.EzPot


# EEFx Membrane parameters
eefx = args.eefx
immx_com = 'resid %s:%s and (name CA)' % (args.helices[0],args.helices[1])  # Center of mass selection for IMMx position.
immx_thickness = args.immx_thickness        # 25.4*0.75 + 27.0*0.25 = 25.8 A for DMPC/POPC bicelle
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
ini_model=structures_in[0]
if method == 'fold':
    import psfGen
    psfGen.pdbToPSF(ini_model)
    protocol.genExtendedStructure()
if method == 'refine':
    protocol.loadPDB(ini_model, deleteUnknownAtoms=True)


# Load EEFX2 topology
if eefx:
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


### DIPOLAR COUPLING RESTRAINTS ###

# Initialize DC tensor
dc_tensor = {}
from varTensorTools import create_VarTensor
oTensor = create_VarTensor('amide_NH')
oTensor.setDaMax(DC_NH_max/2)			      			
oTensor.setDa(DC_NH_max/2)
oTensor.setRh(0)
dc_tensor['amide_NH'] = oTensor

# Load working DC restraints if specified
if DC_NH_work_in:
    dc = PotList('DIPL')
    from rdcPotTools import create_RDCPot
    for (name, file, tensor, scale) in [('amide_NH', DC_NH_work_in, dc_tensor['amide_NH'], 1)]:
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
if CSA_N1_work_in:
    csa = PotList('CS')
    from csaPotTools import create_CSAPot
    for (name, file, scale) in [('amide_N',  CSA_N1_work_in, 1)]:
        CS=create_CSAPot(name=name, file=file, oTensor=csa_tensor['amide_N'])
        CS.setScale(scale)
        CS.setTensorClass("bond")
        CS.setThreshold(0.01)
        CS.setPotType('square')
        CS.setSigma(CSA_N1_tensor)
        CS.setBeta(CSA_N1_beta)
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
if DIHE_in:
    protocol.initDihedrals(DIHE_in)
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

if HBDA_in:
    # hbda - distance/angle bb hbond term
    protocol.initHBDA(HBDA_in)
    pots.append( XplorPot('HBDA') )

    # hbdb - knowledge-based backbone hydrogen bond term
    protocol.initHBDB()
    pots.append( XplorPot('HBDB') )


### NOE DISTANCE RESTRAINTS

if NOE_in:
    dsts = PotList('NOE')
    from noePotTools import create_NOEPot
    for (exp, file, scale) in [('noe', NOE_in, 1)]:
        NOE = create_NOEPot(exp,file,)
        NOE.setScale(scale)
        NOE.setPotType('soft')
        dsts.append(NOE)
    pots.append(dsts)
    rampedParams.append(MultRamp(ini_noe, fin_noe, "dsts.setScale(VALUE)"))
    lowTempParams.append( StaticRamp("dsts.setScale(fin_noe )"))

    
# Ez Potential
if EzPot:
    Ez_pots = []
    from membraneTools import EzPot
    Ezt= EzPot("EZt")
    Ezt.setScale(1)
    Ezt.setXYCenter(0)
    Ezt.setThickness(25)
    Ez_pots.append(Ezt)

# Initialize IMMx membrane
if eefx:
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
if eefx:
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
    
else:
    # Standard settings using torsionDB
    pots.append(XplorPot('VDW'))          # dflt scale [1]
    highTempParams.append(StaticRamp("""protocol.initNBond(cutnb=100,repel=1.2,nbxmod=4,tolerance=45,onlyCA=1)"""))

    # Standard all-atom repel for simulated annealing
    rampedParams.append(StaticRamp("protocol.initNBond(nbxmod=4)"))
    rampedParams.append(MultRamp(ini_vrep,fin_vrep,"xplor.command('param nbonds repel VALUE end end')"))
    
    # Low temperature dynamics
    lowTempParams.append(StaticRamp("protocol.initNBond(nbxmod=4)"))
    lowTempParams.append(StaticRamp(fin_vrep,"xplor.command('param nbonds repel VALUE end end')"))
    

# Set up Internal Variables Module
from ivm import IVM
dyn = IVM()                     # IVM object for torsion-angle dynamics.
dyn.reset()                     # reset ivm topology for torsion-angle dynamics.
protocol.torsionTopology(dyn)
minc = IVM()                    # IVM object for Cartesian minimization.
protocol.cartesianTopology(minc)

from simulationTools import AnnealIVM
cool = AnnealIVM(initTemp=ini_temp,     # Cooling loop.
                 finalTemp=fin_temp,
                 tempStep=ste_temp,
                 ivm=dyn,
                 rampedParams=rampedParams)


### ACCEPTANCE CRITERIA - VEGLIA ###

def accept(pots):
    """
    return True if current structure meets acceptance criteria
    """
    if pots['DIPL']:
        if pots['DIPL'].violations() > 3:
            return False
    if pots['CS']:
        if pots['CS'].violations() > 3:
            return False
    if pots['CDIH']:
        if pots['CDIH'].violations() > 2:
            return False
    if pots['BOND'].rms() > 0.1:
        return False
    if pots['ANGL'].rms() > 2:
        return False
    if pots['IMPR'].rms() > 4:
        return False
    return True


# Calculate structure module for folding - Veglia protocol
def calcOneStructure_fold(loopInfo):
    """
    This function calculates a single structure, performs analysis on the
    structure and then writes out a pdb file with remarks.
    """
    # Generate initial structure and minimize.
    #===========================================================================
    # Generate initial structure by randomizing torsion angles.
    # Then set torsion angles from restraints (this shortens high T dynamics).
    import monteCarlo
    monteCarlo.randomizeTorsions(dyn)
    protocol.fixupCovalentGeom(maxIters=100, useVDW=1)
    if DIHE_in:
        import torsionTools
        torsionTools.setTorsionsFromTable(DIHE_in)
        
    # Import dynamics toolset  
    from protocol import initDynamics
    from simulationTools import AnnealIVM
    ini_dyn  = IVM()
    sim_cool = IVM()
    fin_dyn  = IVM()
    protocol.torsionTopology(ini_dyn)
    protocol.torsionTopology(fin_dyn)
    protocol.torsionTopology(sim_cool)

    # Temperature settings
    ini_t  = 3500   #-------- initial temp -----#
    fin_t  = 25	    #-------- final temp -------#
    stp_t  = 12.5   #-------- step size --------#

    # Stage 1 - Initial torsion minimization
    InitialParams(rampedParams)	       # Initialize parameters for high temp dynamics
    if eefx:
        setCenter(immx_com, Zpos)       # Translate selected center of mass to IMMx Zpos
        InitialParams(highTempParams1) # Reset some rampedParams.
    else:
        InitialParams(highTempParams)  # Reset some rampedParams. No EEFX.

    # Stage 1 - Torsion angle minimization.
    protocol.initMinimize(ini_dyn,
                          potList=pots,
                          numSteps=100,
                          printInterval=50)
    ini_dyn.run()
    
    # Stage 2 - High temperature torsion dynamics
    if eefx:
        # Phase in EEFX step 1
        protocol.initDynamics(dyn,
                              potList=pots,         # potential terms to use.
                              bathTemp=ini_t,       # set bath temperature.
                              initVelocities=1,     # uniform initial velocities.
                              finalTime=3,          # run for finalTime or
                              numSteps=3001,        # numSteps * 0.001, whichever is less.
                              printInterval=100)    # printing rate in steps.
        dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
        ini_dyn.run()

        # Phase in EEFX step 2
        InitialParams(highTempParams2)    
        setCenter(immx_com, Zpos)            # translate selected center of mass to IMMx Zpos.
        protocol.initDynamics(dyn,
                              potList=pots,         # potential terms to use.
                              bathTemp=ini_t,       # set bath temperature.
                              initVelocities=1,     # uniform initial velocities.
                              finalTime=3,          # run for finalTime or
                              numSteps=3001,        # numSteps * 0.001, whichever is less.
                              printInterval=100)    # printing rate in steps.
        dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
        ini_dyn.run()

        
    # High temperature dynamics. With or without EEFX.
    if eefx:
        setCenter(immx_com, Zpos)             # translate selected center of mass to IMMx Zpos.
        
    InitialParams(highTempParams)
    protocol.initDynamics(dyn,
                          potList=pots,      # potential terms to use.
                          bathTemp=ini_t,    # set bath temperature.
                          initVelocities=1,  # uniform initial velocities.
                          finalTime=20,	     # run for finalTime or
                          numSteps=20001,    # numSteps * 0.001, whichever is less.
                          printInterval=500) # printing rate in steps.
    dyn.setETolerance(ini_temp/100)          # used to det. stepsize, dflt [temp/1000].
    ini_dyn.run()

    
    # Stage 3 - Simulated annealing
    if eefx:
        setCenter(immx_com, Zpos)             # translate selected center of mass to IMMx Zpos.
        
    InitialParams(rampedParams)
    cool = AnnealIVM(ivm           = sim_cool,
 		     initTemp      = ini_t,
                     finalTemp     = fin_t,
                     tempStep      = stp_t,
		     rampedParams  = rampedParams)
    protocol.initDynamics(ivm	        = sim_cool,
			  potList	= pots,
			  finalTime	= 0.2,
			  numSteps	= 201,
			  printInterval = 100)
    cool.run()
    
    # Stage 3b - Optional Depth of insertion (EzPot)
    if EzPot:
        from xplor import select
        m = IVM(xplor.simulation)
        protocol.initMinimize(m,numSteps=10)
        #m.setVerbose(m.verbose() | m.printNodeDef)
        m.setStepType("MinimizeCG")
        m.setNumSteps(10)
        m.setDEpred(1)
        m.setETolerance(1e-7)
        m.setPrintInterval(1)
        groupList=m.groupList()
        groupList.append(select('resid 0:35'))
        m.setGroupList(groupList)
        m.setHingeList([('translate', select('resid 0:35')),])
        m.potList().removeAll()
        m.potList().add(Ezt)
        m.run()

    # Stage 4 - Low temperature torsion dynamics
    InitialParams(lowTempParams)
    protocol.initDynamics(ivm          = fin_dyn,
                          potList        = pots, 
                          bathTemp       = fin_t,
                          initVelocities = 1,
                          finalTime      = 20,           
                          numSteps       = 15000, 	  
                          printInterval  = 3000)
    fin_dyn.setETolerance( fin_t/10 )   
    fin_dyn.run()

    # Stage 5 - Powel torsion angle minimization
    protocol.initMinimize(ivm	= fin_dyn,
			  potList	= pots,
			  printInterval	= 500)
    fin_dyn.run()
    pass
    
    
# Calculate structure module (Marassi template)
def calcOneStructure_refine(loopInfo):
    """
    This function calculates a single structure, performs analysis on the
    structure and then writes out a pdb file with remarks.
    """
    InitialParams(rampedParams)     # parameters for SA.
    InitialParams(highTempParams)   # reset some rampedParams.

    # Initial Cartesian minimization.
    #===========================================================================
    protocol.initMinimize(minc,
                          numSteps=500,         # dflt [500 steps]
                          potList=pots,
                          printInterval=50,
                          dEPred=10)
    minc.run()
    
    # High temperature dynamics.
    #===========================================================================
    #setCenter(imm_com, Zpos)            # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,         # potential terms to use
                          bathTemp=ini_temp,    # set bath temperature.
                          initVelocities=1,     # uniform initial velocities.
                          finalTime=10,         # run for finalTime ps or
                          numSteps=10001,       # numSteps * 0.001, whichever is less
                          printInterval=500)
    dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000]
    dyn.run()

    # Initialize integrator and loop for simulated annealing and run.
    #===========================================================================
    InitialParams(rampedParams)
    protocol.initDynamics(dyn,
                          potList=pots,
                          finalTime=0.2,        # run for finalTime ps or
                          numSteps=201,         # numSteps * 0.001, whichever is less
                          printInterval=100)    
    cool.run()                                  # Run cooling loop.
    
    # Final minimization.
    #===========================================================================
    # Torsion angle minimization.
    protocol.initMinimize(dyn,
                          numSteps=500,         # dflt [500 steps]
                          potList=pots,
                          printInterval=50)
    dyn.run()

    # Cartesian all-atom minimization.
    protocol.initMinimize(minc,
                          numSteps=500,         # dflt [500 steps]
                          potList=pots,
                          printInterval=50,
                          dEPred=10)
    minc.run()

    # Rotate coordinates and tensor in EEF/IMMM membrane axis frame.
    #===========================================================================
    #setCenter(IMM_com, Zpos)            # translate selected center of mass to IMMx Zpos.
    setCenterXY()           # translate protein coordinates to XY center.
    
    # Do analysis and write structure when this routine is finished.
    pass


# Loop control for folding stage
if method == 'fold':
    from simulationTools import StructureLoop, FinalParams
    StructureLoop(structLoopAction=calcOneStructure_fold,
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




# Loop control (Marassi template)
#from simulationTools import StructureLoop, FinalParams
#StructureLoop(structLoopAction=calcOneStructure_refine,
#              pdbFilesIn=structures_in,
#              numStructures=nstructures,
#              pdbTemplate=outname,
#              calcMissingStructs=True,
#              doWriteStructures=True,           # analyze and write coords after calc
#              genViolationStats=True,           # print stats file
#              averageContext=FinalParams(rampedParams),
#              averageFitSel="(%s) and name CA" %# selection for bkbn fit and rmsd [CA]
#                              orderedRegion,
#              averageCompSel="not name H*",     # selection for heavy atom rmsd
#              averageTopFraction=0.2,           # Report stats on top 20%
#              averagePotList=pots,              # Terms for stats and avg
#              #averageCrossTerms=[rdcs,csas,rdcsCross,csasCross],    # Cross correlation terms
#              averageSortPots=pots,             # Terms used to sort models
#              ).run()
