# SIMPLE PROTOCOL TO FOLD OR REFINE MEMBRANE PROTEIN STRUCTURES
# WITH OS-ssNMR CSA AND DC RESTRAINTS
#
# Written by D. K. Weber, Veglia Lab (last revised April 30 2021)
#
# DESCRIPTION
# -----------
# The XPLOR-NIH script implements the hybrid method (Shi...Veglia, 2009, DOI: 10.1007/s10858-009-9328-9)
# with several updates to utilize new functions and force fields now in XPLOR-NIH. This script is designed
# to be run using command line arguments from a BASH script and should not require any modification. 
#
# PROTOCOL
# --------
# This follows the refinement protocol from Tian...Marassi, 2015, DOI: 10.1016/j.bpj.2015.06.047
# The template in the XPLOR-NIH tutorials (eefx-membrane) was used predominantly, but modified
# to used DC/CSA restraints according to Veglia lab preferences. Also updated to use RepelPot.
#
# 0. Input PDB structure. Optionally unfold at loading and prior to dynamics 
# 1. Initial torsion angle minimization (100 steps)
# 2. Center protein* to membrane then high temperature torsion dynamics with RepelPot (A K for 3000 steps)
# 3. Center protein* then high temperature torsions dynamics phasing in EEFx parameters (A K for 3000 steps)
# 4. Center protein* then high temperature torsion dynamics with only EEFx paramters (A K for B steps)
# 5. Center protein* again then simulated annealing (A K to C K in D K steps, E steps per increment)
# 6. Minimize insertion depth (Z-position) using knowledge-based Ez-Potential*
# 7. Powell torsion angle minimization (500 steps)
# 8. Powell Cartesian minimization (500 steps)
# * Optional steps

# RECOMMENDED SETTINGS
# --------------------
# Fold: Unfold = 'yes', A=3500, B=25000, C=25, D=12.5, E=201, nstructures = 512 (take top 5), 
#       resetCenter = 'yes', ezPot = 'resid 0:n' (n = last resid), rampeefx = yes
# 
# Refine: Unfold = 'no', A=315, B=10000, C=2.5, D=1.25, E=201, nstructures = 512 for each of top 5 structures from folding, 
#         resetCenter = 'no', ezPot = '', rampeefx = no


# ARGUMENT PARSER
# ===============================================================================
import argparse
def parse_args():
    parser = argparse.ArgumentParser(description='Initial folding for membrane protein structure.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--structure_in', type=str,
        help='Starting structure (PDB file)'
    )
    parser.add_argument(
        '--DC_NH_work', type=str,
        help='DC restraints table - forces applied. Default: None', default=''
    )
    parser.add_argument(
        '--DC_NH_free', type=str,
        help='DC restraints table - back-calculated without forces applied. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_work', type=str,
        help='CSA restraints table - forces applied. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_free', type=str,
        help='CSA restraints table - back-calculated without forces applied. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_gly_work', type=str,
        help='CSA restraints table for glycines - forces applied. Default: None', default=''
    )
    parser.add_argument(
        '--CSA_N1_gly_free', type=str,
        help='CSA restraints table for glycines - back-calculated without forces applied. Default: None', default=''
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
        help='IMMx Membrane thickness. Default: 25.72 (DMPC/POPC bicelle)', default=25.72
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
    parser.add_argument(
        '--seed', type=int,
        help='Random seed. Default: automatic', default=''
    )
    parser.add_argument(
        '--highTempSteps', type=int, 
        help='Number of steps used for high temperature torsion dynamics with full EEFx potential. Default: 25000 (25 ps)', default=25000
    )
    parser.add_argument(
        '--initialTemp', type=float, 
        help='Temperature for high temperature torsion dynamics and starting temperature for simulated annealing. Default: 3500.0 K', default=3500.0
    )
    parser.add_argument(
        '--finalTemp', type=float, 
        help='Final temperature for simulated annealing. Default: 25.0 K', default=25.0
    )
    parser.add_argument(
        '--stepTemp', type=float, 
        help='Temperature step for simulated annealing. Default: 12.5 K', default=12.5
    )
    parser.add_argument(
        '--annealSteps', type=int, 
        help='Number of steps per simulated annealing temperature. Default: 201', default=201
    )
    parser.add_argument(
        '--unfold', type=str, choices=['yes', 'no'], 
        help='Generate extended structure and randomize torsions. Default: yes', default='yes'
    )
    parser.add_argument(
        '--repelStart', type=str, choices=['yes', 'no'], 
        help='Implement two initial high-temperature torsion dynamics stages using VDW. Only valid if using EEFX forcefield. Default: yes', default='yes'
    )
    parser.add_argument(
        '--resetCenter', type=str, choices=['yes', 'no'], 
        help='Reset Z-coordinates of TM domain after high-temp dynamics and before simulated annealing. Default: yes', default='yes'
    )
    parser.add_argument(
        '--ezPot', type=str, 
        help='Use EzPot to position . Default: None', default=''
    )
 #   parser.add_argument(
 #       '--eefx', type=str, choices=['yes', 'no'], 
 #       help='Use EEFX forcefield. If no, then RepelPot is used. Default: yes', default='yes'
 #   )
    parser.add_argument(
        '--rampeefx', type=str, choices=['yes', 'no'], 
        help='Ramp EEFX scale during simulated annealing. Default: yes', default='yes'
    )
    parser.add_argument(
        '--relax', type=str, choices=['yes', 'no'],
        help='Do torsion dynamics after simulated annealing. Default: no', default='no'
    )
    parser.add_argument(
        '--relaxTerms', type=str, nargs='+',
        help='Perform torsion dynamics after annealing with only these restraint terms (full scale). Default: All allowed terms.',
        default='', choices=['CDIH', 'NOE', 'CS_w', 'DIPL_w', 'torsionDB', 'BOND', 'IMPR', 'ANGL', 'RAMA']
    )
    parser.add_argument(
        '--relaxSteps', type=int, 
        help='Number of steps used for torsion dynamics after simulated annealing. Default: 15000 (15 ps)', default=15000
    )
    parser.add_argument(
        '--relaxTemp', type=float, 
        help='Temperature for high torsion dynamics after simulated annealing. Default: 25.0 K', default=25.0
    )
    parser.add_argument(
        '--torsionPot', type=str, choices=['torsionDB', 'RAMA', 'none'],
        help='Apply knowledge-based torsion restraints. Default: yes', default='yes'
    )
    parser.add_argument(
        '--nonbondedPot', type=str, choices=['eefx', 'repel', 'VDW'],
        help='Nonbonded parameters to use. Default: eefx', default='eefx'
    )
    
    args = parser.parse_args()
    return args

# Read arguments
args = parse_args()


import protocol
import numpy as np
import random


# FORCE CONSTANTS FOR RESTRAINTS
# ===============================================================================
# Set force constants in kcal/mol/(err^2) - initial and final values
ini_bond  = 1.0	    ;	fin_bond  = 1.0   # Bond lengths - Marassi settings
ini_angl  = 0.4	    ;	fin_angl  = 1.0   # Bond angles - Marassi settings
ini_impr  = 0.1	    ;	fin_impr  = 1.0   # Improper torsion angles - Marassi settings
ini_dihd  = 10.0    ;	fin_dihd  = 200.0 # Experimental torsions (CDIH)
ini_noe   = 2.0	    ;	fin_noe	  = 40.0  # NOE distances - Veglia settings


# Set intial and final CSA and DC force constants - based on weighting factors to CDIH terms
# See for details: Shi et al (2009) J Biol NMR, DOI: 10.1007/s10858-009-9328-9
w_slf     = args.w_slf                       # SLF/CDIH force scaling 
w_r       = args.w_r                         # DC/CSA force scaling
ini_dc    = 0.00010 * w_slf *  w_r/(w_r+1)   # Dipolar coupling force initial
fin_dc    = fin_dihd * w_slf * w_r/(w_r+1)   # Dipolar coupling force final
ini_cs    = 0.00010 * w_slf * 1/(w_r+1)      # CSA force inital
fin_cs    = fin_dihd * w_slf * 1/(w_r+1)     # CSA force final


# INITIALIZE STRUCTURE
# ===============================================================================
structure_in = args.structure_in
if args.nonbondedPot == 'eefx':
    import eefxPotTools
    eefxPotTools.initEEFx()
ini_model=structure_in
protocol.loadPDB(ini_model, deleteUnknownAtoms=True)
#import psfGen
#psfGen.pdbToPSF(ini_model)
if args.unfold == 'yes':
    protocol.genExtendedStructure()


# INITIALIZE POTENTIALS & PARAMETERS
# ===============================================================================
# Initialize potentials and parameters for high temperature dynamics simulated annealing
pots = PotList()
from simulationTools import StaticRamp, MultRamp
from simulationTools import InitialParams
highTempParams1 = []                  # Settings for high T torsion dynamics - Stage 1 EEFX
highTempParams2 = []                  # Settings for high T torsion dynamics - Stage 2 EEFX
highTempParams = []                   # Settings for high T torsion dynamics - REPEL and Stage 3 EEFX
rampedParams = []                     # Settings for annealing stage
relaxParams = []                    # Settings for low T dynamics (i.e., Powell minimization)


# DIPOLAR COUPLING RESTRAINTS
# ===============================================================================
# Read input restraints and settings
DC_NH_max = args.DC_NH_max     # Maximum N-H dipolar coupling in kHz
DC_NH_work  = args.DC_NH_work  # Input file for N-H DC working restraints
DC_NH_free  = args.DC_NH_free  # Input file for N-H DC restraints without forces applied

# Initialize N-H DC tensor
dc_tensor = {}
from varTensorTools import create_VarTensor
oTensor = create_VarTensor('amide_NH')
oTensor.setDaMax(DC_NH_max/2)			      			
oTensor.setDa(DC_NH_max/2)
oTensor.setRh(0)
dc_tensor['amide_NH'] = oTensor

# Load working DC restraints if specified
if DC_NH_work:
    dc_w = PotList('DIPL_w')
    from rdcPotTools import create_RDCPot
    for (name, file, tensor, scale) in [('amide_NH_work', DC_NH_work, dc_tensor['amide_NH'], 1)]:
        DIPL_w=create_RDCPot(name=name, file=file, oTensor=tensor)
        DIPL_w.setScale(scale)
        DIPL_w.setShowAllRestraints(True)
        DIPL_w.setUseSign(False)
        DIPL_w.setThreshold(0.01)
        DIPL_w.setPotType('square')
        dc_w.append(DIPL_w)
        pass
    pots.append(dc_w)
    #highTempParams.append(StaticRamp("dc.setScale(ini_dc)"))
    rampedParams.append(MultRamp(ini_dc, fin_dc, "dc_w.setScale(VALUE)"))
    if 'DIPL_w' in args.relaxTerms:
        relaxParams.append(StaticRamp("dc_w.setScale(fin_dc)"))
    else:
        relaxParams.append(StaticRamp("dc_w.setScale(0)"))
        
    for t in dc_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass

# Load free DC restraints if specified
if DC_NH_free:
    dc_f = PotList('DIPL_f')
    from rdcPotTools import create_RDCPot
    for (name, file, tensor, scale) in [('amide_NH_free', DC_NH_free, dc_tensor['amide_NH'], 1)]:
        DIPL_f=create_RDCPot(name=name, file=file, oTensor=tensor)
        DIPL_f.setScale(scale)
        DIPL_f.setShowAllRestraints(True)
        DIPL_f.setUseSign(False)
        DIPL_f.setThreshold(0.01)
        DIPL_f.setPotType('square')
        dc_f.append(DIPL_f)
        pass
    pots.append(dc_f)
    #highTempParams.append(StaticRamp("dc.setScale(ini_dc)"))
    rampedParams.append(MultRamp(0, 0, "dc_f.setScale(VALUE)"))
    relaxParams.append(StaticRamp("dc_f.setScale(0)"))
    
    for t in dc_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass
    

# CSA RESTRAINTS
# ===============================================================================
# Load restraint files
CSA_N1_work = args.CSA_N1_work
CSA_N1_free = args.CSA_N1_free
CSA_N1_gly_work = args.CSA_N1_gly_work
CSA_N1_gly_free = args.CSA_N1_gly_free

# Set 15N tensor values (correct to reduced sigma notation and reorder for XPLOR-NIH input)
# Non-glycine
sxx, syy, szz = args.CSA_N1_tensor
iso = np.mean(args.CSA_N1_tensor)
CSA_N1_tensor = (iso-sxx, iso-szz, iso-syy)
sxx, syy, szz = args.CSA_N1_tensor_gly
# Glycine
iso = np.mean(args.CSA_N1_tensor_gly)
CSA_N1_tensor_gly = (iso-sxx, iso-szz, iso-syy)
CSA_N1_beta = args.CSA_N1_beta
CSA_N1_beta_gly = args.CSA_N1_beta_gly

# Initialize 15N CSA tensor
csa_tensor = {}
from varTensorTools import create_VarTensor
oTensor = create_VarTensor('amide_N')
oTensor.setDaMax(DC_NH_max)			      			
oTensor.setDa(DC_NH_max/2)
oTensor.setRh(0)
csa_tensor['amide_N'] = oTensor

# Load working CSA restraints if specified
# Add more CSA tensors here if needed.
CSA_all_work = []
if CSA_N1_work:
    CSA_all_work.append(('amide_N1_work', CSA_N1_work, CSA_N1_tensor, CSA_N1_beta, 1))
if CSA_N1_gly_work:
    CSA_all_work.append(('amide_N1_gly_work', CSA_N1_gly_work, CSA_N1_tensor_gly, CSA_N1_beta_gly, 1))
    
CSA_all_free = []
if CSA_N1_free:
    CSA_all_free.append(('amide_N1_free', CSA_N1_free, CSA_N1_tensor, CSA_N1_beta, 1))
if CSA_N1_gly_free:
    CSA_all_free.append(('amide_N1_gly_free', CSA_N1_gly_free, CSA_N1_tensor_gly, CSA_N1_beta_gly, 1))

# Apply working restraints if files specified
if CSA_all_work:
    csa_w = PotList('CS_w')
    from csaPotTools import create_CSAPot
    for (name, file, atom_sigma, ang_b, scale) in CSA_all_work:
        CS_w=create_CSAPot(name=name, file=file, oTensor=csa_tensor['amide_N'])
        CS_w.setScale(scale)
        CS_w.setTensorClass("bond")
        CS_w.setThreshold(0.01)
        CS_w.setPotType('square')
        CS_w.setSigma(atom_sigma)
        CS_w.setBeta(ang_b)
        CS_w.setAtomOrder("123")
        CS_w.setShowAllRestraints(True)
        CS_w.setVerbose(True)
        CS_w.setDaScale(DC_NH_max*-1)
        csa_w.append(CS_w)
    pots.append(csa_w)
    #highTempParams.append(StaticRamp("csa.setScale(ini_cs)"))
    rampedParams.append(MultRamp(ini_cs, fin_cs, "csa_w.setScale(VALUE)"))
    if 'CS_w' in args.relaxTerms:
        relaxParams.append(StaticRamp("csa_w.setScale(fin_cs)"))
    else:
        relaxParams.append(StaticRamp("csa_w.setScale(0)"))

    for t in csa_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass

# Apply free restraints if files specified
if CSA_all_free:
    csa_f = PotList('CS_f')
    from csaPotTools import create_CSAPot
    for (name, file, atom_sigma, ang_b, scale) in CSA_all_free:
        CS_f=create_CSAPot(name=name, file=file, oTensor=csa_tensor['amide_N'])
        CS_f.setScale(scale)
        CS_f.setTensorClass("bond")
        CS_f.setThreshold(0.01)
        CS_f.setPotType('square')
        CS_f.setSigma(atom_sigma)
        CS_f.setBeta(ang_b)
        CS_f.setAtomOrder("123")
        CS_f.setShowAllRestraints(True)
        CS_f.setVerbose(True)
        CS_f.setDaScale(DC_NH_max*-1)
        csa_f.append(CS_f)
    pots.append(csa_f)
    #highTempParams.append(StaticRamp("csa.setScale(ini_cs)"))
    rampedParams.append(MultRamp(0, 0, "csa_f.setScale(VALUE)"))
    relaxParams.append(StaticRamp("csa_f.setScale(0)"))

    for t in csa_tensor.values():
        t.setFreedom( "fixDa, fixRh, fixAxis" )
        pass


# EXPERIMENTAL DIHEDRAL RESTRAINTS
# ===============================================================================
# Load restraint files
DIHE = args.DIHE

# Apply restraints if files specified
if DIHE:
    protocol.initDihedrals(DIHE)
    from xplorPot import XplorPot
    pots.append(XplorPot('CDIH'))
    highTempParams1.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    highTempParams2.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    highTempParams.append(StaticRamp("pots['CDIH'].setScale(ini_dihd)"))
    rampedParams.append(StaticRamp("pots['CDIH'].setScale(fin_dihd)"))
    if 'CDIH' in args.relaxTerms:
        relaxParams.append(StaticRamp("pots['CDIH'].setScale(fin_dihd)"))
    else:
        relaxParams.append(StaticRamp("pots['CDIH'].setScale(0)"))
    pots['CDIH'].setThreshold(5.00)     # dflt [2.0]


# CHEMISTRY KNOWLEDGE TERMS
# ===============================================================================
# torsionDB torsion-angle database potential - supersedes RAMA
if args.torsionPot == 'torsionDB':
    from torsionDBPotTools import create_TorsionDBPot
    torsionDB = create_TorsionDBPot(name='torsionDB', system='protein')
    pots.append(torsionDB)
    rampedParams.append(MultRamp(0.002, 2, "torsionDB.setScale(VALUE)"))
    if 'torsionDB' in args.relaxTerms:
        relaxParams.append(StaticRamp("pots['torsionDB'].setScale(2)"))
    else:
        relaxParams.append(StaticRamp("pots['torsionDB'].setScale(0)"))

# Old RAMA torsion-angle database potential
if args.torsionPot == 'RAMA':
    protocol.initRamaDatabase()
    pots.append(XplorPot('RAMA'))
    rampedParams.append(MultRamp(.02,1,"pots['RAMA'].setScale(VALUE)"))
    if 'RAMA' in args.relaxTerms:
        relaxParams.append(StaticRamp("pots['RAMA'].setScale(1)"))
    else:
        relaxParams.append(StaticRamp("pots['RAMA'].setScale(0)"))
    
        
# Bonded geometrical parameters
from xplorPot import XplorPot
pots.append(XplorPot('BOND'))          # Dflt scale [1]
pots.append(XplorPot('ANGL'))          # Dflt scale [1]
pots.append(XplorPot('IMPR'))          # Dflt scale [1]
rampedParams.append(MultRamp(ini_bond, fin_bond, "pots['BOND'].setScale(VALUE)"))
rampedParams.append(MultRamp(ini_angl, fin_angl, "pots['ANGL'].setScale(VALUE)"))
rampedParams.append(MultRamp(ini_impr, fin_impr, "pots['IMPR'].setScale(VALUE)"))

if 'BOND' in args.relaxTerms:
    relaxParams.append(StaticRamp("pots['BOND'].setScale(fin_bond)"))
else:
    relaxParams.append(StaticRamp("pots['BOND'].setScale(0)"))
    
if 'ANGL' in args.relaxTerms:
    relaxParams.append(StaticRamp("pots['ANGL'].setScale(fin_angl)"))
else:
    relaxParams.append(StaticRamp("pots['ANGL'].setScale(0)"))
    
if 'IMPR' in args.relaxTerms:
    relaxParams.append(StaticRamp("pots['IMPR'].setScale(fin_impr)"))
else:
    relaxParams.append(StaticRamp("pots['IMPR'].setScale(0)"))

# Set threshold for violations of PotList terms
pots['BOND'].setThreshold(0.05)     # dflt [0.05]
pots['ANGL'].setThreshold(5.00)     # dflt [2.0]
pots['IMPR'].setThreshold(5.00)     # dflt [2.0]


# HYDROGEN BONDING RESTRAINTS
# ===============================================================================
# Load restraint files
HBDA = args.HBDA

# Apply restraints if files specified
if HBDA:
    # hbda - distance/angle bb hbond term
    protocol.initHBDA(HBDA)
    pots.append( XplorPot('HBDA') )
    
    # hbdb - knowledge-based backbone hydrogen bond term
    protocol.initHBDB()
    pots.append( XplorPot('HBDB') )

        
# NOE DISTANCE RESTRAINTS
# ===============================================================================
# Load restraint files
NOE = args.NOE

# Apply restraints if files specified
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
    if 'NOE' in args.relaxTerms:
        relaxParams.append( StaticRamp("dsts.setScale(fin_noe)"))
    else:
        relaxParams.append( StaticRamp("dsts.setScale(0)"))


# IMMX IMPLICIT MEMBRANE MODEL
# ===============================================================================
# Read IMMx Membrane parameters
# Membrane hydrophobic thickness: DMPC (25.4A), DPPC (28.6A), POPC (27.0A), DOPC (29.6A)
# See Marsh, Handbook of lipid bilayers, 2nd Ed. p. 379. (FROM XPLOR EXAMPLES)
immx_com = 'resid %s:%s and (name CA)' % (args.tm_domain[0],args.tm_domain[1])
immx_thickness = args.immx_thickness        # 25.4*0.8 + 27.0*0.2 = 25.72 A for DMPC/POPC bicelle ()
immx_nparameter = args.immx_nparameter      # IMMx n parameter of membrane profile

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
setCenter(immx_com, Zpos)		# Translate selected center of mass to IMMx Zpos.
# only apply restraints if eefx is specified
if args.nonbondedPot == 'eefx':
    pots.append(eefx)
    highTempParams1.append(StaticRamp("eefx.setScale(0)"))
    if args.rampeefx == 'yes':
        highTempParams.append(StaticRamp("eefx.setScale(0.004)"))
        rampedParams.append(MultRamp(0.004,1.0,"eefx.setScale(VALUE)"))
    if args.rampeefx == 'no':
        highTempParams.append(StaticRamp("eefx.setScale(1.0)"))
        rampedParams.append(StaticRamp("eefx.setScale(1.0)"))
    relaxParams.append(StaticRamp("eefx.setScale(1.0)"))
 

# INITIALIZE REPEL FOR FIRST STAGES OF TORSION DYNAMICS
#===============================================================================
if args.nonbondedPot == 'eefx':
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
    relaxParams.append(StaticRamp("pots['VDW'].setScale(0)"))


# INITIALIZE NEW REPELPOT POTENTIAL 
#===============================================================================
repelStart = args.repelStart
if args.nonbondedPot == 'repel':
    # Turn off repelStart (if accidentally set to 'yes') option
    # used to phase in EEFX params
    repelStart = 'no'
    
    from repelPotTools import create_RepelPot,initRepel
    repel = create_RepelPot('repel')
    pots.append(repel)

    if args.torsionPot == 'torsionDB':
        rampedParams.append(StaticRamp("initRepel(repel,use14=False)"))
    else:
        rampedParams.append(StaticRamp("initRepel(repel,use14=True)"))
        
    rampedParams.append(MultRamp(.004,4,"repel.setScale(VALUE)") )
    relaxParams.append(StaticRamp("repel.setScale(4)"))
    # nonbonded interaction only between CA atoms
    highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )

    # 1-4 interactions when using torsionDB
    if args.torsionPot == 'torsionDB':
        # Selected 1-4 interactions.
        import torsionDBPotTools
        repel14 = torsionDBPotTools.create_Terminal14Pot('repel14')
        pots.append(repel14)
        highTempParams.append(StaticRamp("repel14.setScale(0)"))
        rampedParams.append(MultRamp(0.004, 4, "repel14.setScale(VALUE)"))
        relaxParams.append(StaticRamp("repel14.setScale(4)"))


# INITIALIZE OLD REPEL POTENTIAL
#===============================================================================
if args.nonbondedPot == 'VDW':
    # Turn off repelStart (if accidentally set to 'yes') option.
    repelStart = 'no'
    pots.append(XplorPot('VDW'))

    if args.torsionPot == 'torsionDB':
        # Exclude 1-2, 1-3, 1-4 interactions if using torsionDB (nbxmod = 4) 
        highTempParams.append(StaticRamp("""protocol.initNBond(cutnb=100,
                                                           repel=1.2,
                                                           rcon=0.004,
                                                           nbxmod=4,
                                                           tolerance=45,
                                                           onlyCA=1)"""))
        rampedParams.append(StaticRamp("protocol.initNBond(nbxmod=4)"))

    else:
        # Exclude 1-2, 1-3. Compute 1-4 interactions if not using torsionDB (nbxmod = 5) 
        highTempParams.append(StaticRamp("""protocol.initNBond(cutnb=100,
                                                           repel=1.2,
                                                           rcon=0.004,
                                                           nbxmod=5,
                                                           tolerance=45,
                                                           onlyCA=1)"""))
        rampedParams.append(StaticRamp("protocol.initNBond(nbxmod=5)"))
        
    rampedParams.append(MultRamp(0.9, 0.8,"xplor.command('param nbonds repel VALUE end end')"))
    rampedParams.append(MultRamp(0.004, 4.0,"xplor.command('param nbonds rcon VALUE end end')"))

    
# EZ-POTENTIAL
#===============================================================================
Ez_pots = []
if args.ezPot:
    from membraneTools import EzPot
    Ezt= EzPot("EZt")
    Ezt.setScale(1)
    Ezt.setXYCenter(0)
    Ezt.setThickness(args.immx_thickness)
    Ez_pots.append(Ezt)


# SET UP IVM OBJECTS (dyn, minc) THAT PERFORM DYNAMICS AND MINIMIZATION.
#===============================================================================
# IVM (internal variable module) is used to perform dynamics and minimization in
# both torsion-angle and Cartesian space. Bonds, angles and many impropers cannot
# change with internal torsion-angle dynamics.
from ivm import IVM
dyn = IVM()                       # IVM object for torsion-angle dynamics.
dyn.reset()                       # reset ivm topology for torsion-angle dynamics.
protocol.torsionTopology(dyn)

minc = IVM()                      # IVM object for Cartesian minimization.
protocol.cartesianTopology(minc)


# DYNAMICS SETTINGS FOR HIGH T AND ANNEALING STAGES.
#===============================================================================
nstructures = args.nstructures                  # Number of structures to calculate
outname = "SCRIPT_STRUCTURE.sa"                 # File prefix same as python script
# Random seed
if args.seed:
    seed = args.seed
else:
    seed = random.randint(1,10000)

ini_temp = args.initialTemp   ; fin_temp = args.finalTemp  ; step_temp = args.stepTemp
highTempSteps = args.highTempSteps  ; annealSteps = args.annealSteps  ;
protocol.massSetup()              # Give atoms uniform weights except for axes.


# CALCULATE STRUCTURE MODULE
#===============================================================================
def calcOneStructure(loopInfo):
    """
    This function calculates a single structure, performs analysis on the
    structure and then writes out a pdb file with remarks.
    """        
    # Generate initial structure and minimize.
    #===========================================================================
    # Generate initial structure by randomizing torsion angles.
    if args.unfold == 'yes':
        import monteCarlo
        monteCarlo.randomizeTorsions(dyn)

    # Then set torsion angles from restraints (this shortens high T dynamics).
    protocol.fixupCovalentGeom(maxIters=100, useVDW=1)
    import torsionTools
    if DIHE:
        import torsionTools
        torsionTools.setTorsionsFromTable(DIHE)

    # Initialize parameters
    InitialParams(rampedParams)     # Parameters for SA.
    InitialParams(highTempParams1)  # Reset some rampedParams.

    if args.resetCenter == 'yes':
        setCenter(immx_com, Zpos)       # Translate selected center of mass to IMMx Zpos.

    # Torsion angle minimization.
    #==========================================================================
    protocol.initMinimize(dyn,
                          potList=pots,
                          numSteps=100,
                          printInterval=50)
    dyn.run()

    # High temperature dynamics.
    #===========================================================================
    # Start with REPEL to remove clashes then phase out
    if repelStart == 'yes':
        # High temperature dynamics stage 1.
        protocol.initDynamics(dyn,
                              potList=pots,         # potential terms to use.
                              bathTemp=ini_temp,    # set bath temperature.
                              initVelocities=1,     # uniform initial velocities.
                              finalTime=30,          # run for finalTime or
                              numSteps=3001,        # numSteps * 0.001, whichever is less.
                              printInterval=100)    # printing rate in steps.
        dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
        dyn.run()
    
        # High temperature dynamics stage 2.
        InitialParams(highTempParams2)
        if args.resetCenter == 'yes':
            setCenter(immx_com, Zpos)               # translate selected center of mass to IMMx Zpos.
        protocol.initDynamics(dyn,
                              potList=pots,         # potential terms to use.
                              bathTemp=ini_temp,    # set bath temperature.
                              initVelocities=1,     # uniform initial velocities.
                              finalTime=30,          # run for finalTime or
                              numSteps=3001,        # numSteps * 0.001, whichever is less.
                              printInterval=100)    # printing rate in steps.
        dyn.setETolerance(ini_temp/100)             # used to det. stepsize, dflt [temp/1000].
        dyn.run()

    # High temperature dynamics stage 3.
    InitialParams(highTempParams)
    if args.resetCenter == 'yes':
        setCenter(immx_com, Zpos)                        # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,                  # potential terms to use.
                          bathTemp=ini_temp,             # set bath temperature.
                          initVelocities=1,              # uniform initial velocities.
                          numSteps=highTempSteps,        # numSteps * 0.001, whichever is less.
                          finalTime=highTempSteps/100,	 # run for finalTime or
                          printInterval=100)             # printing rate in steps.
    dyn.setETolerance(ini_temp/100)                      # used to det. stepsize, dflt [temp/1000].
    dyn.run()

    # Initialize integrator and loop for simulated annealing and run.
    #===========================================================================
    # Dynamics for annealing.
    if args.resetCenter == 'yes':
        setCenter(immx_com, Zpos)                # translate selected center of mass to IMMx Zpos.
    protocol.initDynamics(dyn,
                          potList=pots,
                          finalTime=0.4,         # run for finalTime or
                          numSteps=annealSteps,  # numSteps*0.001, whichever is less.
                          printInterval=100)    

    # Set up cooling loop and run.
    from simulationTools import AnnealIVM
    AnnealIVM(initTemp=ini_temp,
                       finalTemp=fin_temp,
                       tempStep=step_temp,
                       ivm=dyn,
                       rampedParams=rampedParams
                       ).run()

    # Run Ez-Potential to position Z-axis
    #===========================================================================
    if args.ezPot:
        from xplor import select
        m = IVM(xplor.simulation)
        protocol.initMinimize(m,numSteps=10)
        m.setVerbose(m.verbose() | m.printNodeDef)
        m.setStepType("MinimizeCG")
        m.setNumSteps(10)
        m.setDEpred(1)
        m.setETolerance(1e-7)
        m.setPrintInterval(1)
        groupList=m.groupList()
        groupList.append(select(args.ezPot))
        m.setGroupList(groupList)
        m.setHingeList([('translate', select(args.ezPot)),])
        m.potList().removeAll()
        m.potList().add(Ezt)
        m.run()

    if args.relax == 'yes':
        InitialParams(relaxParams)
        protocol.initDynamics(dyn,
                              potList=pots,                # potential terms to use.
                              bathTemp=args.relaxTemp,     # set bath temperature.
                              initVelocities=1,            # uniform initial velocities.
                              finalTime=30,                # run for finalTime or
                              numSteps=args.relaxSteps,    # numSteps * 0.001, whichever is less.
                              printInterval=100)           # printing rate in steps.
        dyn.setETolerance(args.relaxTemp/100)                    # used to det. stepsize, dflt [temp/1000].
        dyn.run()
        
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
              averageContext=FinalParams(relaxParams),
              averageFitSel="",                 # selection for bkbn rmsd [CA]
              averageCompSel="",                # selection for heavy atom rmsd
              averageTopFraction=0.1,           # Report stats on top 10%
              averagePotList=pots,              # Terms for stats and avg
              averageSortPots=[                 # Terms used to sort models
                  term for term in pots \
                  if term.instanceName() not in ("torsionDB","VDW")],
).run()

