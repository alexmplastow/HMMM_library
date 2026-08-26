########################################################
# Imports
########################################################

import sys
sys.path.append('/home/alpal/projects/vacAstudies/colvarCopycat/scripts/HMMM_library')

import MDAnalysis as mda

import objects
import functions

########################################################
# Parameters
########################################################

charmmSSMexampleColvar = "/home/alpal/projects/vacAstudies/colvarCopycat/charmmColvarFiles/membrane_hmmm_restraint.namd.monoLipid.col"

DCLEsystem = "/home/alpal/projects/vacAstudies/colvarCopycat/structureData/DCLEhmmm/ionLPure_rProtMembrane_0Rot.pdb"

SCSEsystem = "/home/alpal/projects/vacAstudies/colvarCopycat/structureData/SCSEhmmm/hmmm-100SSM.pdb"

SCSEcolvarOutputFile = "/home/alpal/projects/vacAstudies/colvarCopycat/scripts/hmmm-100SSM.col"

########################################################
# Main
########################################################

exampleColvars = functions.col2Colvar(charmmSSMexampleColvar)

dcleUniverse = mda.Universe(DCLEsystem)

dcleHmmmMembraneInstance = objects.hmmmMembrane(dcleUniverse, solventResname = "DCLE")

scseUniverse = mda.Universe(SCSEsystem)

scseHmmmMembraneInstance = objects.hmmmMembrane(scseUniverse, solventResname = "SCSE")

scseColvars = functions.translateCharmmColvarFileToSCSEColvars(exampleColvars, scseHmmmMembraneInstance)

functions.colvarInstances2ColvarFile(scseColvars, SCSEcolvarOutputFile)

########################################################
# TODOnes
########################################################

#TODO 1.0.0: review your colvar file and identify the pertinent atoms which must be constrained ✅

#The constrain is most cetainly applied at the lipid tails, the charmm GUI implementation
	#Includes the type 1 colvar on CST and CFT for SSM 
	#I have SSM's heavy atom end labeled C6F and C6S

	#For DPPC: charmm GUI labeled the constrained atoms labeled C2T and C3T
		# that is most certainly where the constraint should be placed

	#I do not believe CHL1 has any z-axis constraints because it is identical to a typical 
		#Full membrane model

	#DCLE has and upperwall bondary and such applied to C1, CL11, CL12, AND C2
		#I'd apply the rule the all heavy atoms

#TODO 1.0.1: devise a function which simply prints out he relavent indices or use an old method ✅
#functions.colvarAtomIndexSanityCheck(exampleColvars)

#TODO: 1.0.2.-1: devise a function which mines your pdb for data relavent to the colvar ✅

#functions.checkHMMMUniverseConstraintAtomsSanityCheck(hmmmMembraneInstance)

#TODO 1.0.2.0: devise a function which allows those colvar indices to be modified ✅

#You can expose the .atomNumbersString from your colvarInstance, no need for a specialized function

#TODO 1.0.2.1: devise a way to change the name (i.e. ctl24..., scse{resid})

#You can exampose the name using colvarInstance.name, again, there is no need for another modification ✅

#TODO: 1.0.4.-1/1.0.3.-1: Go ahead and generate a function which will first take the type 1 and 2
	#Templates then generate a new collective variable block for each residue in the membrane
	#Except for cholesterol of course ✅

#		
#			________________________
#			|parse_colvars_file	|
#	colvarsText --> |			| --> colvarDict
#			------------------------
#
#	You'd generally create a colvar instance with the following
#		CH_dict = parse_colvars_file(colvarText)
#		colvarDict = CH_dict['colvar']
#		colvarInstances = []
#		#colvarNames = colvarDict.keys()
#		for colvarName in colvarNames:
#			colvarInstances.append(objects.colvar(colvarDict[colvarName]))

#You can simply recycle a type 1 and type 2 instance, plug into a colvarInstance

#TODO 1.0.4/1.0.3: review the new colvar file and determine which atoms should ✅
	#TODO 1.0.4 receive mention in the colvar file

	#NOTE: I opened up the atoms of interest, they seem to match

#TODO 1.0.5: print out the modification ✅

	#NOTE: already done by 1.0.4

#TODO 1.0.6: save the modified colvars to a file ✅

	#NOTE: already done by 1.0.4


