import objects
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import MDAnalysis

def getProteinCOMandLeafletZs(proteinInstance, membraneInstance):

        T = proteinInstance.getTinNs()
        R_z = proteinInstance.getCOM_R_z()
        membraneMax_R_z = np.ones(len(T)) * membraneInstance.maxZ
        membraneMin_R_z = np.ones(len(T)) * membraneInstance.minZ
        return T, R_z, membraneMax_R_z, membraneMin_R_z

def findEuclidianDistance(a, b):
    return np.linalg.norm(a - b)

def checkFrameResidueContacts(proteinInstance, membraneInstance, contactCutoff = 8, collectLipidTypeCounts = False, determineSingleContact = False):
	lipids = membraneInstance.lipids
	residues = proteinInstance.residues

	D_ij = np.zeros((len(lipids), len(residues)))
	C_ij = np.zeros((len(lipids), len(residues)))
	residue_C_j = np.zeros((len(residues)))

	residue_C_j_DPPC = np.zeros((len(residues)))
	residue_C_j_SSM = np.zeros((len(residues)))
	residue_C_j_CHL1 = np.zeros((len(residues)))

	for i, lipid in enumerate(lipids):
		for j, residue in enumerate(residues):

			d = findEuclidianDistance(lipid.r, residue.r)
			D_ij[i][j] = d

			if d < contactCutoff:
				C_ij[i][j] = 1
				residue_C_j[j] += 1

				#NOTE: it's difficult to see, but this routine is for cases where
					#NOTE: I just want to determine if any part of the protein
					#NOTE: is interacting with any part of the membrane
				if determineSingleContact:
					return 1

				if collectLipidTypeCounts:
					if lipid.lipidType == "DPPC":
						residue_C_j_DPPC[j] += 1
					elif lipid.lipidType == "SSM":
						residue_C_j_SSM[j] += 1
					elif lipid.lipidType == "CHL1":
						residue_C_j_CHL1[j] += 1
					else:
						raise Exception("not seeing the lipid type...")

	#NOTE: one final check to determine if none of the residues are actively contacting a
		#NOTE: a lipid
	if determineSingleContact:

		#NOTE: innocent until proven guilty
		PL_contact = False
		for i, _ in enumerate(lipids):
			for j, _ in enumerate(residues):

				if C_ij[i][j] == 1:
					PL_contact == True


		if PL_contact:
			return 1
		else:
			return 0


	if collectLipidTypeCounts:
		return D_ij, C_ij, residue_C_j, residue_C_j_DPPC, residue_C_j_SSM, residue_C_j_CHL1
	else:
		return D_ij, C_ij, residue_C_j


def SanityCheckSimulationResidueContacts(proteinInstance, membraneInstance, contactCutoff = 8, dynamicSanityCheck = False, sanityCheckPlotSaveDir = "/home/alpal/projects/vacAstudies/analyticalSolutions/residueContactSaveDir"):
        trajectoryIndices = proteinInstance.getTrajectoryIndices()
        #NOTE: I am trusting proteinInstance, membraneInstance = 
                #NOTE:proteinInstance(i), membraneInstance(i)
        residueIdentifiers = [f"{residue.resnum}:{residue.aminoAcidLetter}" for residue in proteinInstance.residues]
        for i in tqdm(trajectoryIndices):
                proteinInstance.goto(i)
                D_ij, C_ij, residue_C_j = checkFrameResidueContacts(proteinInstance, membraneInstance, contactCutoff = contactCutoff)

                if dynamicSanityCheck:
                        fig, ax = plt.subplots(figsize=(10, 4))

                        hbars = ax.bar(residueIdentifiers, residue_C_j, align='center')

                        ax.set_title(f"Residue contacts at frame {i}")
                        ax.set_xlabel("Residue Identifier (amino acid + index)")
                        ax.set_ylabel("Total number of lipid contacts")


                        ax.tick_params(axis='x', labelrotation=90, labelsize=6)

                        ax.set_ylim([0, 10])

                        plt.tight_layout()
                        plt.savefig(f"{sanityCheckPlotSaveDir}/frame_{i:04d}.png", dpi=300)
                        plt.close()

#NOTE: the specimenTitle variable is integrated in the line:
	#NOTE: ax.set_title(f"Residue contacts (sampled every {sample_stride} frames)")
	#NOTE: if a specimen Title is entered, then:
	#NOTE: ax.set_title(f"{specimenTitle}: residue contacts (sampled every {sample_stride} frames)")

	
	

def getCompleteContactCount(
		proteinInstance,
		membraneInstance,
		contactCutoff = 8,
		plot = True,
		sample_stride = 10,
		lipidTypeVisualization = False,
		specimenTitle = None
):
	trajectoryIndices = proteinInstance.getTrajectoryIndices()
	sampledTrajectoryIndices = trajectoryIndices[::sample_stride]

	residues = proteinInstance.residues
	residueIdentifiers = [
		f"{res.resnum}:{res.aminoAcidLetter}"
		for res in residues
	]

	total_residue_C_j = np.zeros(len(residues))

	if lipidTypeVisualization:
		total_residue_C_j_DPPC = np.zeros(len(residues))
		total_residue_C_j_SSM = np.zeros(len(residues))
		total_residue_C_j_CHL1 = np.zeros(len(residues))

	for i in tqdm(sampledTrajectoryIndices):
		proteinInstance.goto(i)

		if not lipidTypeVisualization:
			_, _, residue_C_j = checkFrameResidueContacts(
				proteinInstance,
				membraneInstance,
				contactCutoff = contactCutoff,
				collectLipidTypeCounts = False
			)
			total_residue_C_j += residue_C_j

		else:
			_, _, _, residue_C_j_DPPC, residue_C_j_SSM, residue_C_j_CHL1 = \
				checkFrameResidueContacts(
					proteinInstance,
					membraneInstance,
					contactCutoff = contactCutoff,
					collectLipidTypeCounts = True
				)

			total_residue_C_j_DPPC += residue_C_j_DPPC
			total_residue_C_j_SSM += residue_C_j_SSM
			total_residue_C_j_CHL1 += residue_C_j_CHL1

	if plot:
		fig, ax = plt.subplots(figsize=(12, 4))

		if not lipidTypeVisualization:
			ax.bar(residueIdentifiers, total_residue_C_j)
			ax.set_ylabel("Total contacts")

		else:
			colors = {
				"DPPC": "tab:blue",
				"SSM": "tab:orange",
				"CHL1": "tab:green"
			}

			bottom = np.zeros(len(residues))

			ax.bar(
				residueIdentifiers,
				total_residue_C_j_DPPC,
				bottom = bottom,
				color = colors["DPPC"],
				label = "DPPC"
			)
			bottom += total_residue_C_j_DPPC

			ax.bar(
				residueIdentifiers,
				total_residue_C_j_SSM,
				bottom = bottom,
				color = colors["SSM"],
				label = "SSM"
			)
			bottom += total_residue_C_j_SSM

			ax.bar(
				residueIdentifiers,
				total_residue_C_j_CHL1,
				bottom = bottom,
				color = colors["CHL1"],
				label = "CHL1"
			)

			ax.set_ylabel("Contacts (stacked by lipid type)")
			ax.legend(title="Lipid Type")

		if specimenTitle == None:
			ax.set_title(f"Residue contacts (sampled every {sample_stride} frames)")
		else:
			ax.set_title(f"{specimenTitle}: residue contacts (sampled every "
					f"{sample_stride} frames)")

		ax.set_xlabel("Residue (letter:index)")
		ax.tick_params(axis='x', labelrotation=90, labelsize=6)

		plt.tight_layout()
		plt.show()

	if lipidTypeVisualization:
		return (
			total_residue_C_j_DPPC,
			total_residue_C_j_SSM,
			total_residue_C_j_CHL1
		)
	else:
		return total_residue_C_j

def removeEquilibration(topology, trajectory):

        u = MDAnalysis.Universe(topology, trajectory)

        u.select_atoms("all").write(f"{trajectory}.nonEquilibrated.dcd", frames = u.trajectory[500:])

#NOTE: this is a slight departure from previous functions which employ the use of membrane
	#NOTE: and protein objects. I think this is sensable in situations where I am making
	#NOTE: comparisons accross mutliple trajectories

	#NOTE: I do get to employ such instances all the same
def getInteractionFraction(topology, trajectory, contactCutoff = 8):
	
	u = MDAnalysis.Universe(topology, trajectory)
	
	proteinInstance = objects.protein(u)
	membraneInstance = objects.membrane(u)

	trajectoryIndices = proteinInstance.getTrajectoryIndices()
	totalFrames = len(trajectoryIndices)

	contactCount = 0

	for i in tqdm(trajectoryIndices):
		proteinInstance.goto(i)
		contactBoolean = checkFrameResidueContacts(proteinInstance, membraneInstance, contactCutoff = contactCutoff, determineSingleContact = True)
		contactCount+=contactBoolean

	return contactCount/totalFrames

def getInteractionFraction(topology, trajectory, contactCutoff=8, stride=1):

	u = MDAnalysis.Universe(topology, trajectory)

	proteinInstance = objects.protein(u)
	membraneInstance = objects.membrane(u)

	trajectoryIndices = proteinInstance.getTrajectoryIndices()

	sampledIndices = trajectoryIndices[::stride]
	totalFrames = len(sampledIndices)

	contactCount = 0

	for i in tqdm(sampledIndices):
		proteinInstance.goto(i)
		contactBoolean = checkFrameResidueContacts(proteinInstance, membraneInstance,
							contactCutoff=contactCutoff,
							determineSingleContact=True)
		contactCount += contactBoolean

	contactFraction = contactCount / totalFrames
	n = totalFrames

	return contactFraction, n

#TODO: define a function for defining contacts using the paper Emad so thoughtfully shared with you
#TODO: devise a way to plug the distance matrices into this function
#TODO: output the resulting file for analysis in R

#NOTE: in the OG paper by Best-Hummer-Eaton β is the an arbitrary value, but is typically
    #NOTE: 5 Å^-1, λ is set to 1.8 for the paper, but Hale went with 4 so......

#NOTE: ref is "Native contacts determine protein folding mechanisms in atomistic simulations"
    #NOTE: or https://www.pnas.org/doi/10.1073/pnas.1311599110

def computeContact(d):
    return 1/(1+np.exp(5*(d-4)))
import re
import pandas as pd

import colvarObjects

_TOKEN_RE = re.compile(r"[{}()]|,|[^\s{}(),]+")

#Author: chatGPT
def tokenize_colvars(text):
	tokens = []

	for line in text.splitlines():
		line = line.split("#", 1)[0].strip()

		if not line:
			continue

		tokens.extend(_TOKEN_RE.findall(line))
		tokens.append("\n")

	return tokens

#Author: chatGPT
def atom(token):
	if re.fullmatch(r"[-+]?\d+", token):
		return int(token)

	if re.fullmatch(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?", token):
		return float(token)

	return token

#Author: chatGPT
def clean_values(values):
	values = [v for v in values if v not in ("\n", ",", "(", ")")]
	values = [atom(v) for v in values]

	if len(values) == 1:
		return values[0]

	return values

#Author: chatGPT
def add_entry(d, key, value):
	if key in d:
		if not isinstance(d[key], list):
			d[key] = [d[key]]

		d[key].append(value)
	else:
		d[key] = value

#Author: chatGPT
def find_matching_brace(tokens, start):
	depth = 0

	for i in range(start, len(tokens)):
		if tokens[i] == "{":
			depth += 1
		elif tokens[i] == "}":
			depth -= 1

			if depth == 0:
				return i

	raise ValueError("Unclosed brace")

#Author: chatGPT
def block_is_value_list(tokens, start):
	end = find_matching_brace(tokens, start - 1)
	inner = [t for t in tokens[start:end] if t != "\n"]

	return inner and "{" not in inner and "}" not in inner

#Author: chatGPT
def parse_block(tokens, index):
	result = {}

	while index < len(tokens):
		token = tokens[index]

		if token == "\n":
			index += 1
			continue

		if token == "}":
			return result, index + 1

		key = token
		index += 1

		while index < len(tokens) and tokens[index] == "\n":
			index += 1

		if index < len(tokens) and tokens[index] == "{":
			index += 1

			if block_is_value_list(tokens, index):
				values = []

				while tokens[index] != "}":
					values.append(tokens[index])
					index += 1

				index += 1
				add_entry(result, key, clean_values(values))
			else:
				value, index = parse_block(tokens, index)
				add_entry(result, key, value)

		else:
			values = []

			while index < len(tokens) and tokens[index] not in ("\n", "}"):
				values.append(tokens[index])
				index += 1

			add_entry(result, key, clean_values(values))

	return result, index

#Author: chatGPT
def parse_colvars_file(text):
	tokens = tokenize_colvars(text)
	parsed = {}
	index = 0

	while index < len(tokens):
		if tokens[index] == "\n":
			index += 1
			continue

		block_type = tokens[index]
		index += 1

		while index < len(tokens) and tokens[index] == "\n":
			index += 1

		if index >= len(tokens) or tokens[index] != "{":
			raise ValueError(f"Expected '{{' after {block_type!r}; got {tokens[index]!r}")

		index += 1
		block, index = parse_block(tokens, index)

		if block_type == "colvar" and "name" in block:
			parsed.setdefault("colvar", {})[block["name"]] = block

		elif block_type == "harmonic" and "colvars" in block:
			parsed.setdefault("harmonic", {})[block["colvars"]] = block

		else:
			parsed.setdefault(block_type, []).append(block)

	return parsed

def col2Colvar(colvarFile, debug1 = False, debug2 = False):

	colvarText = open(colvarFile).read()
	CH_dict = parse_colvars_file(colvarText)

	'''
	#NOTE: output is: 
		#colvar
		#harmonic
	for key in colvarDict.keys():
		print(key)
	'''

	colvarDict = CH_dict['colvar']
	if debug1:
	
		exampleKey = list(colvarDict.keys())[0]
			
		print("////////////////////////////////////////")
		print("ColvarDict")
		print(colvarDict[exampleKey])
		print("////////////////////////////////////////")

	colvarInstances = []
	colvarNames = colvarDict.keys()
	for colvarName in colvarNames:
		colvarInstances.append(colvarObjects.colvar(colvarDict[colvarName]))

	if debug2:
		exampleColvar = colvarInstances[0]

		print("////////////////////////////////////////")
		print('colvar attributes')
		print(dir(exampleColvar))
		print("////////////////////////////////////////")

	harmonicInstances = []
	harmonicDict = CH_dict['harmonic']
	for harmonicDictKey in harmonicDict.keys():
		harmonicInstances.append(colvarObjects.harmonic(harmonicDict[harmonicDictKey]))


	for colvarInstance in colvarInstances:
		for harmonicInstance in harmonicInstances:
			if colvarInstance.name == harmonicInstance.colvarName:
				colvarInstance.addHarmonic(harmonicInstance)

	return colvarInstances


def colvarTypeSanityCheck(colvarInstances):
        
        print("***************************************")

        for colvarInstance in colvarInstances:
                
                colvarInstance.determineType()

                if colvarInstance.type == 1:
                        
                        colvarInstance.determineTypeSanityCheck()

                        break

        print("***************************************")

        for colvarInstance in colvarInstances:
                
                colvarInstance.determineType()

                if colvarInstance.type == 2:

                        colvarInstance.determineTypeSanityCheck()

                        break

        print("***************************************")

def dataStructureSanityCheck(colvarInstances):
        for colvarInstance in colvarInstances:
                colvarType = colvarInstance.determineType()
                if colvarType == 1:
                        print("********************************")
                        print(colvarInstance.harmonic.forceConstant)
                        print(type(colvarInstance.harmonic.forceConstant))
                        print("********************************")
                elif colvarType == 2:
                        print("********************************")
                        print(type(colvarInstance.lowerwallconstant))
                        print(type(colvarInstance.upperwallconstant))
                        print("********************************")
                        print(dir(colvarInstance))
                        raise Exception("Stop here")


'''
def modColvarSprings(colvarFile,
                                forceConstant,
                                upperWallConstant,
                                lowerWallConstant,
                                sanityCheck = False):

                colvarInstances = col2Colvar(colvarFile)

                for colvarInstance in colvarInstances:
                        colvarType = colvarInstance.determineType()

                        if colvarType == 1:

                                colvarInstance.harmonic.forceConstant = forceConstant

                                if sanityCheck:
                                        print("type 1 exec")

                        elif colvarType == 2:

                                colvarInstance.lowerwallconstant = lowerWallConstant
                                colvarInstance.upperwallconstant = upperWallConstant

                                if sanityCheck:
                                        print("!!!")
                                        print("type 2 exec")
                                        raise Exception("No point in proceeding")
'''

def modColvarSprings(colvarInstance,
                                forceConstant,
                                upperWallConstant,
                                lowerWallConstant):

    colvarType = colvarInstance.determineType()

    if colvarType == 1:

        colvarInstance.harmonic.forceConstant = forceConstant

    elif colvarType == 2:

        colvarInstance.lowerwallconstant = lowerWallConstant
        colvarInstance.upperwallconstant = upperWallConstant


def colvarInstances2ColvarFile(colvarInstances, outputColvarFile):

        blockStrings = []

        for colvarInstance in colvarInstances:

                blockStrings.append(colvarInstance.constructBlock())

        with open(outputColvarFile, "w") as outputFile:
                outputFile.write("\n".join(blockStrings) + "\n")


def modifyColvarParametersWithCSVParameters(springParametersCSV, colvarFile, sanityCheck = False):


	df = pd.read_csv(springParametersCSV)
	colvarInstances = col2Colvar(colvarFile)
        
	for i in df.index:
		forceConstant = df.at[i, 'forceConstant']
		upperWallConstant = df.at[i, 'upperWallConstant']
		lowerWallConstant = df.at[i, 'lowerWallConstant']

		for colvarInstance in colvarInstances:

			modColvarSprings(colvarInstance,
					forceConstant,
					upperWallConstant,
					lowerWallConstant)

		colvarInstances2ColvarFile(colvarInstances, 
			f"membrane_hmmm_restraint_fixed_{forceConstant}_"
			f"{upperWallConstant}_{lowerWallConstant}.namd.col")



