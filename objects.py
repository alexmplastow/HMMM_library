import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

#Courtesy of chatGPT
#* over a code block indicates the block was at least partially written by an NLP model
from pathlib import Path
import re

#TODO: write a sanity check to ensure you are not botching positions for either lipids or
	#TODO: proteins
	#TODO: I may have accidentally delcared some variable as static in a method
	#TODO: in a method, in a method

#NOTE: this looks like a liability, r is not updating as needed
	#NOTE: this will just hold my data in a structure which allows me
	#NOTE: to pull usefull information later

#*
threeLetterAA2oneLetterAA = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",

    # Common alternate/ambiguous codes seen in MD and PDB files
    "ASX": "B",   # Aspartic acid or Asparagine (D/N)
    "GLX": "Z",   # Glutamic acid or Glutamine (E/Q)
    "XLE": "J",   # Leu or Ile
    "SEC": "U",   # Selenocysteine
    "PYL": "O",   # Pyrrolysine

    # Sometimes present if protonation states are encoded:
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
}


#TODO: adapt this code to accomodate both lipid and protein
#NOTE: you have haphazardly declared each atom to have a residue
	#NOTE: but describe a residue as clasically protein
	#NOTE: so when you plug in an atom which corresponds to a lipid, your code
	#NOTE: gets kinda kooky
class atom:
	def __init__(self, MDAatom):
		self.MDAatom = MDAatom

		if MDAatom.residue.resname in list(threeLetterAA2oneLetterAA.keys()):
			
			self.residue = residue(MDAatom.residue)

		#Usefull features in case I find myself needing to map back to the
			#Particular lipid
		self.resname = MDAatom.resname
		self.resID = MDAatom.resid
		self.name = MDAatom.name
		self.index = MDAatom.index

	@property
	def r(self):
		return self.MDAatom.position



#TODO: devise a contacts constructor
#TODO: add get heavy atoms distances

class residue:
	def __init__(self, universeResidue):
		self.universeResidue = universeResidue
		self.resID = self.universeResidue.resid
		self.resnum = self.universeResidue.resnum
		#NOTE: self.aminoAcid is most certainly a string
		self.aminoAcid = self.universeResidue.resname
		self.aminoAcidLetter = threeLetterAA2oneLetterAA[self.aminoAcid]
	@property
	def r(self):
		for atom in self.universeResidue.atoms:
			if atom.name == "CA":
				r = atom.position
				break
		return r

	@property
	def atoms(self):
		return self.universeResidue.atoms

	@property
	def heavyAtoms(self):
		heavyAtoms = [atom for atom in self.universeResidue.atoms if atom.name[0] != 'H']
		return heavyAtoms
		
class protein:

	def __init__(self, universe):
		self.universe = universe

	@property
	def atoms(self):
		return [atom(MDAatom) for MDAatom in self.universe.select_atoms("protein")]

	#UPDATE this function to reflect the new class
	@property
	def residues(self):
		universeResidues = self.universe.select_atoms("protein").residues
		residues = [residue(universeResidue) 
				for universeResidue 
				in universeResidues]
		return residues

	
	@property
	def frameNum(self):
		return self.universe.trajectory.frame

	@property
	def t_inNs(self):
		return self.frameNum * 0.2

	#NOTE: this will need an update now that we've moved to a custom
		#NOTE: definition .center_of_mass() will fail
	@property
	def COM(self):
		atoms = self.universe.select_atoms("protein")

		#print("/////////////////////////////////////////")
		#print("Number of protein atoms:", atoms.n_atoms)
		#print("Total protein mass:", atoms.total_mass())
		#print("Unique masses:", np.unique(atoms.masses))
		#print("NaN coordinates:", np.isnan(atoms.positions).any())

		#print("Sample names:", atoms.names[:20])
		#print("Sample types:", atoms.types[:20])
		#print("Sample elements:", atoms.elements[:20])
		#print("/////////////////////////////////////////")

		COM = atoms.center_of_mass()
		return COM


	@property
	def geometricCenter(self):
		atoms = self.universe.select_atoms("protein")
		return atoms.center_of_geometry()


	@property
	def residuePositions(self):
		R = np.zeros((0,3))
		for residue in self.residues:
			R = np.vstack((R, residue.r))
		return R

	@property
	def maxZ(self):

		R = self.residuePositions
		
		return max(R[:, 2])
	
	@property
	def minZ(self):

		R = self.residuePositions

		return min(R[:, 2])


	def getTrajectoryIndices(self):
		return list(range(0, len(self.universe.trajectory)))

	def goto(self, frameNum):
		self.universe.trajectory[frameNum]

	def getCOM_R_z(self):
		
		trajectoryIndices = self.getTrajectoryIndices()

		COM_R_z = np.array([])
		for i in trajectoryIndices:
			self.goto(i)
			COM_R_z = np.append(COM_R_z, self.COM[-1])

		return COM_R_z

	def getTinNs(self):

		trajectoryIndices = self.getTrajectoryIndices()

		T = np.array([])
		
		for i in trajectoryIndices:
			self.goto(i)
			T = np.append(T, self.t_inNs)

		return T

	def getRMSF(self):

		protein = self.atoms.select_atoms("protein")

		align.AlignTraj(self.universe,
						self.universe,
						select="protein and name CA",
						in_memory=True).run()

		reference_coordinates = self.universe.trajectory.timeseries(
						asel=protein).mean(axis=1)

		reference = mda.Merge(protein).load_new(
						reference_coordinates[:, None, :],
						order="afc")

		align.AlignTraj(self.universe,
						reference,
						select="protein and name CA",
						in_memory=True).run()

		calphas = protein.select_atoms("name CA")

		rmsfer = RMSF(calphas).run()

		self.rmsf = rmsfer.results.rmsf
		return rmsfer.results.rmsf

	#TODO: modify this function so you can use the in-built atom object
	def getHeavyAtomMembraneDistances(self, membrane, sanityCheck = False):
		
		protein_cols = [
			f"{atom.name}_{atom.residue.aminoAcid}_resID{atom.residue.resID}"
			for atom in self.atoms
			]

		lipid_rows = [
 			f"{atom.name}_{atom.resname}_resID{atom.resID}"
			for atom in membrane.heavySurfaceAtoms
			]

		protein_xyz = np.array([atom.r for atom in self.atoms])

		lipid_xyz = np.array([atom.r for atom in membrane.heavySurfaceAtoms])

		if sanityCheck:
			self.protein_xyz = protein_xyz
			self.lipid_xyz = lipid_xyz

		distances = np.linalg.norm(
			lipid_xyz[:, None, :] - protein_xyz[None, :, :],
			axis=2
		)

		return pd.DataFrame(
				distances,
				index=lipid_rows,
				columns=protein_cols)

	def getHeavyAtomMembraneDistanceTensor(self, membrane, 
						saveFramesAsCSVs = True, 
						sanityCheck = False):

		#*
		outdir = Path("tmp")
		outdir.mkdir(exist_ok=True)

		#*
		frame_re = re.compile(r"frame(\d{4})\.csv$")


		#*
		completed_frames = []
		for path in outdir.glob("frame*.csv"):
			match = frame_re.match(path.name)
			if match:
				completed_frames.append(int(match.group(1)))

		#*
		last_completed = max(completed_frames, default=-1)
		frames_to_run = [frame for frame in self.getTrajectoryIndices()
				if frame > last_completed]

		for frame in tqdm(frames_to_run, desc = "Frames elapsed"):

			self.goto(frame)

			Δd_df = self.getHeavyAtomMembraneDistances(membrane,
								sanityCheck = sanityCheck)

			if sanityCheck:

				print("You are generating position CSVs, rememeber to \n"
					"remember to use diff *csv to verify \n"
					"your positions attribute is dynamic")
				
				pd.DataFrame(self.protein_xyz).to_csv(
				f'proteinPositionTmp/proteinFrame{frame:04d}.csv')

				pd.DataFrame(self.lipid_xyz).to_csv(
				f'lipidPositionTmp/lipidFrame{frame:04d}.csv')


			if saveFramesAsCSVs:
				Δd_df.to_csv(f'tmp/frame{frame:04d}.csv')

            #Strictly speaking, I can engineer this method to run without the membrane instance
                #;however, I still find it usefull for appreciating the overall workflow
            
                #I'd prefer my code be easily comprehended than efficient, given the code is already
                #complex, anyone can understand that I put a protein and a membrane into a method
                #to get contact data, but it would be hard to know this without combing through my 
                #scripts or MDAnalysis documentation

            #NOTE: the cutoff is not a literal cutoff for the contact, I merely need to funnel
                #the resulting "contacts" into an additional function for soft cutoff classificaiton

            #NOTE: I don't technical use it to define a contact, cutoff is more like a parameter for
                #search depth
        
        
	def getProteinMembraneAtomContacts(self, membrane, cutoff = 8, savePath = "./tmp"):
		
		mdaProtein = self.universe.select_atoms("protein")
		mdaMembrane = membrane.universe.select_atoms("resname DPPC SSM CHL1 PSM")
		
		ts = self.u.trajectory[self.frameNum]
		
		pairs, distances = capped_distance(
			protein.positions,
			membrane.positions,
			max_cutoff=cutoff,
			box=self.u.trajectory[self.frameNum],
			return_distances=True
			)
			
		proteinContactIndices = pairs[:, 0]
		membraneContactIndices = pairs[:, 1]
		distances = distances
		#NOTE: the function here utilizes the Best-Hummer_Eaton approximation
		softContacts = [functions.computeContact(d) for d in distances]
		
		return proteinContactIndices, membraneContactIndices, distances, softContacts
		
		#TODO: finish fleshing out this function
		def getProteinMembraneAtomContactsAcrossFrames(self, membrane, cutoff=8):
		
			for i in trajectoryIndices:
			
				self.goto(i)

				self.getProteinMembraneAtomContacts

	 #NOTE: I'm not defining an axis of rotation, I'm defining the extent to 
                        #Which the protein should be rotated in the three dimensions of
                        #a classical cartesian plane
                def rotate(self, θ_x, θ_y, θ_z):

                        #NOTE: there is no sense in reinventing Rodrigues' rotation
                                #NOTE: trick, so I will simply implement the rotation
                                #NOTE: in mdanalyis and redefine atom
                                #NOTE: positions accordingly
                                #NOTE: you will need to update the atoms
                                #NOTE: attribute to check if this rotation has
                                #NOTE: been performed, also, print a warning
                                #NOTE: so you don't find yourself doing this with
                                #NOTE: trajectory data
                                #NOTE: be careful to ensure the residues handle the
                                #NOTE: update too, or your code will be 
                                #NOTE: pointelessly compromised

                        #TODO: add a sanity check for atom positions ☑️
                        #TODO: add a santiy check for residue positions ☑️

                        #I'm encoding the intended rotation operations
                        if θ_x != 0:
                                d_x_bool = 1
                        else:
                                d_x_bool = 0

                        if θ_y != 0:
                                d_y_bool = 1
                        else:
                                d_y_bool = 0

                        if θ_z != 0:
                                d_z_bool = 1
                        else:
                                d_z_bool = 0


                        proteinSelection = self.universe.select_atoms("protein")

                        #Executing the intended rotations here
                        if d_x_bool == 1:
                                mda.transformations.rotate(
                                                θ_x,
                                                direction = [d_x_bool, 0, 0],
                                                proteinSelection)

                        if d_y_bool == 1:
                                mda.transformations.rotate(
                                                θ_y,
                                                direction = [0, d_y_bool, 0],
                                                proteinSelection)

                        if d_z_bool == 1:
                                mda.transformations.rotate(
                                                θ_z,
                                                direction = [0, 0, d_z_bool],
                                                proteinSelection)

                        warnings.warn("rotation methods should not be executed for trajectory"
                                        " filled universes, do not use this method"
                                        "unless simply handling a PDB"
                                        "if additional metadata is associated"
                                        "with the object, unexepected behavior will happen")




#NOTE: the lipid subtypes are DPPC, SSM, CHL1, avoid DCLE
	#NOTE: use the phosphate groups in the DPPC and SSMD and the 03 in the CHL1

class lipid:
	def __init__(self, universeLipidAtoms):
		self.universeLipidAtoms = universeLipidAtoms

		for atom in universeLipidAtoms:
			
			if atom.resname == "CHL1" and atom.name == "O3":
				
				self.representativeAtom = atom
				self.lipidType = atom.resname
				break
			
			if atom.resname == "SSM" and atom.name == "P":
			
				self.representativeAtom = atom
				self.lipidType = atom.resname
				break
			
			if atom.resname == "DPPC" and atom.name == "P":
				
				self.representativeAtom = atom
				self.lipidType = atom.resname
				break

			if atom.resname == "PSM" and atom.name == "P":
				
				self.representativeAtom = atom
				self.lipidType = atom.resname
				break

	@property
	def atoms(self):
		return [atom(mdaAtom) for mdaAtom in self.universeLipidAtoms]

	@property
	def r(self):
		return self.representativeAtom.position

	

class membrane:

	def __init__(self, universe):
		self.universe = universe

	@property
	def frameNum(self):
		return self.universe.trajectory.frame

	@property
	def t_inNs(self):
		return self.frameNum * 0.2


	#This only works under the assumption of harmonic constraints
		#On the lipids, there would not be much stopping lipid from disolving from
		#The membrane otherwise

	@property
	def atoms(self):

		DPPCatoms = self.universe.select_atoms("resname DPPC")
		SSMatoms = self.universe.select_atoms("resname SSM")
		CHL1atoms = self.universe.select_atoms("resname CHL1")
		PSMatoms = self.universe.select_atoms("resname PSM")
		
		self.DPPCatoms = DPPCatoms
		self.SSMatoms = SSMatoms
		self.CHL1atoms = CHL1atoms
		self.PSMatoms = PSMatoms
		
		atoms = self.DPPCatoms + self.SSMatoms + self.CHL1atoms + self.PSMatoms

		atoms = [atom(MDAatom) for MDAatom in atoms]

		return atoms

	@property
	def surfaceAtoms(self):

		#Manipulates the atoms property into yielding some helpful attributes
		atoms = self.atoms

		surfaceAtoms = self.SSMatoms + self.CHL1atoms + self.PSMatoms

		#NOTE: I'd rather use custom atom objects
		surfaceAtoms = [atom(MDAatom) for MDAatom in surfaceAtoms]

		return surfaceAtoms

	@property
	def heavySurfaceAtoms(self):

		heavySurfaceAtoms = []
		for atom in self.surfaceAtoms:
			if 'H' not in atom.name:
				heavySurfaceAtoms.append(atom)

		return heavySurfaceAtoms

	#NOTE: it seems I was a little overzealous with the properties
		#NOTE: most of the associated meta-data will be static now
	@property
	def lipids(self):
		
		lipids = []
		residues = self.universe.residues
		for residue in residues:
			
			if residue.resname == "DPPC" or \
			residue.resname == "SSM" or \
			residue.resname == "CHL1" or \
			residue.resname == "PSM":
			
				lipids.append(lipid(residue.atoms))

		#First checking if there are colvarIndices within the membrane object's memory
		if hasattr(self, "lipidColvarIndicesList"):

			for scopeLipid, colvarIndices in zip(lipids, self.lipidColvarIndicesList):

				scopeLipid.colvarIndices = colvarIndices

		return lipids

	@property
	def lipid_Rs(self):

		R = np.zeros((0,3))

		for lipid in self.lipids:
			R = np.vstack((R, lipid.r))

		return R
		

	def returnMembraneAtomZ_positions(self):

		R = np.array([])
		for atom in self.atoms:
			r = atom.r[-1]
			R = np.append(R, r)

		return R

	@property
	def maxZ(self):

		R = self.returnMembraneAtomZ_positions()
		
		return max(R)
	
	@property
	def minZ(self):

		R = self.returnMembraneAtomZ_positions()

		return min(R)

	def getTrajectoryIndices(self):
		return list(range(0, len(self.universe.trajectory)))
	
	def goto(self, frameNum):
		self.universe.trajectory[frameNum]

	def lipidPositionSanityCheck(self, lipidIndex = 0):
	
		X = np.array([])
		Y = np.array([])

		for frame in tqdm(self.getTrajectoryIndices()):
			self.goto(frame)
			r = self.lipid_Rs[lipidIndex]

			X = np.append(X, r[0])
			Y = np.append(Y, r[1])

		plt.plot(X, Y)
		#Latex text courtesy of chatGPT
		plt.xlabel(r"$r_y \;(\mathrm{\AA})$")
		plt.ylabel(r"$r_x \;(\mathrm{\AA})$")
		plt.title(r"$\text{Single Lipid} r(t)$")
		plt.show()

	def membranePositionSanityCheck(self,  visualize = False):
		τ = np.array([])
		Z = np.array([])
		for frame in self.getTrajectoryIndices():
			#Changing frame
			self.goto(frame)

			#Acquiring temporal data
			t = self.t_inNs
			τ = np.append(τ, t)

			#Acquiring positionData
			Z_max = self.maxZ
			Z = np.append(Z, Z_max)

		plt.title('The standard deviation should be not more than 5')
		plt.xlabel('time (ns)')
		plt.ylabel('position along the z axis (Å)')
		plt.plot(τ, Z)
		plt.show()


class solventResidue:
	def __init__(self, universeResidue):
		self.universeResidue = universeResidue
		self.resname = self.universeResidue.resname
		self.resID = self.universeResidue.resid

	@property
	def atoms(self):
		return [atom(mdaAtom) for mdaAtom in self.universeResidue.atoms]

#The common solvent name or resnames in VMD speak are SCSE and DCLE/DCLED	
class hmmmMembrane(membrane):

	def __init__(self, universe, solventResname = "SCSE"):
		super().__init__(universe)
		self.solventResname = solventResname


	#I would simply treat this as entirely separate from your lipid attributes
	@property
	def solventResidues(self):
		
		solventResidues = [
			solventResidue(mdaResidue) for mdaResidue in
			self.universe.select_atoms(f"resname {self.solventResname}").residues]

		#TODO: This is just a template, from the lipid example, it will need updating
			#TODO: make sure it is updated when the time comes
		#First checking if there are colvarIndices within the membrane object's memory
		if hasattr(self, "solventColvarIndicesList"):

			for scopeSolventResidue, colvarIndices in zip(solventResidues, 
							self.solventColvarIndicesList):

				scopeSolventResidue.colvarIndices = colvarIndices

		return solventResidues

	def solventAtoms(self):

		solventAtoms = [
			atom(mdaAtom) for mdaAtom in
			self.universe.select_atoms(f"resname {self.solventResname}")
		]
		return solventAtoms

	#DCLE
	#Charm GUI gives ssm's terminal carbons [CST, CFT] 
		# and dppc's terminal carbon [C2T,C3T]	

	#SCSE
	#My/Muyun's implementation uses [C6S, C6F] for SSM
		#dppc's terminal carbon is [C26, C36]

	def addConstraintIndicesToEachLipid(self, 
						SSMname = ['CST', 'CFT'], 
						DPPCname = ['C2T', 'C3T']):

		indexFindingBool = False
		#NOTE: I do not like using lists in this way, but to keep 
			#NOTE: the colvar indices glued to my residue, I must
			#NOTE: first generate one such instance which is logged into
			#NOTE: the membrane objects memory

		colvarIndicesList = []
		#Adding the colvar indices for each lipid
		for lipid in self.lipids:
			lipid.colvarIndices = []
			for atom in lipid.atoms:
				if atom.name == SSMname[0] or atom.name == SSMname[-1]:
					lipid.colvarIndices.append(atom.index)
					indexFindingBool = True

				elif atom.name == DPPCname[0] or atom.name == DPPCname[-1]:
					lipid.colvarIndices.append(atom.index)
					indexFindingBool = True

			colvarIndicesList.append(lipid.colvarIndices)

		self.lipidColvarIndicesList = colvarIndicesList

			



				#TODO add one for DPPC

		if not indexFindingBool:
			raise Exception(f"Not a single {SSMname[0]},{SSMname[-1]},"
					f"{DPPCname[0]}, or {DPPCname[-1]} was found")

		print('Colvar Indices attribute added to lipids')

	#NOTE: the other parameter is SCSE
	def addConstraintIndicesToEachSolventResidue(self):
		
		if self.solventResname == 'DCLE':
			atomNames = ['C1', 'CL11', 'CL12', 'C2']
		
		elif self.solventResname == 'SCSE':
			atomNames = ['C1', 'C2']

		colvarIndicesList = []
		for solventResidue in self.solventResidues:
			
			solventResidue.colvarIndices = []
			for atom in solventResidue.atoms:
				if atom.name in atomNames:
					solventResidue.colvarIndices.append(atom.index)

			colvarIndicesList.append(solventResidue.colvarIndices)

		self.solventColvarIndicesList = colvarIndicesList

		print('Colvar Indices attribute added to solventResidues')

							

	#@property
	#def residues(self):	

#TODO: modify this to reflect the change in class names 
#TODO: add a heavy atoms property to your residue
#TODO: devise an atom constructor

#TODO: debug this
	#NOTE: you might want to set this function to save the csvs


#NOTE: be mindful of the dummyAtom text or you will botch the file regen
class distanceZ:
	def __init__(self, atomNumbers, dummyAtom, axis):
		self.atomNumbers = atomNumbers
		self.dummyAtom = dummyAtom
		self.axis = axis

#NOTE: this script presumes you are passing the colvar dict into the
	#NOTE: constructor
class colvar:
	def __init__(self, blockDict):
		try:
			self.name = blockDict['name']
		except:
			pass

		try:
			self.distanceZInstance = distanceZ(blockDict['distanceZ']['main']['atomNumbers'], blockDict['distanceZ']['ref'][1:], blockDict['distanceZ']['axis'])
		except:
			pass
		
		try:
			self.upperWall = blockDict['upperWall']
		except:
			pass
		
		try:
			self.upperBoundary = blockDict['upperBoundary']
		except:
			pass
		
		try:
			self.upperwallconstant = blockDict['upperwallconstant']
		except:
			pass
		
		try:
			self.lowerWall = blockDict['lowerWall']
		except:
			pass
		
		try:
			self.lowerBoundary = blockDict['lowerBoundary']
		except:
			pass
		
		try:
			self.lowerwallconstant = blockDict['lowerwallconstant']
		except:
			pass

	def addHarmonic(self, harmonic):
		self.harmonic = harmonic
	
	#NOTE: There are two types of colvars included in the *col file
		#NOTE: this determines which kind I have loaded
	def determineType(self):
		try:
			self.upperWall
			self.upperBoundary
			self.upperwallconstant
			self.lowerWall
			self.lowerBoundary
			self.lowerwallconstant
			
			self.type = 2

		except:

			self.type = 1

		return self.type

	#Note: Static, but best treated as an attribute
	@property
	def colvarType(self):
		return self.determineType()

	def determineTypeSanityCheck(self):

		self.determineType()

		print(f"the type is: {self.type}")
		print("The other attributes are:")
		for attribute in dir(self):
			print(attribute)

	#NOTE: static
	@property
	def atomNumbersString(self):

		if self.colvarType == 1:

			atomNumbersString = (
				f"{self.distanceZInstance.atomNumbers[0]} "
				f"{self.distanceZInstance.atomNumbers[-1]}"
				)

		if self.colvarType == 2:

			atomNumbersString = (
				" ".join(str(atomNumber) for atomNumber in self.distanceZInstance.atomNumbers)
				)

		return atomNumbersString

	def constructBlock(self):

		self.determineType()

		if self.type == 1:

			blockString = (
				f"colvar {{\n"
				f"    name {self.name}\n"
				f"    distanceZ {{\n"
				f"        main {{ atomNumbers {{ {self.atomNumbersString} }} }}\n"
				f"        ref {{ dummyAtom ( {self.distanceZInstance.dummyAtom[0]}, {self.distanceZInstance.dummyAtom[1]}, {self.distanceZInstance.dummyAtom[2]} ) }}\n"
				f"        axis ({self.distanceZInstance.axis[0]}, {self.distanceZInstance.axis[1]}, {self.distanceZInstance.axis[2]})\n"
				f"    }}\n"
				f"}}\n"
				f"harmonic {{\n"
				f"    colvars {self.harmonic.colvarName}\n"
				f"    centers {self.harmonic.centers}\n"
				f"    forceConstant {self.harmonic.forceConstant}\n"
				f"}}")

			return blockString

		elif self.type == 2:

			blockString = (
				f"colvar {{\n"
				f"    name {self.name}\n"
				f"    upperWall         {self.upperWall}\n"
				f"    upperBoundary     {self.upperBoundary}\n"
				f"    upperwallconstant {self.upperwallconstant}\n"
				f"    lowerWall         {self.lowerWall}\n"
				f"    lowerBoundary     {self.lowerBoundary}\n"
				f"    lowerwallconstant {self.lowerwallconstant}\n"
				f"    distanceZ {{\n"
				f"        main {{ atomNumbers {{ {self.atomNumbersString} }} }}\n"
				f"        ref {{ dummyAtom ( {self.distanceZInstance.dummyAtom[0]}, {self.distanceZInstance.dummyAtom[1]}, {self.distanceZInstance.dummyAtom[2]} ) }}\n"
				f"        axis ({self.distanceZInstance.axis[0]}, {self.distanceZInstance.axis[1]}, {self.distanceZInstance.axis[2]})\n"
				f"    }}\n"
				f"}}"
			)

			return blockString

		else:
			raise Exception("Alex WTF!?")

				
		
			


class harmonic:
	def __init__(self, harmonicDict):
		self.colvarName = harmonicDict['colvars']
		self.centers = harmonicDict['centers']
		self.forceConstant = harmonicDict['forceConstant']

#NOTE: this class is courtesy of chatGPT..... and a complete pain in my butt
class block:
	def __init__(self, name, entries=None):
		self.name = name
		self.entries = entries if entries is not None else []

	
