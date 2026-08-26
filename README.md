# HMMM_library

A collection of objects and functions which wraps mostly MDAnalysis objects to provide trejectory awary representation of protines...sorry I'm getting a phone call

| Object           | Description                                                                                                                                                                                                             |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `atom`           | Lightweight wrapper around an MDAnalysis atom. Exposes atom identity, residue metadata, trajectory index, and current position (`r`).                                                                                   |
| `residue`        | Represents a protein residue. Provides residue identifiers, one-letter amino-acid notation, Cα position, and access to all or heavy atoms.                                                                              |
| `protein`        | Represents the protein selection in an MDAnalysis universe. Provides dynamic atom/residue collections, centers, axial bounds, trajectory navigation, RMSF analysis, and protein–membrane distance/contact calculations. |
| `lipid`          | Represents one lipid residue. Uses the phosphate atom for DPPC/SSM/PSM or O3 for CHL1 as its representative position.                                                                                                   |
| `membrane`       | Represents a DPPC/SSM/PSM/CHL1 membrane. Provides lipid and atom collections, surface/heavy-atom selections, membrane bounds, trajectory navigation, and position sanity checks.                                        |
| `solventResidue` | Represents one HMMM organic-solvent residue and exposes its constituent atoms.                                                                                                                                          |
| `hmmmMembrane`   | Extends `membrane` with SCSE or DCLE solvent handling and methods for assigning NAMD colvar atom indices to individual lipids and solvent molecules.                                                                    |
| `distanceZ`      | Stores the atom group, dummy reference point, and axis defining a NAMD `distanceZ` collective variable.                                                                                                                 |
| `harmonic`       | Stores the colvar name, restraint center, and force constant for a NAMD harmonic bias.                                                                                                                                  |
| `colvar`         | Represents a parsed NAMD colvar definition. Supports harmonic restraints and upper/lower wall restraints and can reconstruct the corresponding configuration block.                                                     |
| `block`          | Generic named container used to represent parsed configuration blocks and their entries.                                                                                                                                |
