# This script is meant to write all the output we need.

def write_topology_atoms(pep_atoms):
    # This function is used when creating a new topology (from a peptide, as example).
    # The outputs should be pasted into the topology file.

    file = open("output/topology_atoms", "w")
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(pep_atoms.to_string(index = False)))
    file.close()


def write_atomtypes_atp(pep_atomtypes):
    # This function is used to create the atomtypes.atp during the preparation of the peptide.
    # Save into the peptide_ff
    file = open("output/ff_peptide/atomtypes.atp", "w")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()
    # Save into the fibril_ff
    file = open("output/ff_fibril/atomtypes.atp", "w")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()
    # Save into the merge_ff
    file = open("output/ff_merge/atomtypes.atp", "w")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()


def write_topology_bonds(top_bonds):
    # This function creates the bonds file to paste into the topology
    file = open("output/topology_bonds", "w")
    file.write("[ bonds ]")
    file.write("\n")
    file.write(str(top_bonds.to_string(index = False)))
    file.close()


def write_topology_angles(top_angles):
    # This function creates the angles file to paste into the topology
    file = open("output/topology_angles", "w")
    file.write("[ angles ]")
    file.write("\n")
    file.write(str(top_angles.to_string(index = False)))
    file.close()


def write_topology_dihedrals(top_dihedrals):
    # This function creates the dihedrals file to paste into the topology
    file = open("output/topology_dihedrals", "w")
    file.write("[ dihedrals ]")
    file.write("\n")
    file.write(str(top_dihedrals.to_string(index = False)))
    file.close()


def write_pep_ffbonded(bonds, angles, dihedrals):
    # This function creates the peptide ffnonbonded
    file = open("output/ff_peptide/ffbonded.itp", "w")
    file.write("[ bondtypes ]\n")
    file.write(str(bonds.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ angletypes ]\n")
    file.write(str(angles.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ dihedraltypes ]\n")
    file.write(str(dihedrals.to_string(index = False)))
    file.close()


def write_pep_ffnonbonded(atomtypes, pairs):
    file = open("output/ff_peptide/ffnonbonded.itp", "w")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(pairs.to_string(index = False)))
    file.close()


def write_fib_ffbonded(bonds, angles, dihedrals):
    # This function creates the peptide ffnonbonded
    file = open("output/ff_fibril/ffbonded.itp", "w")
    file.write("[ bondtypes ]\n")
    file.write(str(bonds.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ angletypes ]\n")
    file.write(str(angles.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ dihedraltypes ]\n")
    file.write(str(dihedrals.to_string(index = False)))
    file.close()


def write_fib_ffnonbonded(atomtypes, pairs):
    file = open("output/ff_fibril/ffnonbonded.itp", "w")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(pairs.to_string(index = False)))
    file.close()


def write_merge_ffbonded(pep_ff_bonds, pep_ff_angles, merge_dihedrals):
    # This function creates the peptide ffnonbonded
    file = open("output/ff_merge/ffbonded.itp", "w")
    file.write("[ bondtypes ]\n")
    file.write(str(pep_ff_bonds.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ angletypes ]\n")
    file.write(str(pep_ff_angles.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ dihedraltypes ]\n")
    file.write(str(merge_dihedrals.to_string(index = False)))
    file.close()


def write_merge_ffnonbonded(atomtypes, merge_pairs):
    file = open("output/ff_merge/ffnonbonded.itp", "w")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(merge_pairs.to_string(index = False)))
    file.close()


def write_smog_to_gromos_dihedrals(propers_to_gro):
    # This function creates the angles file to paste into the topology
    file = open("output/smog_to_gromos_dihedrals", "w")
    file.write("[ dihedrals ]")
    file.write("\n")
    file.write(str(propers_to_gro.to_string(index = False)))
    file.close()