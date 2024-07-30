import pandas as pd

"""## STARTING VARIABLES- EDIT FOR EACH UNIQUE RUN ##"""
starting_atom_type = 2 # ensures proper atom typing (e.g. {Pt: 1, Ni: 2, C:3, H:4})
atoms_list = ['C', 'H']
xyz_file = "../Old_Files/C25.xyz"
txt_file = "../Old_Files/C25_polym.txt"
num_skip_rows = 2

##################################################################################################

"""## XYZ DATA HANDLING ##"""
# read dataframe
xyz_df = pd.read_csv(xyz_file, delim_whitespace=True, skiprows=num_skip_rows, names=["atom", "x", "y", "z"])

# rearrange columns
df_cols = xyz_df.columns.tolist()
df_cols = df_cols[1:] + [df_cols[0]]
xyz_df = xyz_df[df_cols]

# adding atom num label to each C and H
atom_dict = {}
new_atom_labels = []
for atom in xyz_df["atom"]:
    if atom in atom_dict:
        atom_dict[atom] += 1
        new_atom_labels += ["# " + atom + str(atom_dict[atom])]
    else:
        atom_dict[atom] = 1
        new_atom_labels += ["# " + atom + str(atom_dict[atom])]
xyz_df["atom"] = new_atom_labels
print(xyz_df)

##################################################################################################
"""WRITING TO A TEXT FILE"""
num_atoms = xyz_df.shape[0]

with open(txt_file, 'w') as f:
    # First line is a comment line
    f.write('# {} Polymer Molecule File\n# header section\n'.format(xyz_file[:-4]))

    # Specify number of atoms and atom types
    f.write('{} atoms\n'.format(num_atoms))
    f.write('\n')

    f.write('# body section\n')

    # Coords section
    f.write('Coords\n\n')
    for index, row in xyz_df.iterrows():
        f.write(str(index + 1) + '\t')
        for col in xyz_df.columns:
            f.write('{}\t'.format(row[col]))
        f.write('\n')
    f.write('\n')

    # Types Section
    f.write('Types\n\n')
    atom_type = xyz_df['atom'].str.extract(r'(\w)', expand=False)
    for index, row in xyz_df.iterrows():
        f.write(str(index + 1) + '\t')
        for type_num in range(len(atoms_list)):
            if atom_type[index] == atoms_list[type_num]:
                f.write('{}\t# {}'.format(type_num + starting_atom_type, atoms_list[type_num]))
        f.write('\n')
    f.write('\n')

    # Charges Section
    f.write('Charges\n\n')
    for index, row in xyz_df.iterrows():
        f.write('{}\t0.000\n'.format(index + 1))
