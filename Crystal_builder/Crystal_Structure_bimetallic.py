import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def cell_generator(lattice_param, lattice_type, atom_types, create_file, system_size):
    """
    Single atom cell_generator
    :param (int) lattice_param:  lattice spacing (A) (e.g. PtNi- 3.859A)
    :param (str) lattice_type: fcc or bcc
    :param (list) atom_types: Pt3Ni=[1, 1, 1, 2], PtNi= [1, 2]
    :param (str) create_file: file_name.data
    :param (int) system_size: number of lattice centers
    :return: <np.array[[x, y, z]] atom positions
    """
    lattice_parameter = lattice_param

    scaled_basis = np.array([[1.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0]]) * lattice_parameter

    fcc_atoms = np.array([[0.0, 0.0, 0.0],
                          [0.5, 0.5, 0.0],
                          [0.5, 0.0, 0.5],
                          [0.0, 0.5, 0.5]]) * lattice_parameter

    bcc_atoms = np.array([[0.0, 0.0, 0.0],
                          [0.5, 0.5, 0.5]]) * lattice_parameter

    if lattice_type == 'fcc':
        atoms = fcc_atoms
    elif lattice_type == 'bcc':
        atoms = bcc_atoms


    positions = []
    for i in range(system_size):
        for j in range(system_size):
            for k in range(system_size):
                cell_position = np.array([i,j,k])
                cartesian_position = np.inner(cell_position.T, scaled_basis)
                for atom in atoms:
                    positions.append(cartesian_position + atom)

    # writing lattice to lammps
    with open(create_file, 'w') as fdata:
        fdata.write(create_file[:-5] + '\n\n')
        fdata.write('{} atoms\n'.format(len(positions)))
        fdata.write('{} atom types\n'.format(2))
        fdata.write('{} {} xlo xhi\n'.format(0.0, system_size * lattice_parameter))
        fdata.write('{} {} ylo yhi\n'.format(0.0, system_size * lattice_parameter))
        fdata.write('{} {} zlo zhi\n'.format(0.0, system_size * lattice_parameter))
        fdata.write('\n')
        fdata.write('Masses\n\n')
        fdata.write('1 {} # Pt in g/mol\n'.format(195.08))
        fdata.write('2 {} # Ni in g/mol\n'.format(58.69))
        fdata.write('Atoms\n')
        for i, pos in enumerate(positions):
            fdata.write('{} {} {} {} {}\n'.format(i + 1, atom_types[i % len(atom_types)], *pos))


def plot_atoms(atoms):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2], c='b', marker='o', )

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set title
    ax.set_title('3D Plot of FCC Atoms')

    plt.show()

Pt3Ni = cell_generator(3.924, 'fcc', [1, 1, 1, 2], 'Pt3Ni_FCC', 15)
PtNi = cell_generator()
# PtNi3 = cell_generator(3.924, 'fcc', [1, 1, 1, 2], 'PtNi3_FCC', 15) not much information on this
plot_atoms(np.array(PtNi3))


# Platinum Lattice Parameter (https://www.princeton.edu/~maelabs/mae324/glos324/platinum.htm#:~:text=Platinum%2C%20Pt&text=At%20room%20temperature%20the%20metal,a%20fracture%20strain%20of%200.4.)

