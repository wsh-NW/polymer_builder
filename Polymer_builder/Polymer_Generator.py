import numpy as np
import PolymerClass as Polymer
import Utility_Functions as utl

# TODO: poly_tries is technically redundant!

poly_tries = 100

def npolymer_creator(npolymers=10, polymer_length=70, box_dimensions=[0, 40, 0, 40, 0, 50], plot=True):
    """
    Creates nPolymers and plots them
    :param int npolymers: number of Polymers
    :param int polymer_length: number of CH2 units
    :param list[int] box_dimensions: dictated by catalyst surface generation in LAMMPS_converter
    :param bool plot: toggle on and off for plotting polymers
    :return: carbon and hydrogen positions
    """
    atom_positions = []
    carbon_positions = []
    hydrogen_positions = []
    length_distribution = []
    for p in range(npolymers):
        for t in range(poly_tries):
            poly = Polymer.Polymer(polymer_length, box_dimensions, atom_positions)
            length_distribution += [poly.length_distribution]
            if poly.check_poly_build():
                atom_positions += poly.return_carbons() + poly.return_hydrogens()
                carbon_positions.append(poly.return_carbons())
                hydrogen_positions.append(poly.return_hydrogens())
                break
            if t == poly_tries - 1:
                failed_mean = np.mean(length_distribution)
                failed_stdev = np.std(length_distribution)
                print(f"unsuccessful: built {p} polymer, with a mean {failed_mean} and stdev {failed_stdev} of unsuccessful tries")
                return
    print("all polymers built successfully!")
    if plot:
        utl.plot_multiple_polymers(carbon_positions, hydrogen_positions)
    return carbon_positions, hydrogen_positions
