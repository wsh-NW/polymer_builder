import numpy as np
import PolymerClass as Polymer
import matplotlib.pyplot as plt
import Utility_Functions as utl

atom_positions = []
carbon_positions = []
hydrogen_positions = []
poly_tries = 100

def npolymer_creator(npolymers=10, polymer_length=70, box_dimensions=[0, 40, 0, 40, 0, 50]):
    length_distribution = []
    for p in range(npolymers):
        for t in range(poly_tries):
            poly = Polymer.Polymer(polymer_length, box_dimensions, atom_positions)
            length_distribution += [poly.length_distribution]
            if poly.check_poly_build():
                addPolymer(poly)
                break
            if t == poly_tries - 1:
                failed_mean = np.mean(length_distribution)
                failed_stdev = np.std(length_distribution)
                print(f"unsuccessful: built {p} polymer, with a mean {failed_mean} and stdev {failed_stdev} of unsuccessful tries")
                return
    print("all polymers built successfully!")

def addPolymer(poly: Polymer):
    """
    adds Polymer to atom_positions
    :param poly: Polymer object
    :return: NONE
    """
    global atom_positions, carbon_positions, hydrogen_positions
    atom_positions += poly.return_carbons() + poly.return_hydrogens()
    carbon_positions.append(poly.return_carbons())
    hydrogen_positions.append(poly.return_hydrogens())

########################################################################################################################
""" RUN/TESTING SECTION """

npolymer_creator()
utl.plot_multiple_polymers(carbon_positions, hydrogen_positions)