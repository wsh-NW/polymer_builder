import numpy as np
import matplotlib.pyplot as plt
import Utility_Functions as uf
import math

# TODO BUGS: fix terminal hydrogens, understand overlap cutoff- will it be accurate with multiple polymer chains
# TODO CONTINUE: add new class to build multiple polymers and a script to print out positions for lammps
# TODO FURTHER: include Pt and make realistic for box

class Polymer:
    """
    Simple Polymer Class
    """
    # Angle Source: https://www.sciencedirect.com/science/article/pii/0371195165801989?via%3Dihub
    CH_len = 1.07  # A (class attribute)
    CC_len = 1.53
    CCC_angle = 112
    HCH_angle = 107
    tetrahedral = 109.5 # angle used for terminal hydrogens
    angle_dev = 20 # deviation from perfect CCC angle
    attempts = 50
    overlap_cutoff = 1.5

########################################################################################################################
    """INITIALIZING PARAMETERS FOR CARBON CHAIN"""
    def __init__(self, Cn: int, box_boundaries: list):
        """
        Initialize polymer chain
        :param int Cn: carbon count
        :param list box_boundaries: box limits [-x, x, -y, y, -z, z]
        """
        # initial variables
        self.blim = box_boundaries
        self.boxLengths = np.array([self.blim[i + 1] - self.blim[i] for i in range(0, len(self.blim), 2)])
        self.Cn = Cn
        self.Rg = 0
        self.weight = Cn * (12.011 + 2 * 1.008)
        # continuously updated lists
        self.carbon_list = []
        self.hydrogen_list = []
        self.CH_positions = []
        self.polymerObjects = []
        # CREATING POLYMER
        self._build_carbon_chain()
        self._add_hydrogens()

    def _build_carbon_chain(self):
        """
        picks random point within box limits then adds carbons within an overlap_cutoff
        :return null:
        """
        C_start_pos = np.random.rand(3) * self.boxLengths + np.array([self.blim[0], self.blim[2], self.blim[4]])
        self.carbon_list.append(C_start_pos)
        curr_C = C_start_pos
        for n in range(self.Cn):
            if not self._add_carbon(curr_C):
                print(f"Max attempts to add C{n} reached. Please try again")
                return
            else:
                curr_C = self.carbon_list[-1]

        print("carbon chain successfully built")


    def _add_carbon(self, curr_pos: np.array):
        for retry in range(self.attempts):
            if len(self.carbon_list) == 1:
                v = uf.generate_random_unit_vector()
            else:
                v = self._generate_random_unit_vector_equil_angle(self.carbon_list[-2:])
            v_CC = v * self.CC_len
            new_pos = curr_pos + v_CC
            if self._check_newC(curr_pos, new_pos):
                self.carbon_list.append(new_pos)
                return True
        return False

    def _check_newC(self, curr_pos: np.array, new_pos: np.array) -> bool:
        def check_inbounds() -> bool:
            x, y, z = curr_pos
            x1, x2, y1, y2, z1, z2 = self.blim
            return (x1 < x < x2) and (y1 < y < y2) and (z1 < z < z2)

        def no_prevC_overlap() -> bool:
            return np.linalg.norm(new_pos - curr_pos) > self.overlap_cutoff

        return check_inbounds() and no_prevC_overlap()

########################################################################################################################
    """ADDING HYDROGENS"""
    def _add_hydrogens(self):
        """
        generates a hydrogen list overlayed on the carbon backbone
        :return:
        """
        def _add3h(ccc: list[np.ndarray]) -> list[np.ndarray]:
            """
            generates three hydrogens at end or beginning of chain with least sterics
            ** to generate hydrogen with least sterics, uses c1c2c3 to generate ortho ccc vector w
            :param ccc: first last three carbons starting with terminal carbon
            :return: list of hydrogens
            """
            c1, c2, c3 = ccc
            ccc_point, rot_axis, v = uf.unit_basis_3pts(ccc)
            unit_h1 = uf.unit_v(uf.rodrigues_rotation(self.tetrahedral, c2 - c1, rot_axis))
            if np.dot(unit_h1, ccc_point) < 0: # test if it is correct orientation, otherwise flip
                unit_h1 = -unit_h1
            h1 = unit_h1 * self.CH_len + c1
            # other hydrogens rotate h1 around c1-c2 bond
            h2 = uf.unit_v(uf.rodrigues_rotation(self.tetrahedral, unit_h1, c2 - c1)) * self.CH_len + c1
            h3 = uf.unit_v(uf.rodrigues_rotation(-self.tetrahedral, unit_h1, c2 - c1)) * self.CH_len + c1
            return [h1, h2, h3]

        def _add2h(c_list: list[np.ndarray]) -> list[np.ndarray]:
            c1, c2, c3 = c_list
            ccc_arrow, w, rot_axis = uf.unit_basis_3pts(c_list)
            unit_h1 = uf.unit_v(uf.rodrigues_rotation(self.HCH_angle / 2, ccc_arrow, rot_axis))
            h1 = unit_h1 * self.CH_len + c_list[1]
            unit_h2 = uf.unit_v(uf.rodrigues_rotation(-self.HCH_angle / 2, ccc_arrow, rot_axis))
            h2 = unit_h2 * self.CH_len + c_list[1]
            return [h1, h2]

        hydrogen_list = []
        for i, carbon in enumerate(self.carbon_list):
            if i == 0:
                c_list = self.carbon_list[0:3]
                hydrogen_list += _add3h(c_list)
            elif i == len(self.carbon_list) - 1:
                c_list = self.carbon_list[-1:-4:-1]
                hydrogen_list += _add3h(c_list)
            else:
                c_list = self.carbon_list[i - 1: i + 2]
                hydrogen_list += _add2h(c_list)

        self.hydrogen_list = hydrogen_list
        print("hydrogens successfully added")


########################################################################################################################
    """UTILITY FUNCTIONS- NONACCESSIBLE"""

    def _generate_random_unit_vector_equil_angle(self, c1c2: list[np.ndarray]) -> np.ndarray:
        """
        generates unit vector that satisfies equilibrium angle +∆error away from current C-C
        ** start with a ghost carbon to generate random w_hat (ortho CCC)
        ** ∆error specified in class variable self.HCH_angle and angle_dev
        -> then rotate original vector along that plane with an equilibrium angle + ∆error
        Rotation Source: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
        :param c1c2: carbons used to generate third carbon unit vector. Rotation around C2
        :return: random unit vector within angle_error range
        """
        angle_werror = np.random.uniform(self.CCC_angle - self.angle_dev, self.CCC_angle + self.angle_dev)
        ghost_carbon = uf.generate_random_unit_vector() * self.CC_len + self.carbon_list[-1]
        k, rotate_axis, v = uf.unit_basis_3pts(self.carbon_list[-2:] + [ghost_carbon])
        rotate_this = c1c2[-2] - c1c2[-1] # rotates about the last Catom on the axis ortho CCC
        new_carbon = uf.rodrigues_rotation(angle_werror, rotate_this, rotate_axis)
        return uf.unit_v(new_carbon)

########################################################################################################################
    """RETURNING FUNCTIONS- ACCESSIBLE TO USER"""

    def plot_carbon_positions(self) -> list[np.array]:
        """
        returns n list of carbon positions: array[x,y,z]
        :return: list of 3D np.arrays
        """
        array3d = np.array(self.carbon_list)
        x_label, y_label, z_label = "x (A)", "y (A)", "z (A)"
        x, y, z = array3d[:, 0], array3d[:, 1], array3d[:, 2]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x, y, z, marker='o')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_zlabel(z_label)
        plt.title("carbon backbone")
        plt.show()
        return self.carbon_list

    def plot_Polymer_CH2(self):
        """
        plots carbons and hydrogens
        :return: none
        """
        # carbon backbone plotted
        x_label, y_label, z_label = "x (A)", "y (A)", "z (A)"
        x, y, z = uf.convert_xyz(self.carbon_list)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x, y, z, marker='o')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_zlabel(z_label)
        plt.title("CH2 Backbone")

        # hydrogens plotting (units of CH2 and CH3 connected)
        h_index = 3 # start on
        for i, c in enumerate(self.carbon_list):
            if i == 0:
                for h in self.hydrogen_list[0:3]:
                    x, y, z = uf.convert_xyz([h] + [c])
                    plt.plot(x, y, z, marker="o")
            elif i == len(self.carbon_list) - 1:
                for h in self.hydrogen_list[-3:]:
                    x, y, z = uf.convert_xyz([h] + [c])
                    plt.plot(x, y, z, marker="o")
            else:
                x, y, z = uf.convert_xyz([self.hydrogen_list[h_index], c, self.hydrogen_list[h_index + 1]])
                plt.plot(x, y, z, marker="o")
                h_index += 2

        plt.show()



########################################################################################################################
"""RUN CODE BELOW"""

if __name__ == "__main__":
    num_units = 10
    box_boundaries = [-18, 18, -18, 18, 7, 45]
    polyethylene = Polymer(num_units, box_boundaries)
    # c_posns = polyethylene.return_and_plot_carbon_positions()
    polyethylene.plot_carbon_positions()
    polyethylene.plot_Polymer_CH2()
