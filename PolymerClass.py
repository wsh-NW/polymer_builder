import numpy as np
import matplotlib.pyplot as plt

class Polymer:
    """
    Simple Polymer Class (not concerned about overlap yet)
    """
    CH_len = 1.07  # A (class attribute)
    CC_len = 1.53
    HCH_angle = 109.5
    attempts = 50
    overlap_cutoff = 0.5

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
        self.CH_positions = []
        self.polymerObjects = []
        # CREATING POLYMER
        self._build_carbon_chain()

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

        print("Polymer successfully built")


    def _add_carbon(self, curr_pos: np.array):
        for retry in range(self.attempts):
            v = self._generate_random_unit_vector()
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

########################################################################################################################
    """RETURNING FUNCTIONS- ACCESSIBLE TO USER"""
    def return_and_plot_carbon_positions(self) -> list[np.array]:
        """
        returns n list of carbon positions: array[x,y,z]
        :return: list of 3D np.arrays
        """
        self._plot_3d_array(self.carbon_list, "Carbon Chain Plot")
        return self.carbon_list

########################################################################################################################
    """UTILITY FUNCTIONS- NONACCESSIBLE"""
    def _generate_random_unit_vector(self) -> np.array:
        """
        generates random unit vector (r=1) with phi and theta
        Source: https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
        :return: random unit vector
        """
        phi = np.random.random() * 2 * np.pi # phi = random(0,2pi)
        costheta = np.random.random() * 2 + -1 # costheta = random(-1,1)
        theta = np.arccos(costheta)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return np.array([x, y, z])

    def _plot_3d_array(self, array3d: list[np.array], title, xaxis="x (A)", yaxis="y (A)", zaxis="z (A)"):
        points = [np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9]), np.array([10, 11, 12])]

        array3d = np.array(array3d)
        x, y, z = array3d[:, 0], array3d[:, 1], array3d[:, 2]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x, y, z, marker='o')
        ax.set_xlabel(xaxis)
        ax.set_ylabel(yaxis)
        ax.set_zlabel(zaxis)
        plt.show()

if __name__ == "__main__":
    num_units = 10
    box_boundaries = [-18, 18, -18, 18, 7, 45]
    polyethylene = Polymer(num_units, box_boundaries)
    print(f'Coordinates for {num_units} unit(s) of polyethylene:')
    c_posns = polyethylene.return_and_plot_carbon_positions()