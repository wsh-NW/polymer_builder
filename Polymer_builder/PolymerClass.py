import numpy as np
import matplotlib.pyplot as plt
import Utility_Functions as uf
import math

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
    attempts = 100
    overlap_cutoff = 1.5 # 1.5 # applies to internal C-C (chain formation) and intra C-C/H interactions

    """TROUBLE SHOOTING"""
    length_distribution = []

########################################################################################################################
    """BUILDING FUNCTIONS FOR CARBON CHAIN"""
    def __init__(self, Cn: int, box_boundaries: list, atom_positions=None):
        """
        Initialize polymer chain
        :param int Cn: carbon count
        :param list box_boundaries: box limits [-x, x, -y, y, -z, z]
        :param list atom_positions: all atoms within box
        """
        # initial variables
        self.blim = [box_boundaries[i] + self.overlap_cutoff * (-1) ** i for i in range(6)] # box limits with overlap cutoff
        self.boxLengths = np.array([self.blim[i + 1] - self.blim[i] for i in range(0, len(self.blim), 2)])
        self.Cn = Cn
        self.Rg = 0
        self.weight = Cn * (12.011 + 2 * 1.008)
        if atom_positions is None:
            atom_positions = []
        else:
            self.atom_positions = atom_positions

        # continuously updated lists and variables
        self.good_poly_build = True
        self.carbon_list = []
        self.hydrogen_list = []

        """BUILDING POLYMER"""
        self._build_carbon_chain()
        self._add_hydrogens()


    def _build_carbon_chain(self):
        """
        picks random point within box limits then adds carbons within an overlap_cutoff
        ** stops adding if carbon is not able to be added within max attempts
        :return null:
        """
        for retry in range(self.attempts):
            c_start_pos = np.array([np.random.uniform(self.blim[i], self.blim[i+1]) for i in range(0, 6, 2)])
            temp_c_list = [c_start_pos]
            curr_c = c_start_pos
            for n in range(self.Cn):
                add_successful = self._add_carbon(curr_c, temp_c_list)
                if (add_successful):
                    curr_c = temp_c_list[-1]
                else:
                    break
            self.length_distribution += [len(temp_c_list)]
            if len(temp_c_list) == self.Cn:
                self.carbon_list = temp_c_list
                return
        if not self.carbon_list:
            self.good_poly_build = False

        # print("carbon chain successfully built")


    def _add_carbon(self, curr_pos: np.array, temp_c_list: list[np.ndarray]):
        """
        adds carbon if it fulfills the requirements in check_newC. If max attempts reached, makes good_poly_build false
        :param curr_pos:
        :return:
        """
        for retry in range(np.floor_divide(self.attempts, 5)):
            if len(temp_c_list) == 1:
                v = uf.generate_random_unit_vector()
            else:
                v = self._generate_random_unit_vector_equil_angle(temp_c_list[-2:])
            v_CC = v * self.CC_len
            new_pos = curr_pos + v_CC
            if self._check_newC(curr_pos, new_pos):
                temp_c_list.append(new_pos)
                return True
        return False

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
        ghost_carbon = uf.generate_random_unit_vector() * self.CC_len + c1c2[-1]
        k, rotate_axis, v = uf.unit_basis_3pts(c1c2[-2:] + [ghost_carbon])
        rotate_this = c1c2[-2] - c1c2[-1] # rotates about the last Catom on the axis ortho CCC
        new_carbon = uf.rodrigues_rotation(angle_werror, rotate_this, rotate_axis)
        return uf.unit_v(new_carbon)

########################################################################################################################
    """ADDING HYDROGENS"""
    def _add_hydrogens(self):
        """
        generates a hydrogen list overlayed on the carbon backbone if carbon backbone is complete
        :return:
        """
        if not self.check_poly_build():
            return

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
        print("Polymer successfully built!")

########################################################################################################################
    """RETURN AND VALIDATION FUNCTIONS - RETURN FUNCTIONS ACCESSIBLE TO USER"""
    """ **check if Polymer was built correctly under specific constraints** """

    def check_poly_build(self):
        return self.good_poly_build

    def return_carbons(self):
        return self.carbon_list

    def return_hydrogens(self):
        return self.hydrogen_list

    def _check_newC(self, prev_c: np.array, new_c: np.array) -> bool:
        """
        checks new carbon against box boundaries, previous carbon, and all existing atoms in simulation box
        :param prev_c: 
        :param new_c: 
        :return: 
        """
        def check_inbounds() -> bool:
            x, y, z = prev_c
            x1, x2, y1, y2, z1, z2 = self.blim
            return (x1 < x < x2) and (y1 < y < y2) and (z1 < z < z2)

        """def no_prevC_overlap() -> bool: # no longer need after angle correction
            return np.linalg.norm(new_c - prev_c) > self.overlap_cutoff"""

        def check_allatoms() -> bool:
            for atom in self.atom_positions:
                dist_atom_newc = np.linalg.norm(atom - new_c)
                if dist_atom_newc < self.overlap_cutoff:
                    return False
            return True

        return check_inbounds() and check_allatoms()


########################################################################################################################
"""RUN CODE BELOW"""

if __name__ == "__main__":
    num_units = 10
    box_boundaries = [-18, 18, -18, 18, 7, 45]
    polyethylene = Polymer(num_units, box_boundaries)
    uf.plot_carbon_positions(polyethylene.return_carbons())
    uf.plot_Polymer_CH2(polyethylene.return_carbons(), polyethylene.return_hydrogens())
