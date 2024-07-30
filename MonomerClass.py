import numpy as np

import MonomerClass


class Monomer:
    CH_len = 107 # pm (class attribute)
    CC_len = 153
    HCH_angle = 109.5

    def __init__(self, carbonPosition: list, adjacentCarbons: [MonomerClass.Monomer, MonomerClass.Monomer]):
        """
        :param list carbonPosition: list of carbon positions
        :param list of Monomer adjacentCarbons:
        """
        self.weights = [12.011, 1.008, 1.008]
        self.C1 = np.array(carbonPosition)
        self.C2obj, self.C3obj = adjacentCarbons
        self.H1, self.H2 = self.generate_hydrogen_positions(np.array(adjacentCarbons))

    def generate_hydrogen_positions(self, C2obj: MonomerClass.Monomer, C3obj: MonomerClass.Monomer) -> (np.array, np.array):
        adjC = [C2obj.C1, C3obj.C1]
        if adjC[2]:
            c1c2 = adjC[0] - self.C1
            c1c3 = adjC[1] - self.C1
            z = -(c1c2 + c1c3)
            zi = z/np.linalg.norm(z) # z axis to rotate H's off of
            y = np.cross(c1c2, c1c3)
            yi = y/np.linalg.norm(y) # y axis orthogonal to plane of carbons
            h1 = zi * np.cos(self.HCH_angle) +  yi * np.sin(self.HCH_angle) # new H: zcos(ø) + ysin(ø)
            h2 = h1 - 2 * zi * np.cos(self.HCH_angle)
            return h1 * self.CH_len + self.C1, h2 * self.CH_len + self.C1

        elif len(adjC) == 1:



        else:

        return [[], []]

    def add_C1_to_previous_carbonC2(self):
        self.C2obj.assign_C3(self)

    def assign_C3(self, C3: MonomerClass.Monomer):
        self.C3obj = C3

    def get_coordinates(self):
        return self.C1 + self.H1 + self.H2

    def get_weight(self):
        return sum(self.weights)

if __name__ == "__main__":
    CH2 = Monomer([0, 0])
