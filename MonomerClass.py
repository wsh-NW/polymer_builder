import numpy as np

class Monomer:
    CH_len = 107 # pm (class attribute)
    CC_len = 153
    HCH_angle = 109.5

    def __init__(self, carbonPosition: list, hydrogenPositions: list):
        """
        simple monomer class with Carbon and hydrogen positions
        """
        self.weights = [12.011, 1.008, 1.008]
        self.C1 = np.array(carbonPosition)
        self.H1, self.H2 = hydrogenPositions

    def return_C1


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


if __name__ == "__main__":
    CH2 = Monomer([0, 0])
