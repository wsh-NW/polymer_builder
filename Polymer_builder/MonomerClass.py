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

if __name__ == "__main__":
    print("bulgogi")
