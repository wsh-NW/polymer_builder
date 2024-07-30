import MonomerClass as Monomer
class Polymer:
    def __init__(self, numMonomers):
        self.weight = 0
        self.Rg = 0
        self.numMonomers = numMonomers
        self.positions = []
        self.weights = []

    def build_polymer(self, numMonomers):
        for n in range(numMonomers):
            mon = Monomer()

    def get_coordinates(self):
        return []

    def auto_rotate(self):
        pass

if __name__ == "__main__":
    num_units = 10
    polyethylene = Polymer(num_units)
    coords = polyethylene.get_coordinates()
    print(f'Coordinates for {num_units} unit(s) of polyethylene:')
    print(coords)

