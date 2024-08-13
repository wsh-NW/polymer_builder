class GeneralCrystalGenerator:
    """
    General Crystal Class
    """
    lattice_spacing = 3.859
    n_unit_cells = 0

    def __init__(self, lattice_type, crystal_facet, atom_list):
        self.lattice_type = lattice_type
        self.crystal_facet = crystal_facet
        self.atom_list = atom_list
        self.crystal_posns = self.crystal_generation()

    def crystal_generation(self):
        """
        generates desired crystal positions and returns in an array
        :return: list[np.ndarray] crystal xyz positions
        """
        # 1) build new rotated basis
        #   - axis = cross_product(v normal to facet plane,
        return []

    def box_boundaries(self):

    def create_dictionary(self):
        """
        creates dictionary to make scripting easier
        :return dict: dictionary of parameters
        """
        dictionary = {}
        dictionary["lattice type"] = self.lattice_type
        dictionary["crystal facet"] = self.crystal_facet
        dictionary["atom list"] = self.atom_list
        return dictionary


