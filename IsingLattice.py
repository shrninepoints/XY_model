import numpy as np
from Lattice import Lattice

class IsingLattice(Lattice):
    def __init__(self, size, J=1, T=1):
        super().__init__(size, J, T)
        self.vortices = np.random.choice([1, -1], size=(self.size, self.size))

    def set_vortices(self, vortices):
        if not all(np.isin(vortices, [1, -1])):
            raise ValueError("Vortices values should be either +1 or -1")
        self.vortices = vortices

    # You can override other methods or add new ones specific to the IsingLattice class.
