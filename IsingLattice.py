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

    def visualize_lattice(self):
        pass
    
    def save_to_txt(self, filename):
        with open(filename, 'w') as f:
            # Save vortices
            f.write("# Vortices\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{self.vortices[i, j]} ")
                f.write("\n")
            
            # Save horizontal bonds
            f.write("# Horizontal Bonds\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{int(self.bonds_horizontal[i, j])} ")
                f.write("\n")
            
            # Save vertical bonds
            f.write("# Vertical Bonds\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{int(self.bonds_vertical[i, j])} ")
                f.write("\n")
            
            # Save parents
            f.write("# Parents\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{self.parent[i, j][0]} {self.parent[i, j][1]} ")
                f.write("\n")
            
            # Save sizes
            f.write("# Sizes\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{self.size[i, j]} ")
                f.write("\n")

    def load_from_txt(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
            
            # Load vortices
            idx = lines.index("# Vortices\n") + 1
            for i in range(self.N):
                values = list(map(int, lines[idx + i].split()))
                for j in range(self.N):
                    self.vortices[i, j] = values[j]
            
            # Load horizontal bonds
            idx = lines.index("# Horizontal Bonds\n") + 1
            for i in range(self.N):
                values = list(map(int, lines[idx + i].split()))
                for j in range(self.N):
                    self.bonds_horizontal[i, j] = bool(values[j])
            
            # Load vertical bonds
            idx = lines.index("# Vertical Bonds\n") + 1
            for i in range(self.N):
                values = list(map(int, lines[idx + i].split()))
                for j in range(self.N):
                    self.bonds_vertical[i, j] = bool(values[j])
            
            # Load parents
            idx = lines.index("# Parents\n") + 1
            for i in range(self.N):
                values = list(map(int, lines[idx + i].split()))
                for j in range(0, len(values), 2):
                    self.parent[i, j//2] = (values[j], values[j+1])
            
            # Load sizes
            idx = lines.index("# Sizes\n") + 1
            for i in range(self.N):
                values = list(map(int, lines[idx + i].split()))
                for j in range(self.N):
                    self.size[i, j] = values[j]
