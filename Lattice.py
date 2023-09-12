import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import copy
from collections import defaultdict


class Lattice:
    def __init__(self, N):
        self.N = N
        self.vortices = np.array([[(np.cos(theta), np.sin(theta)) for theta in np.random.uniform(0, 2*np.pi, N)] for _ in range(N)])
        #bonds_horizontal[i.j] is between point (i,j) and (i,j+1)
        #bonds_vertical[i.j] is between point (i,j) and (i+1,j)
        self.bonds_horizontal = np.full((N, N), False)
        self.bonds_vertical = np.full((N, N), False)
        
        # Initialize parent and rank arrays for Disjoint Set
        self.parent = np.array([[(i,j) for j in range(N)] for i in range(N)])
        self.size = np.ones((N, N), dtype=int)

    def copy(self):
        return copy.deepcopy(self)

    def get_vector(self, x, y):
        return self.vortices[x % self.N, y % self.N]
    
    def get_horizontal_bond(self, x, y):
        return self.bonds_horizontal[x % self.N, y % self.N]
    
    def get_vertical_bond(self, x, y):
        return self.bonds_vertical[x % self.N, y % self.N]
    
    def set_vector(self, x, y, theta):
        self.vortices[x % self.N, y % self.N] = (np.cos(theta), np.sin(theta))
    
    def set_horizontal_bond(self, x, y, value):
        self.bonds_horizontal[x % self.N, y % self.N] = value
        if value:  # If the bond is set to True, merge the sets
            self.union(x, y, x, (y+1) % self.N)
    
    def set_vertical_bond(self, x, y, value):
        self.bonds_vertical[x % self.N, y % self.N] = value
        if value:  # If the bond is set to True, merge the sets
            self.union(x, y, (x+1) % self.N, y)
        
    def set_all_bonds(self, value):
        self.bonds_horizontal.fill(value)
        self.bonds_vertical.fill(value)

    def set_all_bonds_random(self):
            self.bonds_horizontal = np.random.choice([True, False], (self.N, self.N))
            self.bonds_vertical = np.random.choice([True, False], (self.N, self.N))

    def find(self, x, y):
        """Find the representative element of the set containing node (x, y)."""
        parent_x, parent_y = self.parent[x, y]
        if (x, y) != (parent_x, parent_y):
            # Path compression
            self.parent[x, y] = self.find(parent_x, parent_y)
        return self.parent[x, y]

    def union(self, x1, y1, x2, y2):
        """Merge the sets containing nodes (x1, y1) and (x2, y2)."""
        root1 = self.find(x1, y1)
        root2 = self.find(x2, y2)

        if root1.all() != root2.all():
            if self.size[root1[0],root1[1]] < self.size[root2[0],root2[1]]:
                self.parent[root1[0],root1[1]] = root2
                self.size[root2[0],root2[1]] += self.size[root1[0],root1[1]]
            else:
                self.parent[root2[0],root2[1]] = root1
                self.size[root1[0],root1[1]] += self.size[root2[0],root2[1]]

    def find_all_clusters(self):
        """Find all clusters in the lattice using the Disjoint Set."""
        clusters = defaultdict(list)
        for i in range(self.N):
            for j in range(self.N):
                root = self.find(i, j)
                clusters[tuple(root)].append((i, j))
        return list(clusters.values())

    def find_largest_cluster(self):
        """Find the size of the largest cluster using the Disjoint Set."""
        clusters = self.find_all_clusters()
        return max(len(cluster) for cluster in clusters) if clusters else 0

    def rotate_all_vortices(self):
        #Rotate all vortices in the lattice by the same random angle.
        # Generate a random angle between -pi and pi
        alpha = np.random.uniform(-np.pi, np.pi)

        # Calculate rotation matrix
        rotation_matrix = np.array([[np.cos(alpha), -np.sin(alpha)],
                                    [np.sin(alpha), np.cos(alpha)]])

        # Apply rotation to all vortices
        for i in range(self.N):
            for j in range(self.N):
                self.vortices[i, j] = np.dot(rotation_matrix, self.vortices[i, j])
        
    def reset_lattice(self):
        self.bonds_horizontal = np.full((self.N, self.N), False)
        self.bonds_vertical = np.full((self.N, self.N), False)
        
        # Initialize parent and rank arrays for Disjoint Set
        self.parent = np.array([[(i,j) for j in range(self.N)] for i in range(self.N)])
        self.size = np.ones((self.N, self.N), dtype=int)
                            
    def visualize_lattice(self):
        fig, ax = plt.subplots(figsize=(10, 10))
        cmap = plt.get_cmap('hsv')  # Color map for angles

        for i in range(self.N):
            for j in range(self.N):
                x, y = i, j
                angle = np.arctan2(self.vortices[i, j][1], self.vortices[i, j][0])
                color = cmap((angle + np.pi) / (2 * np.pi))  # Convert angle to [0, 1] for color mapping

                # Draw vertex
                circle = plt.Circle((x, y), 0.2, color=color)
                    
                ax.add_artist(circle)

                # Draw right bond
                linewidth = 3 if self.bonds_horizontal[i, j] else 0.5
                ax.plot([x, (x + 1) ], [y, y], color='black', linewidth=linewidth)

                # Draw down bond
                linewidth = 3 if self.bonds_vertical[i, j] else 0.5
                ax.plot([x, x], [y, (y + 1) ], color='black', linewidth=linewidth)

        ax.set_xlim(-1, self.N)
        ax.set_ylim(-1, self.N)
        ax.set_aspect('equal')
        ax.axis('off')

        # Add colorbar as legend for angles
        norm = mcolors.Normalize(vmin=-np.pi, vmax=np.pi)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='horizontal', pad=0.01)
        cbar.set_label('Angle (from -π to π)')

        plt.show()

    def save_to_txt(self, filename):
        with open(filename, 'w') as f:
            # Save vortices
            f.write("# Vortices\n")
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{self.vortices[i, j][0]} {self.vortices[i, j][1]} ")
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
                values = list(map(float, lines[idx + i].split()))
                for j in range(0, len(values), 2):
                    self.vortices[i, j//2] = (values[j], values[j+1])
            
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

if __name__ == '__main__':
    lattice = Lattice(10)
    lattice.load_from_txt("./temp.txt")

    lattice.set_all_bonds(False)
    lattice.set_vertical_bond(5,5,True)
    lattice.set_vertical_bond(5,6,True)
    lattice.set_horizontal_bond(5,5,True)
    lattice.set_horizontal_bond(3,4,True)

    print(lattice.find_largest_cluster())

    lattice.save_to_txt("./temp.txt")
    lattice.visualize_lattice()