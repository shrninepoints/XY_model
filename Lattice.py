import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import copy


class Lattice:
    def __init__(self, N, filename=None):
        self.N = N
        #bonds_horizontal[i.j] is between point (i,j) and (i,j+1)
        #bonds_vertical[i.j] is between point (i,j) and (i+1,j)
        self.vortices = np.array([[(np.cos(theta), np.sin(theta)) for theta in np.random.uniform(0, 2*np.pi, N)] for _ in range(N)])
        self.bonds_horizontal = np.full((N, N), True)
        self.bonds_vertical = np.full((N, N), True)

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
    
    def set_vertical_bond(self, x, y, value):
        self.bonds_vertical[x % self.N, y % self.N] = value
        
    def set_all_bonds(self, value):
        self.bonds_horizontal.fill(value)
        self.bonds_vertical.fill(value)

    def set_all_bonds_random(self):
            self.bonds_horizontal = np.random.choice([True, False], (self.N, self.N))
            self.bonds_vertical = np.random.choice([True, False], (self.N, self.N))

    def find_largest_cluster(self):
        visited = np.zeros((self.N, self.N), dtype=bool)
        cluster_sizes = []

        def dfs_iterative(x, y):
            stack = [(x, y)]
            size = 0

            while stack:
                x, y = stack.pop()
                if 0 <= x < self.N and 0 <= y < self.N and not visited[x, y]:
                    visited[x, y] = True
                    size += 1

                    # Check 4 connected bonds
                    if self.get_horizontal_bond(x, y) and not visited[(x + 1) % self.N, y]:
                        stack.append(((x + 1) % self.N, y))
                    if self.get_horizontal_bond(x - 1, y) and not visited[x - 1, y]:
                        stack.append((x - 1, y))
                    if self.get_vertical_bond(x, y) and not visited[x, (y + 1) % self.N]:
                        stack.append((x, (y + 1) % self.N))
                    if self.get_vertical_bond(x, y - 1) and not visited[x, y - 1]:
                        stack.append((x, y - 1))

            return size

        for i in range(self.N):
            for j in range(self.N):
                if not visited[i, j]:
                    cluster_sizes.append(dfs_iterative(i, j))

        return max(cluster_sizes) if cluster_sizes else 0


    def find_all_clusters(self):
        #Find all clusters in the lattice
        visited = np.zeros((self.N, self.N), dtype=bool)
        clusters = []

        def dfs(i, j, current_cluster):
            if visited[i, j]:
                return
            visited[i, j] = True
            current_cluster.append((i, j))

            # Check neighbors
            neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
            for x, y in neighbors:
                # Apply periodic boundary conditions
                x %= self.N
                y %= self.N

                # If the bond exists and the site hasn't been visited, continue the DFS
                if ((x == (i+1)%self.N and self.get_horizontal_bond(i, j)) or
                    (x == (i-1)%self.N and self.get_horizontal_bond(x, y)) or
                    (y == (j+1)%self.N and self.get_vertical_bond(i, j)) or
                    (y == (j-1)%self.N and self.get_vertical_bond(x, y))):
                    dfs(x, y, current_cluster)

        # Iterate over all sites and start a DFS if the site hasn't been visited
        for i in range(self.N):
            for j in range(self.N):
                if not visited[i, j]:
                    current_cluster = []
                    dfs(i, j, current_cluster)
                    if current_cluster:
                        clusters.append(current_cluster)

        return clusters

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
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{self.vortices[i, j][0]} {self.vortices[i, j][1]} ")
                f.write("\n")
            
            # Save horizontal bonds
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{int(self.bonds_horizontal[i, j])} ")
                f.write("\n")
            
            # Save vertical bonds
            for i in range(self.N):
                for j in range(self.N):
                    f.write(f"{int(self.bonds_vertical[i, j])} ")
                f.write("\n")

    def load_from_txt(self, filename):
        with open(filename, 'r') as f:
            lines = f.readlines()
            
            # Load vortices
            for i in range(self.N):
                values = list(map(float, lines[i].split()))
                for j in range(0, len(values), 2):
                    self.vortices[i, j//2] = (values[j], values[j+1])
            
            # Load horizontal bonds
            for i in range(self.N):
                values = list(map(int, lines[self.N + i].split()))
                for j in range(self.N):
                    self.bonds_horizontal[i, j] = bool(values[j])
            
            # Load vertical bonds
            for i in range(self.N):
                values = list(map(int, lines[2*self.N + i].split()))
                for j in range(self.N):
                    self.bonds_vertical[i, j] = bool(values[j])

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