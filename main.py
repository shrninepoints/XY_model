from Lattice import Lattice
import numpy as np

def hamiltonian(lattice, J = 1):
    H = 0.0
    N = lattice.N

    for i in range(N):
        for j in range(N):
            # Right neighbor
            if lattice.get_horizontal_bond(i, j):
                H -= J * (lattice.vortices[i, j][0] * lattice.vortices[(i + 1) % N, j][0] +
                                lattice.vortices[i, j][1] * lattice.vortices[(i + 1) % N, j][1])

            # Down neighbor
            if lattice.get_vertical_bond(i, j):
                H -= J * (lattice.vortices[i, j][0] * lattice.vortices[i, (j + 1) % N][0] +
                                lattice.vortices[i, j][1] * lattice.vortices[i, (j + 1) % N][1])

    return H




# Example usage:
lattice = Lattice(10)
print(hamiltonian(lattice))
