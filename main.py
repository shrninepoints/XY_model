from Lattice import Lattice
import numpy as np
import math
import matplotlib.pyplot as plt

def hamiltonian(lattice, J = 1):
    H = 0.0
    N = lattice.N

    for i in range(N):
        for j in range(N):
            # Nearest neighbors: Right and Down due to periodic boundary conditions
            neighbors = [(i+1, j), (i, j+1)]
            
            for (x, y) in neighbors:
                # Apply periodic boundary conditions
                x %= N; y %= N
                # Dot product of spin vectors
                H -= J * np.dot(lattice.get_vector(i, j), lattice.get_vector(x, y))                
    return H

def calculate_coupling_constants(lattice):
    N = lattice.N
    K_horizontal = np.zeros((N, N))
    K_vertical = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            # For horizontal bond
            # if  s_{i,j,x}s_{i+1,j,x} > 0
            if lattice.get_vector(i, j)[0] * lattice.get_vector((i+1)%N, j)[0] > 0:
                p_ij = np.random.uniform(0, 1)
                K_horizontal[i, j] = -math.log(1 - p_ij) / (2 * lattice.get_vector(i, j)[0] * lattice.get_vector((i+1)%N, j)[0])

            # For vertical bond
            if lattice.get_vector(i, j)[0] * lattice.get_vector(i, (j+1)%N)[0] > 0:
                p_ij = np.random.uniform(0, 1)
                K_vertical[i, j] = -math.log(1 - p_ij) / (2 * lattice.get_vector(i, j)[0] * lattice.get_vector(i, (j+1)%N)[0])

    return K_horizontal, K_vertical

def sort_coupling_constants(K_horizontal, K_vertical, epsilon=1e-9):
    K_values = []
    # Extract values from K_horizontal
    for i in range(K_horizontal.shape[0]):
        for j in range(K_horizontal.shape[1]):
            if abs(K_horizontal[i, j]) > epsilon:  # Discard values close to 0
                K_values.append((K_horizontal[i, j], (i, j), 'horizontal'))

    # Extract values from K_vertical
    for i in range(K_vertical.shape[0]):
        for j in range(K_vertical.shape[1]):
            if abs(K_vertical[i, j]) > epsilon:  # Discard values close to 0
                K_values.append((K_vertical[i, j], (i, j), 'vertical'))

    # Sort the values based on the K value
    sorted_K = sorted(K_values, key=lambda x: x[0])

    return sorted_K

def BondInsert(lattice):
    # a. Clear all bonds and reset parent/size
    lattice.reset_lattice()
    
    # b. Rotate all vortices for a random angle
    lattice.rotate_all_vortices()
    
    # c,d. Calculate and sort K value for all bonds
    K_horizontal, K_vertical = calculate_coupling_constants(lattice)
    sorted_bonds = sort_coupling_constants(K_horizontal, K_vertical)

    # e. recursively insert bonds
    # List to record results
    results = []
    index = 0
    
    # Iterate over the sorted bonds and insert them
    for K_value, (i, j), direction in sorted_bonds:
        # Insert the bond
        if direction == 'horizontal':
            lattice.set_horizontal_bond(i, j, True)
        else:
            lattice.set_vertical_bond(i, j, True)
                             
        # Calculate the size of the largest cluster
        largest_cluster_size = lattice.find_largest_cluster()

        # Record the current state and results
        results.append({
            'index':index,
            'K_value': K_value,
            'bond_position': (i, j),
            'direction': direction,
            'largest_cluster_size': largest_cluster_size,
            'lattice_state': lattice.copy()  # Assuming there's a copy method in the Lattice class
        })
        
        index += 1

    # Find the step with the largest change in cluster size
    largest_change = max(results[1:], key=lambda x: x['largest_cluster_size'] - results[results.index(x) - 1]['largest_cluster_size'])

    return largest_change

def SpinReset(lattice):

    # 1. Save the current state

    # 1. Get all clusters
    clusters = lattice.find_all_clusters()

    # 2. For each cluster, generate a random number and decide whether to flip the x-component
    for cluster in clusters:
        flip = np.random.choice([0, 1])
        if flip:
            for i, j in cluster:
                lattice.vortices[i, j, 0] = -lattice.vortices[i, j, 0]
    return lattice

def Simulation(iterations, system_size):

    # Initialize a random Lattice
    lattice = Lattice(system_size)  # Assuming a size of 10x10 for the lattice, but this can be changed

    # List to save the states of the lattice
    all_results = []

    for _ in range(iterations):
        
        # Apply BondInsertion step
        t_max_state = BondInsert(lattice)

        # Apply SpinReset step
        lattice = SpinReset(t_max_state['lattice_state'])
        
        all_results.append(t_max_state.copy())

    return all_results

def plot(all_results):
    # Extract K values from each t_max_state
    K_values = [state['K_value'] for state in all_results]
    print(K_values)
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(K_values, marker='o', linestyle='-')
    plt.xlabel('Iteration')
    plt.ylabel('K Value')
    plt.title('K Values for Each Iteration')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    iterations = 10  
    system_size = 10
    all_results = Simulation(iterations, system_size)
    print("iterations = ", iterations, "system_size = ", system_size)
    print(all_results[-1])
    all_results[-1]['lattice_state'].save_to_txt("./temp.txt")
    plot(all_results)
    pass