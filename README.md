# XY_model

###Lattice.py

The Lattice class characterize an XY model on square N*N lattice. The `vortices` is an array of 2 tuple $(\sigma_x,\sigma_y)$, the `bonds_horizontal/vertical` is an array of bool values denoting whether the bonds are connected. The save and load via txt file should be improved to npz file if needed. The current `__main__` code visualize the state, with connected bonds in bold line, and disconnected bonds in thin line. Color on vortices denotes angle.


###main.py

This file reflect the BondInsertion section of the note. The current `__main__` code test the program with a randomly initialized lattice. The printed values are the size of largest cluster in each step, and the state of the lattice at t_max. The SpinReset part haven't been implemented yet.
