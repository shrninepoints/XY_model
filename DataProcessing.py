import matplotlib.pyplot as plt

# Reading the contents of the file
quantity = 'K_values'
with open(f"./SimulationResult/size_64_time_20230914004654/{quantity}.txt", "r") as file:
    content = file.readlines()

# Convert the content to a list of floats
data = [float(line.strip()) for line in content]

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(data, label=f'Data from {quantity}.txt', color='blue')
plt.title(f'Data from {quantity}.txt')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()
