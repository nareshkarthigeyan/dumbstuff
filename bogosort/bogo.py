import matplotlib.pyplot as plt
import random
import numpy as np

def bogosort(arr):
    sorted_array = sorted(arr)
    it = 0
    while arr != sorted_array:
        random.shuffle(arr)
        it += 1
    return it

# arr = [1, 5, 3, 2]
sort_times = []
for i in range(100000):
    arr = random.sample(range(100), 4)
    num = bogosort(arr)
    sort_times.append(num)
    print(num)

plt.hist(sort_times, bins=30, density=True, alpha=0.6, color='b')

mean, std = np.mean(sort_times), np.std(sort_times)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
pdf = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
plt.plot(x, pdf, 'k', linewidth=2)

# Plot settings
plt.title(f"Distribution of Bogosort Times for Array Size {4}")
plt.xlabel("Iterations")
plt.ylabel("Probability Density")
plt.show()




