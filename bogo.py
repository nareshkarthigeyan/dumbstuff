import matplotlib.pyplot as plt
import random

def bogosort(arr):
    sorted_array = sorted(arr)
    it = 0
    while arr != sorted_array:
        random.shuffle(arr)
        it += 1
    return it

maxTrials = 10
maxSize = 10

avgBogoSortTimes = []

for i in range(maxSize):
    sort_times = []
    for j in range(maxTrials):
        arr = random.sample(range(100), i)
        sort_times.append(bogosort(arr))
    avgBogoSortTimes.append(sum(sort_times) / len(sort_times))

plt.plot(range(1, maxSize + 1), avgBogoSortTimes, 
marker='o', linestyle='-', color='b')
plt.title("Bogosort Times for Increasing Array Sizes")
plt.xlabel("Array Size")
plt.ylabel("Iterations")
plt.grid()
plt.show()



