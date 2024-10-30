import random

datasize = 1000

with open("input.txt", 'w') as f:
    for _ in range(datasize):
        number = random.randint(-datasize, datasize)
        f.write(f"{number}\n")
        f.flush()
        print(f"writing {number} to file")

with open("reverse.txt", 'w') as f:
    for i in range(datasize, -1, -1):
        f.write(f"{i}\n")
        f.flush()

with open("sorted.txt", 'w') as f:
    for i in range(datasize):
        f.write(f"{i}\n")
        f.flush()
