import random

with open("input.txt", 'a') as f:
    for _ in range(10000):
        number = random.randint(-10000, 10000)
        f.write(f"{number}\n")
        f.flush()
        print(f"writing {number} to file")