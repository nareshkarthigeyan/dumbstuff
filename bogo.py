import random
arr = [0, 5, 9, -4, 3]
i = 0;
while True:
    i += 1;
    random.shuffle(arr)
    if arr == sorted(arr):
        break;
print(i)