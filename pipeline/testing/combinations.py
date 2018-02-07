import pandas as pd

def find_combinations(list, sum):
    if not list:
        if sum == 0:
            return [[]]
        return []
    return find_combinations(list[1:], sum) + \
        [[list[0]] + tail for tail in
         find_combinations(list[1:], sum - list[0])]

x = 0
result_counter = 0
target = 96
numbers = [5, 53, 19, 2, 17, 2, 23, 90, 75, 94, 15]

combinations = []
reactions = []

while x == 0:
    print("next round")
    result = find_combinations(numbers, target)
    print(result)
    if len(result) > 0:
        for comb in result:
            if len(comb) < 4:
                combinations.append(comb)
                reactions.append(target)
                result_counter += 1
        else:
            target -= 1
    else:
        target -= 1

    if target == 50:
        break

solution = pd.DataFrame({
    "Reactions" : reactions,
    "Combinations" : combinations
})

print(solution)
