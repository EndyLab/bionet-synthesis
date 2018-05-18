
def option_list(options):
    return options[int(input(''.join([str(counter) + ". " + str(option) + "\n" for counter, option in list(enumerate(options, start=1))]))) - 1]

