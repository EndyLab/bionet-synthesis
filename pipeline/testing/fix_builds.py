import pandas as pd

def well_addresses():
    '''Generates a list of well address A1-H12'''
    letter = ["A","B","C","D","E","F","G","H"]
    number = ["1","2","3","4","5","6","7","8","9","10","11","12"]
    target_well = []
    temp_well = 0
    for n in number:
        for l in letter:
            temp_well = l + n
            target_well.append(temp_well)
    return target_well

part1 = pd.read_csv('../../builds/build006/build006_2018-02-09 11:53:22-1.csv')
part2 = pd.read_csv('../../builds/build006/build006_2018-02-12 17:31:46-2.csv')
part1 = part1[['Gene','Destination']]
part2 = part2[['Gene']]

well_addresses = well_addresses()

part2['Destination'] = well_addresses[48:]

complete = pd.concat([part1,part2])
complete.to_csv('../../builds/build006_2018-02-12 17:31:46-comp.csv')
