from db_config import *

# prev_builds = ['build000','build003','build004','build005','build006','build007']

# for build in session.query(Build).order_by(Build.id):
#     # if build.build_name in prev_builds:
#     #     build.status = 'complete'
#     for plate in build.plates:
#         print(build.build_name,build.status,plate.plate_name)

# session.commit()

for frag in session.query(Fragment).order_by(Fragment.id):
    print(frag.fragment_name)
    if "_link_" in frag.fragment_name:
        print("Linked")
        continue
    elif int(frag.fragment_name[-1]) > 1:
        print("Multi")
        continue
    if frag.seq[8:12] == "GGAG":
        frag.cloning_method = 'entry_1'
    elif frag.seq[8:12] == "AATG":
        frag.cloning_method = 'entry_2'
    else:
        print(frag.seq)
        print(frag.seq[8:12])
        input("Doesn't fit either method")

commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()

#
