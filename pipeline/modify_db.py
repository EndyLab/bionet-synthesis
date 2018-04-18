from db_config import *

prev_builds = ['build000','build003','build004','build005','build006','build007']

for build in session.query(Build).order_by(Build.id):
    # if build.build_name in prev_builds:
    #     build.status = 'complete'
    for plate in build.plates:
        print(build.build_name,build.status,plate.plate_name)

# session.commit()



#
