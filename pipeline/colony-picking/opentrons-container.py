from opentrons import containers

# Colony plate can index up to 41496 locations
containers.create(
    'Colony_Plate',                    # name of you container
    grid=(312, 532),                    # specify amount of (columns, rows)
    spacing=(0.025, 0.025),               # distances (mm) between each (column, row)
    diameter=0.2,                     # diameter (mm) of each well on the plate
    depth=5)                       # depth (mm) of each well on the plate

custom_plate = containers.load('Colony_Plate', 2)

for well in custom_plate.wells():
    print(well)
