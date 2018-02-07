
# Colony plate can index up to 41496 locations
containers.create(
    'Colony_Plate',                    # name of you container
    grid=(156, 266),                    # specify amount of (columns, rows)
    spacing=(0, 0),               # distances (mm) between each (column, row)
    diameter=0.5,                     # diameter (mm) of each well on the plate
    depth=5)                       # depth (mm) of each well on the plate

