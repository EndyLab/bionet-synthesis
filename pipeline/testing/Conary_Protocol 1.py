"""
Conary Protocol 1
@Robot OT2
@Author Opentrons
"""

from opentrons import labware, instruments, robot

# Labware Set-up
tiprack = labware.load('tiprack-200ul', '1')
tiprack2 = labware.load('tiprack-200ul', '4')

deep_plate = labware.load('96-deep-well', '2')

plate1 = labware.load('96-flat', '3')
plate2 = labware.load('96-flat', '6')

glycerol = labware.load('96-flat', '7')

liquid_trash = labware.load('trash-box', '10')
# Instrument Set-up
pipette = instruments.P300_Multi(mount='left', tip_racks=[tiprack, tiprack2])

plates = [plate1, plate2]
pipette.set_pick_up_current(0.8)
pipette.pick_up_tip()
for plate in plates:
    pipette.transfer(100, glycerol, plate.cols(), new_tip='never')
pipette.drop_tip()

for num, col in enumerate(deep_plate.cols()):
    pipette.pick_up_tip()
    pipette.mix(5, 100, col)
    pipette.aspirate(210, col)
    for plate in plates:
        pipette.dispense(100, plate.cols(num))

    pipette.dispense(10, liquid_trash)
    pipette.drop_tip()
