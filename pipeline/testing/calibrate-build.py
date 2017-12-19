from opentrons import robot, containers, instruments

p10_tiprack = containers.load('tiprack-10ul', 'E2')
p200_tiprack = containers.load('tiprack-200ul', 'A3')

centrifuge_tube = containers.load('tube-rack-2ml','A1')
tubes = containers.load('96-PCR-tall','C2')
strip = containers.load('PCR-strip-tall','C3')

source1 = containers.load('96-flat','D2')
source2 = containers.load('96-flat','D1')
source3 = containers.load('96-flat','D3')
source4 = containers.load('96-flat','B2')


trash = containers.load('point', 'B1', 'holywastedplasticbatman')

p10 = instruments.Pipette(
    axis='a',
    max_volume=10,
    min_volume=0.5,
    tip_racks=[p10_tiprack],
    trash_container=trash,
    channels=8,
    name='p10-8'
)

p200 = instruments.Pipette(
    axis='b',
    max_volume=200,
    min_volume=20,
    tip_racks=[p200_tiprack],
    trash_container=trash,
    channels=1,
    name='p200-1'
)

p200.pick_up_tip()
p200.aspirate(150,centrifuge_tube['A1'].bottom())
p200.dispense(source1.wells('A2'))
p200.drop_tip()


p10.pick_up_tip()

p10.aspirate(9, strip['A1'].bottom())
p10.dispense(tubes.wells('A1').bottom())

p10.aspirate(9, source1.wells('A1').bottom())
p10.dispense(source2.wells('A2').bottom())

p10.aspirate(9, source3.wells('A1').bottom())
p10.dispense(source4.wells('A2').bottom())
p10.drop_tip()







