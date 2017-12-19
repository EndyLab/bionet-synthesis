from opentrons import robot, containers, instruments

#p10_tiprack_single = containers.load('tiprack-10ul', 'E3')

p10_tiprack_single = containers.load('tiprack-10ul', 'E1')


p200_tiprack = containers.load('tiprack-200ul', 'A3')

tubes_single = containers.load('96-PCR-tall','C2')
strip_single = containers.load('PCR-strip-tall','C3')

source1_single = containers.load('96-flat','D2')
source2_single = containers.load('96-flat','D1')
source3_single = containers.load('96-flat','B2')
source4_single = containers.load('96-flat','D3')


trash = containers.load('point', 'B1', 'holywastedplasticbatman')

p10single = instruments.Pipette(
    axis='a',
    max_volume=10,
    min_volume=0.5,
    tip_racks=[p10_tiprack_single],
    trash_container=trash,
    channels=1,
    name='p10-8s'
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

p10single.pick_up_tip()

p10single.aspirate(9, strip_single['A1'].bottom())
p10single.dispense(tubes_single.wells('A1').bottom())

p10single.aspirate(9, source1_single.wells('A1').bottom())
p10single.dispense(source2_single.wells('A2').bottom())

p10single.aspirate(9, source3_single.wells('A1').bottom())
p10single.dispense(source4_single.wells('A2').bottom())



p10single.drop_tip()







