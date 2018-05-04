from opentrons import robot, containers, instruments

def touch_tip(self, location=None, radius=1.0):
    """
    Touch the :any:`Pipette` tip to the sides of a well,
    with the intent of removing left-over droplets
    Notes
    -----
    If no `location` is passed, the pipette will touch_tip
    from it's current position.
    Parameters
    ----------
    location : :any:`Placeable` or tuple(:any:`Placeable`, :any:`Vector`)
        The :any:`Placeable` (:any:`Well`) to perform the touch_tip.
        Can also be a tuple with first item :any:`Placeable`,
        second item relative :any:`Vector`
    radius : float
        Radius is a floating point number between 0.0 and 1.0, describing
        the percentage of a well's radius. When radius=1.0,
        :any:`touch_tip()` will move to 100% of the wells radius. When
        radius=0.5, :any:`touch_tip()` will move to 50% of the wells
        radius.
    Returns
    -------
    This instance of :class:`Pipette`.
    Examples
    --------
    ..
    >>> p200 = instruments.Pipette(name='p200', axis='a', max_volume=200)
    >>> p200.aspirate(50, plate[0]) # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    >>> p200.dispense(plate[1]).touch_tip() # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    """

    _description = 'Touching tip'
    self.robot.add_command(_description)

    height_offset = 0

    if helpers.is_number(location):
        height_offset = location
        location = None

    # if no location specified, use the previously
    # associated placeable to get Well dimensions
    if location:
        self.move_to(location, strategy='arc')
    else:
        location = self.previous_placeable

    v_offset = (0, 0, height_offset)

    well_edges = [
        location.from_center(x=radius, y=0, z=1),       # right edge
        location.from_center(x=radius * -1, y=0, z=1),  # left edge
        location.from_center(x=0, y=radius, z=1),       # back edge
        location.from_center(x=0, y=radius * -1, z=1)   # front edge
    ]

    # Apply vertical offset to well edges
    well_edges = map(lambda x: x + v_offset, well_edges)

    [self.move_to((location, e), strategy='direct') for e in well_edges]

    return self

def drop_tip(self, location=None, home_after=True):
    """
    Drop the pipette's current tip
    Notes
    -----
    If no location is passed, the pipette defaults to its `trash_container`
    (see :any:`Pipette`)
    Parameters
    ----------
    location : :any:`Placeable` or tuple(:any:`Placeable`, :any:`Vector`)
        The :any:`Placeable` (:any:`Well`) to perform the drop_tip.
        Can also be a tuple with first item :any:`Placeable`,
        second item relative :any:`Vector`
    Returns
    -------
    This instance of :class:`Pipette`.
    Examples
    --------
    ..
    >>> robot.reset() # doctest: +ELLIPSIS
    <opentrons.robot.robot.Robot object at ...>
    >>> tiprack = containers.load('tiprack-200ul', 'A1')
    >>> trash = containers.load('point', 'A1')
    >>> p200 = instruments.Pipette(axis='a', trash_container=trash)
    >>> p200.pick_up_tip(tiprack[0]) # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    >>> # drops the tip in the trash
    >>> p200.drop_tip() # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    >>> p200.pick_up_tip(tiprack[1]) # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    >>> # drops the tip back at its tip rack
    >>> p200.drop_tip(tiprack[1]) # doctest: +ELLIPSIS
    <opentrons.instruments.pipette.Pipette object at ...>
    """

    if not location and self.trash_container:
        location = self.trash_container

    if isinstance(location, Placeable):
        # give space for the drop-tip mechanism
        location = location.bottom(self._drop_tip_offset)

    _description = "Drop_tip {}".format(
        ('at ' + humanize_location(location) if location else '')
    )
    self.robot.add_command(_description)

    if location:
        self.move_to(location, strategy='arc')

    self.motor.move(self._get_plunger_position('drop_tip'))
    if home_after:
        self.motor.home()

    self.motor.move(self._get_plunger_position('bottom'))

    self.current_volume = 0
    self.current_tip(None)
    return self








#
