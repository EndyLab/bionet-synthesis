    if args.run:
        robot_name = str(os.environ["ROBOT_DEV"][-5:])
        robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
        if robot_number == "1":
            print("Proceeding with run")
        else:
            sys.exit("Run . robot.sh while in the /opentrons/robots directory to change the robot")

