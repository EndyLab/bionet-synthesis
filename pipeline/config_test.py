#import sys
#import os

#scriptpath = "./config.py"

# Add the directory containing your module to the Python path (wants absolute paths)
#sys.path.append(os.path.abspath(scriptpath))

# Do the import
#import config

#print(BASE)

config = open("./config.py","r+")
print(config)
print(BASE)
