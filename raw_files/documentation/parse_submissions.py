import md2py
import os
import sys


# Gets submission.md and converts it to a python object
inFile = sys.argv[1]
with open(inFile,"r") as file:
    markdown = file.read()
    toc = md2py.md2py(markdown)


print("Gene Name:")
print(toc.h1)


