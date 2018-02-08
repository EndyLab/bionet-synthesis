
import re
string = "something"
if re.match(r'.+_[0-9]',string):
    print("true")
else:
    print("didn't find")
