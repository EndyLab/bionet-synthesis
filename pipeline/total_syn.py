from db_config import *

total = 0
for frag in session.query(Fragment).order_by(Fragment.id):
    print(frag.fragment_name, total)
    total += len(frag.seq)

print("Final: ",total)
