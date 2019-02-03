import os
import sys

user_dir = sys.argv[1]

with open('monomer_posfiles.txt', 'w') as derp:
    for file in os.listdir(os.path.join(os.getcwd(), user_dir)):
        if '_0.pos' in file:
            derp.write(file)
            derp.write('\n')