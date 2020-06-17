# python script that process the examples in particular it pack the examples, it generates the vtk files,
# the images out of the vtk files, the mardown wiki pages

import argparse
import os
from os import listdir
from os.path import isdir, join
from pack_examples import pack_example

parser = argparse.ArgumentParser(description='Pack the examples')
parser.add_argument('directory', help='directory where are located the examples')
args = parser.parse_args()

# List all the folder in the directory

onlydir = [ f for f in listdir(args.directory) if isdir(join(args.directory,f)) ]

p = pack_example()
p.pack(args.directory,args.directory+os.sep+"examples")

