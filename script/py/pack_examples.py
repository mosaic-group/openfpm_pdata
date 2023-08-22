__author__ = 'i-bird'

import os
import shutil
import configparser
from os import listdir
from os.path import isdir, join

class pack_example:
    temporary=os.sep + "tmp" + os.sep + "pack_example"
    base_dir=""

    # Get the files to copy
    def get_files(self,dir):
        config = configparser.ConfigParser()
        config.read(dir+os.sep+"config.cfg")
        files = []
        if "pack" in config:
            if "files" in config["pack"]:
                files = config["pack"]["files"]
                files = str.split(files)
                # for each filr add the prefix
                for i in range(0,len(files)):
                    files[i] = dir+os.sep+files[i]
        return files

    # If the directory is an example
    def is_example(self,dir):
        return os.path.isfile(dir+os.sep+"config.cfg")

    # Get the files to pack
    def copy_example(self,dir):
        print("Processing: " + dir)
        files=self.get_files(dir)

        try:
            # Create the temporary folder
            os.makedirs(self.temporary)
        except FileExistsError:
            pass

        # make the example directory

        if len(files) != 0:
            head, tail = os.path.split(files[0])
            try:
                os.makedirs(self.temporary+os.sep+head[len(self.base_dir):])
            except FileExistsError:
                pass

        for file in files:
            print("copy src:"+file+"    dst:"+self.temporary+os.sep+file[len(self.base_dir):])
            os.makedirs(os.path.dirname(self.temporary+os.sep+file[len(self.base_dir):]), exist_ok=True)
            shutil.copy(file,self.temporary+os.sep+file[len(self.base_dir):])


    # Add examples in a recursive way
    def add_r(self, dir):
        if self.is_example(dir):
            self.copy_example(dir)
            return

        # Create the directory in tmp
        os.makedirs(self.temporary+os.sep+dir[len(self.base_dir):]+os.sep)
        makefile = "Makefile" if os.path.isfile(dir+os.sep+"Makefile") else "CMakeLists.txt"

        # copy Makefile
        print("copy src:"+dir+os.sep+makefile+"    dst:"+self.temporary+os.sep+dir[len(self.base_dir):]+os.sep+makefile)
        shutil.copy(dir+os.sep+makefile,self.temporary+os.sep+dir[len(self.base_dir):]+os.sep+makefile)

        # List all the folder in the directory
        onlydir = [ f for f in listdir(dir) if isdir(join(dir,f)) ]
        for example_dir in onlydir:
            self.add_r(dir+os.sep+example_dir)

    # Add an example path
    def pack(self, dir, file):

        # Remove temporary folders and recreate it
        if os.path.exists(self.temporary):
            shutil.rmtree(self.temporary)
        os.makedirs(self.temporary)

        if self.is_example(dir):
            self.copy_example(dir)

        self.base_dir = dir

        # copy Makefile
        print("copy src:"+dir+os.sep+"Makefile"+"    dst:"+self.temporary+os.sep+"Makefile")
        shutil.copy(dir+os.sep+"Makefile",self.temporary+os.sep+"Makefile")

        # List all the folder in the directory
        onlydir = [ f for f in listdir(dir) if isdir(join(dir,f)) ]
        for example_dir in onlydir:
            self.add_r(dir+os.sep+example_dir)

        # make an archive
        shutil.make_archive(file,"bztar",self.temporary)
