import argparse
import json
import os
import sys
import shutil
from pathlib import Path

# export PATH="/home/i-bird/Desktop/MOSAIC/OpenFPM_project/HIP/hipify-clang/dist/:$PATH"
# python script/py/code-transform.py --json "build/compile_commands.json" --compile_command "/usr/local/cuda/bin/nvcc -ccbin=/home/i-bird/MPI/bin/mpic++" --build_path "/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_bianucci_flush_hip_test/openfpm_pdata/build" --add_includes " -I/home/i-bird/MPI/include -I/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_bianucci_flush_hip_test/openfpm_pdata/openfpm_numerics/src/" --directory "/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_bianucci_flush_hip_test/openfpm_pdata" --exclude_directories="/home/i-bird/Desktop/MOSAIC/OpenFPM_project/openfpm_bianucci_flush_hip_test/openfpm_pdata/openfpm_data/src/util/cuda/cub"

cudaExtension = ".cu"
cudaFileToParse=[]
commands=[]
compile_command_to_substitute=""
build_path=""

# function to get unique values 
def unique(list1,list2): 
  
    # intilize a null list 
    unique_list1 = []
    unique_list2 = []
      
    # traverse for all elements 
    for i in range(len(list1)): 
        # check if exists in unique_list or not 
        if list1[i] not in unique_list1:
            unique_list1.append(list1[i])
            unique_list2.append(list2[i])

    return unique_list1, unique_list2



parser = argparse.ArgumentParser(description='Transform the code to HIP or Cupla')
parser.add_argument('--directory',help='convert all files in the directory recursively')
parser.add_argument('--exclude_directories',help='exclude directories from hipify')
parser.add_argument('--extensions',help='extension used to reconize Cuda file examples: --extensions "cu cuh"')
parser.add_argument('--json', help='json compile file')
parser.add_argument('--compile_command', help='The command used to compile (it is going to be substituted with hipify-clang make sure is in your PATH)')
parser.add_argument('--build_path',help='The path where is the compile_commands.json')
parser.add_argument('--add_includes',help='Add additional includes directories (or compilation options) examples --add_includes= " -I/include/dir1 -I/include/dir2" Warning the space at the beginning is fundamental')
parser.add_argument('--remove-hipyfied',help='Remove hipyfied files')
parser.add_argument('--restore',help='Restore original files',action='store_true')


args = parser.parse_args()

print('Starting')
print('compile_command: ' + str(args.compile_command))
print('add_includes: ' + str(args.add_includes))
print('exclude_directories: ' + str(args.exclude_directories))
print("Restore: " + str(args.restore))

if args.compile_command == None:
    sys.exit('Error! you must pass the option --compile_command')

if args.build_path == None:
    sys.exit('Error! you must pass the option --build_path')

if args.json == None:
    sys.exit('Error! you must give the option --json')

# Parse json

if args.json is not None:

    data = json.load(open(args.json))
    for p in data:
        if str(p['file']).endswith(cudaExtension):
            cudaFileToParse.append(p['file'])
            commands.append(p['command'])

# explore directory

files_to_convert = []

if args.directory is not None:

    print("Check directory: " + args.directory)

    arr_exclude = []
    if args.exclude_directories is not None:
        arr_exclude = args.exclude_directories.split()

    for filename in Path(args.directory).rglob('*.cuh'):
        if len(arr_exclude) is not 0:
            for exclude in arr_exclude:
                if not str(filename).startswith(exclude):
                    files_to_convert.append(filename)
        else:
            files_to_convert.append(filename)

# construct a compile command for header files

header_command = set([])

for command in commands:
    command = command.replace(args.compile_command," ",1)
    options = command.split()
    options = [option for option in options if option.startswith("-I") ]
    header_command.update(options)


for i in range(len(files_to_convert)):
    print('[' + str(i+1) + '/' + str(len(files_to_convert)) + ']   ' + 'Converting CUDA File:' + str(files_to_convert[i]))

    #Safe the original

    if not os.path.exists(str(files_to_convert[i]) + "_original"):
        shutil.copy(files_to_convert[i],str(files_to_convert[i]) + "_original")
    else:
        print("SKIPPING Original saving")

    #if os.path.exists(os.path.splitext(files_to_convert[i])[0] + "original.cuh"):
    #    os.remove(os.path.splitext(files_to_convert[i])[0] + "original.cuh")

    if args.restore is True:
        if os.path.exists(str(files_to_convert[i]) + "_original"):
            shutil.copy(str(files_to_convert[i]) + "_original",files_to_convert[i])
        continue

    hip_com = "hipify-clang -D__HIPIFY__  -D__NVCC__ -D__CUDACC_VER_MAJOR__=10 -D__CUDACC_VER_MINOR__=1 -D__CUDACC_VER_MINOR__ -D__CUDACC_VER_BUILD__=243 " + str(files_to_convert[i])
    hip_com += " -o " + os.path.splitext(files_to_convert[i])[0] + '_hipified.hpp'
    hip_com += " -p " + args.build_path
    if args.add_includes != None :
        hip_com += "  " + args.add_includes 
    header_command_str = " ".join(str(hc) for hc in header_command)
    hip_com += "  -- " + header_command_str
    print('Launching ' + hip_com)
    stream = os.popen(hip_com)
    output = stream.read()
    print('----------------------------------------------------------------------------')

# Parse Json file


[cudaFileToParse,commands] = unique(cudaFileToParse,commands)


for i in range(len(cudaFileToParse)):
    print('[' + str(i+1) + '/' + str(len(cudaFileToParse)) + ']   ' + 'Converting CUDA File:' + cudaFileToParse[i])

    #Safe the original
    if not os.path.exists(str(cudaFileToParse[i]) + "_original"):
        shutil.copy(cudaFileToParse[i],str(cudaFileToParse[i]) + "_original")
    else:
        print("SKIPPING Original saving")

    if args.restore is True:
        if os.path.exists(str(cudaFileToParse[i]) + "_original"):
            shutil.copy(str(cudaFileToParse[i]) + "_original",cudaFileToParse[i])
        continue


    hip_com = "hipify-clang -D__HIPIFY__  -D__NVCC__ -D__CUDACC_VER_MAJOR__=10 -D__CUDACC_VER_MINOR__=1 -D__CUDACC_VER_MINOR__ -D__CUDACC_VER_BUILD__=243 " + cudaFileToParse[i]
    hip_com += " -o " + os.path.splitext(cudaFileToParse[i])[0] + '_hipified.cpp'
    hip_com += " -p " + args.build_path
    if args.add_includes != None :
        hip_com += "  " + args.add_includes
    hip_com += "  -- "
    commands[i] = commands[i].replace(args.compile_command,hip_com,1)
    print('Launching ' + commands[i])
    stream = os.popen(commands[i])
    output = stream.read()
    print('----------------------------------------------------------------------------')

if args.restore is True:
    sys.exit(0)

###############3 Replace

for i in range(len(files_to_convert)):
    if os.path.exists(os.path.splitext(files_to_convert[i])[0] + '_hipified.hpp'):
        shutil.move(os.path.splitext(files_to_convert[i])[0] + '_hipified.hpp',files_to_convert[i])
    print('----------------------------------------------------------------------------')

# Parse Json file


[cudaFileToParse,commands] = unique(cudaFileToParse,commands)


for i in range(len(cudaFileToParse)):
    if os.path.exists(os.path.splitext(cudaFileToParse[i])[0] + '_hipified.cpp'):
        shutil.move(os.path.splitext(cudaFileToParse[i])[0] + '_hipified.cpp',cudaFileToParse[i])
    print('----------------------------------------------------------------------------')



