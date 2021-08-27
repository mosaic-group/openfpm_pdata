import numpy as np
from conduit import Node
from pathlib import Path
import os

# setup (compile) lib
from distutils.core import Extension, setup
import sys
sys.argv.extend(["build_ext", "-i"])

conduit_folder = Path().home() / 'scratch' / 'conduit' / 'install-debug'
conduit_includes = [
    'include/conduit',
    'python-modules/conduit'
]
conduit_includes = ' '.join(map(
    lambda x: '-I' + str(conduit_folder / x),
    conduit_includes
))
pybind_includes = '-I/home/stefano/Work/MPI/OpenFPM/pybind11-2.6.2.tar/pybind11-2.6.2/include'
openfpm_includes = '-I -Wno-deprecated-declarations   -I.  -I/home/stefano/opt/openfpm/install/openfpm_numerics/include -I/home/stefano/opt/openfpm/install/openfpm_pdata/include/config -I/home/stefano/opt/openfpm/install/openfpm_pdata/include -I/home/stefano/opt/openfpm/install/openfpm_data/include -I/home/stefano/opt/openfpm/install/openfpm_vcluster/include -I/home/stefano/opt/openfpm/install/openfpm_io/include -I/home/stefano/opt/openfpm/install/openfpm_devices/include -I/home/stefano/opt/openfpm/deps/VCDEVEL/include  -I/home/stefano/opt/openfpm/deps/METIS/include -I/home/stefano/opt/openfpm/deps/PARMETIS/include -I/home/stefano/opt/openfpm/deps/BOOST/include -I/home/stefano/opt/openfpm/deps/HDF5/include -I/home/stefano/opt/openfpm/deps/LIBHILBERT/include'
compile_args = '-std=c++11 -fPIC'  # -Wall -Wextra
extra_compile_args = ' '.join([
    conduit_includes,
    pybind_includes,
    openfpm_includes,
    compile_args
])

conduit_libs_folder = conduit_folder / 'lib'
conduit_libs = ' '.join(map(
    lambda x: '-l' + x,
    ['conduit', 'conduit_blueprint', 'conduit_relay']
))
openfpm_libs_folders = '-L/home/stefano/opt/openfpm/install/openfpm_devices/lib -L/home/stefano/opt/openfpm/install/openfpm_pdata/lib  -L/home/stefano/opt/openfpm/install/openfpm_vcluster/lib -L/home/stefano/opt/openfpm/deps/VCDEVEL/lib  -L/home/stefano/opt/openfpm/deps/METIS/lib -L/home/stefano/opt/openfpm/deps/PARMETIS/lib  -L/home/stefano/opt/openfpm/deps/BOOST/lib -L/home/stefano/opt/openfpm/deps/HDF5/lib -L/home/stefano/opt/openfpm/deps/LIBHILBERT/lib'
openfpm_libs = '-lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc -lrt -ldl -lboost_filesystem -lboost_system -lboost_context -lboost_fiber'
extra_link_args = ' '.join([
    '-Wl,-Bsymbolic',
    '-L' + str(conduit_libs_folder),
    openfpm_libs_folders,
    conduit_libs,
    openfpm_libs
])
module_source_cpp = '/home/stefano/opt/openfpm_pdata/src/lib/pdata_python.cpp'

os.environ['CC'] = 'mpicc'  # auto-include mpi
os.environ['CXX'] = 'mpic++'
os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + '://home/stefano/opt/openfpm/install/openfpm_devices/lib://home/stefano/opt/openfpm/install/openfpm_vcluster/lib:/home/stefano/opt/openfpm/deps/METIS/lib:/home/stefano/opt/openfpm/deps/PARMETIS/lib:/home/stefano/opt/openfpm/deps/BOOST/lib:/home/stefano/opt/openfpm/deps/HDF5/lib:/home/stefano/opt/openfpm/deps/LIBHILBERT/lib'
os.environ['PURE_PYTHON'] = '1'

setup(
    ext_modules=[
        Extension(
            'openfpm',
            sources=[
                module_source_cpp
            ],
            extra_compile_args=extra_compile_args.split(),
            extra_link_args=extra_link_args.split(),
            language='c++'
        )
    ]
)
import openfpm

node = Node()
node['inp'] = 8.0
print(node)

openfpm.f(node)
print(node)

openfpm.openfpm_init()

grid_node = openfpm.create_grid(
    3,  # dim
    2,  # # props
    4,  # gh
    5,  # size[0]
    5,  # size[1]
    5,  # size[2]
    7,  # p[0]
    7,  # p[1]
    7,  # p[2]
    3.14,  # p1[0]
    3.14,  # p1[1]
    3.14,  # p1[2]
    40.2,  # p2[0]
    40.2,  # p2[1]
    40.2,  # p2[2]
)
print(grid_node)

openfpm.openfpm_finalize()
