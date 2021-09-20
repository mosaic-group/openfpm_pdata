import os
from pathlib import Path

import sys
from distutils.core import Extension, setup


def install_openfpm(openfpm_install_dir, openfpm_dep_dir):
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

    openfpm_includes = '-I -Wno-deprecated-declarations   -I.  -I__OFPM__INSTALL__/openfpm_numerics/include -I__OFPM__INSTALL__/openfpm_pdata/include/config -I__OFPM__INSTALL__/openfpm_pdata/include -I__OFPM__INSTALL__/openfpm_data/include -I__OFPM__INSTALL__/openfpm_vcluster/include -I__OFPM__INSTALL__/openfpm_io/include -I__OFPM__INSTALL__/openfpm_devices/include -I__OFPM__DEP__/VCDEVEL/include  -I__OFPM__DEP__/METIS/include -I__OFPM__DEP__/PARMETIS/include -I__OFPM__DEP__/BOOST/include -I__OFPM__DEP__/HDF5/include -I__OFPM__DEP__/LIBHILBERT/include'\
        .replace('__OFPM__INSTALL__', openfpm_install_dir)\
        .replace('__OFPM__DEP__', openfpm_dep_dir)
    compile_args = '-std=c++11 -fPIC -O0'  # -Wall -Wextra
    extra_compile_args = ' '.join([
        conduit_includes,
        openfpm_includes,
        compile_args
    ])

    conduit_libs_folder = conduit_folder / 'lib'
    conduit_libs = ' '.join(map(
        lambda x: '-l' + x,
        ['conduit', 'conduit_blueprint', 'conduit_relay']
    ))
    openfpm_libs_folders = '-L__OFPM__INSTALL__/openfpm_devices/lib -L__OFPM__INSTALL__/openfpm_pdata/lib -L__OFPM__INSTALL__/openfpm_vcluster/lib -L__OFPM__DEP__/VCDEVEL/lib -L__OFPM__DEP__/METIS/lib -L__OFPM__DEP__/PARMETIS/lib -L__OFPM__DEP__/BOOST/lib -L__OFPM__DEP__/HDF5/lib -L__OFPM__DEP__/LIBHILBERT/lib'\
        .replace('__OFPM__INSTALL__', openfpm_install_dir)\
        .replace('__OFPM__DEP__', openfpm_dep_dir)
    openfpm_libs = '-lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc -lrt -ldl -lboost_filesystem -lboost_system -lboost_context -lboost_fiber'
    extra_link_args = ' '.join([
        '-Wl,-Bsymbolic',
        '-L' + str(conduit_libs_folder),
        openfpm_libs_folders,
        conduit_libs,
        openfpm_libs
    ])
    this_folder = os.path.dirname(os.path.realpath(__file__))
    module_source_cpp = this_folder + '/pdata_python.cpp'

    os.environ['CC'] = 'mpicc'  # auto-include mpi
    os.environ['CXX'] = 'mpic++'
    os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':/__OFPM__INSTALL__/openfpm_devices/lib:/__OFPM__INSTALL__/openfpm_vcluster/lib:__OFPM__DEP__/METIS/lib:__OFPM__DEP__/PARMETIS/lib:__OFPM__DEP__/BOOST/lib:__OFPM__DEP__/HDF5/lib:__OFPM__DEP__/LIBHILBERT/lib'
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
