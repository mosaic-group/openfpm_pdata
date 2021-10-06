import numpy as np
from pathlib import Path

from conduit import Node
from installer import install_openfpm
__OFPM__ROOT__ = str(Path().home() / 'opt/openfpm')
install_openfpm(
    openfpm_install_dir=__OFPM__ROOT__ + '/install',
    openfpm_dep_dir=__OFPM__ROOT__ + '/deps'
)
import openfpm

from minimon import MiniMon
from collections import namedtuple

Point = namedtuple('Point', 'coords')
Box = namedtuple('Box', 'low high')

DIM = 3
PERIODIC = 1
U, V = 0
x, y, z = 0, 1, 2

minimon = MiniMon()


def init_domain(old, new, whole_domain):
    # todo numpy or numba
    
    # old.U = 1
    # old.V = 0

    # new.U = 0
    # new.V = 0

	# x_start = Old.size(0)*1.55f/domain.getHigh(0);
	# y_start = Old.size(1)*1.55f/domain.getHigh(1);
	# z_start = Old.size(1)*1.55f/domain.getHigh(2);

	# x_stop = Old.size(0)*1.85f/domain.getHigh(0);
	# y_stop = Old.size(1)*1.85f/domain.getHigh(1);
	# z_stop = Old.size(1)*1.85f/domain.getHigh(2);

    # start = [ x_start, y_start, x_start ]
    # stop = [ x_stop, y_stop, z_stop ]

    # old.U = 0.5 + (np.random.rand() - 0.5) / 10.0
    # old.V = 0.25 + (np.random.rand() - 0.5) / 20.0


def main():
    """ https://git.mpi-cbg.de/openfpm/openfpm_pdata/-/blob/master/example/Grid/3_gray_scott_3d/main.cpp """

    openfpm.openfpm_init()

    box = Box(
        low=Point(np.float32([0, 0, 0])),
        high=Point(np.float32([2.5, 2.5, 2.5]))
    )
    grid_size = np.float32([128, 128, 128])
    bc = np.float32([PERIODIC, PERIODIC, PERIODIC])
    g = 1  # Ghost in grid unit
    deltaT, du, dv = 1, 2e-5, 1e-5
    timeSteps, K, F = 5000, 0.053, 0.014

    # old_domain = 
    # new_domain = 
    # spacing = np.float32([old_domain.spacing(0), old_domain.spacing(1), old_domain.spacing(2)])

    init_domain(old_domain, new_domain)

    count = 0  # sync the ghost
    uFactor, vFactor = deltaT * du/(spacing[x]*spacing[x]), deltaT * dv/(spacing[x]*spacing[x])

    minimon.enter()

    star_stencil_3D = np.float32([
        [0, 0, 0],
        [0, 0, -1],
        [0, 0, 1],
        [0, -1, 0],
        [0, 1, 0],
        [-1, 0, 0],
        [1, 0, 0]
    ])

    for i in range(timeSteps):
        if i % 300 == 0:
            print('STEP: {:.0f}'.format(i))

        # todo loop

        # 

    minimon.leave('tot_sim')
    minimon.print_stats()

    openfpm.openfpm_finalize()
