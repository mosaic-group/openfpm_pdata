import numpy as np
from numba import stencil

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
minimon = MiniMon()

from collections import namedtuple

Point = namedtuple('Point', 'coords')
Box = namedtuple('Box', 'low high')
Domain = namedtuple('Domain', 'U V spacing size')

DIM = 3
PERIODIC = 1
U, V = 0, 0
x, y, z = 0, 1, 2
K, F = 0.053, 0.014  # K and F (Physical constant in the equation)

np.random.seed(42)


def convert_between_ranges(x, domain_x, domain_y):
    ratio = (x - domain_x[0]) / (domain_x[1] - domain_x[0])
    return domain_y[0] + ratio * (domain_y[1] - domain_y[0])


def domain_point2array_index(point, domain, grid_size, to_int=False):
    """ point ~ (1 x) D with coordinates of point in domain space
        domain ~ Domain (with .low and .high)
        grid_size ~ (1 x) D with size of grid in every dimension

        return ~ (1 x) D of type float
    """

    f = lambda d: convert_between_ranges(
        point[d], [domain.low[d], domain.high[d]], [0, grid_size[d]]
    )
    dims = len(point)
    raw = map(f, range(dims))

    if to_int:
        raw = map(np.around, raw)

    return np.float32(list(raw))


def init_domain(old, new, whole_domain):
    old.U = np.ones_like(old.U)
    old.V = np.ones_like(old.V)

    new.U = np.zeros_like(new.U)
    new.V = np.zeros_like(new.V)

    # x_start = old.size[0] * 1.55 / whole_domain.high[0]
    # y_start = old.size[1] * 1.55 / whole_domain.high[1]
    # z_start = old.size[1] * 1.55 / whole_domain.high[2]

    # x_stop = old.size[0] * 1.85 / whole_domain.high[0]
    # y_stop = old.size[1] * 1.85 / whole_domain.high[1]
    # z_stop = old.size[1] * 1.85 / whole_domain.high[2]

    # start = [x_start, y_start, z_start]
    # stop = [x_stop, y_stop, z_stop]

    in_between = 1.2, 2.3  # todo: get from above
    grid_size = [old.U.shape[0], ] * 3  # assuming square grid
    start_i, start_j, start_k = domain_point2array_index(
        [in_between[0], ] * 3, whole_domain, grid_size, to_int=True
    )
    stop_i, stop_j, stop_k = domain_point2array_index(
        [in_between[1], ] * 3, whole_domain, grid_size, to_int=True
    )

    old.U[start_i:stop_i, start_j:stop_j, start_k:stop_k] =\
        0.5 + (np.random.rand(old.U[start_i:stop_i, start_j:stop_j, start_k:stop_k].shape) - 0.5) / 10.0
    old.V[start_i:stop_i, start_j:stop_j, start_k:stop_k] =\
        0.25 + (np.random.rand(old.V[start_i:stop_i, start_j:stop_j, start_k:stop_k].shape) - 0.5) / 20.0


def get_stencil(uFactor, vFactor, dT, format='numpy'):
    stencil = [
        [0, 0, 0],  # Cp
        [0, 0, -1],  # mx
        [0, 0, +1],  # px
        [0, -1, 0],  # my
        [0, +1, 0],  # py
        [-1, 0, 0],  # mz
        [-1, 0, 0]  # pz
    ]

    if format == 'numpy':
        .
    elif format == 'numba':
        def kernel_U(U, V):
            Cp = U[stencil[0]]
            Cpv = V[stencil[0]]

            mx = U[stencil[1]]
            px = U[stencil[2]]

            my = U[stencil[3]]
            py = U[stencil[4]]

            mz = U[stencil[5]]
            pz = U[stencil[6]]

            return Cp +\
                uFactor * (
                    mz + pz + my + py + mx + px - 6 * Cp
                ) -\
                dT * Cp * Cpv * Cpv -\
                dT * F * (Cp - 1.0)  # update based on Eq 2

        def kernel_V(U, V):
            Cp = V[stencil[0]]
            Cpu = U[stencil[0]]

            mx = V[stencil[1]]
            px = V[stencil[2]]

            my = V[stencil[3]]
            py = V[stencil[4]]

            mz = V[stencil[5]]
            pz = V[stencil[6]]

            return Cp +\
                vFactor * (
                    mz + pz + my + py + mx + px - 6 * Cp
                ) -\
                dT * Cpu * Cp * Cp -\
                dT * (F + K) * Cp  # update based on Eq 2




def loop_ofpm(old, new, whole_domain):  # todo compare VS numpy and numba
    pass


def loop_numpy():
    pass


def loop_numba():
    pass


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
    timeSteps = 5000

    # todo calc in .cpp
    old_domain = Domain(U=None, V=None, size=None, spacing=0.0)
    new_domain = Domain(U=None, V=None, size=None, spacing=0.0)
    spacing = np.float32([
        old_domain.spacing[0], old_domain.spacing[1], old_domain.spacing[2]
    ])

    init_domain(old_domain, new_domain)

    count = 0  # sync the ghost
    uFactor = deltaT * du / (spacing[x] * spacing[x])
    vFactor = deltaT * dv / (spacing[x] * spacing[x])

    minimon.enter()
    for i in range(timeSteps):
        if i % 300 == 0:
            print('STEP: {:.0f}'.format(i))

        # todo loop

    minimon.leave('tot_sim')
    minimon.print_stats()

    openfpm.openfpm_finalize()
