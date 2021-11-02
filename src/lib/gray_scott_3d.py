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
from dataclasses import dataclass
from copy import deepcopy

Point = namedtuple('Point', 'coords')
Box = namedtuple('Box', 'low high')


@dataclass
class Domain:
    U: np.ndarray
    V: np.ndarray
    spacing: np.ndarray
    size: np.ndarray


DIM = 3
PERIODIC = 1
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
        point[d], [domain.low.coords[d], domain.high.coords[d]], [0, grid_size[d]]
    )
    dims = len(point)
    raw = map(f, range(dims))

    if to_int:
        raw = map(np.around, raw)
        return np.int64(list(raw))

    return np.float64(list(raw))


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
        0.5 + (np.random.rand(*old.U[start_i:stop_i, start_j:stop_j, start_k:stop_k].shape) - 0.5) / 10.0
    old.V[start_i:stop_i, start_j:stop_j, start_k:stop_k] =\
        0.25 + (np.random.rand(*old.V[start_i:stop_i, start_j:stop_j, start_k:stop_k].shape) - 0.5) / 20.0


def get_stencils(uFactor, vFactor, dT, format='numpy'):
    if format == 'numpy':
        return None, None
    elif format == 'numba':  # todo
        @stencil
        def kernel_U(U, V, dT):  # todo without dT
            Cp = U[0, 0, 0]  # todo avoid repeating locations
            Cpv = V[0, 0, 0]

            mx = U[0, 0, -1]
            px = U[0, 0, +1]

            my = U[0, -1, 0]
            py = U[0, +1, 0]

            mz = U[-1, 0, 0]
            pz = U[-1, 0, 0]

            return Cp +\
                uFactor * (
                    mz + pz + my + py + mx + px - 6 * Cp
                ) -\
                dT * Cp * Cpv * Cpv -\
                dT * F * (Cp - 1.0)  # update based on Eq 2

        @stencil
        def kernel_V(U, V, dT):
            Cp = V[0, 0, 0]  # todo avoid repeating locations
            Cpu = U[0, 0, 0]

            mx = V[0, 0, -1]
            px = V[0, 0, +1]

            my = V[0, -1, 0]
            py = V[0, +1, 0]

            mz = V[-1, 0, 0]
            pz = V[-1, 0, 0]

            return Cp +\
                vFactor * (
                    mz + pz + my + py + mx + px - 6 * Cp
                ) -\
                dT * Cpu * Cp * Cp -\
                dT * (F + K) * Cp  # update based on Eq 2

        return kernel_U, kernel_V


def make_loop(timeSteps, debug_every=300, save_every=500):
    should_debug = lambda timeStep: timeStep % debug_every == 0
    should_save = lambda timeStep: timeStep % save_every == 0

    def _f(stencils, old_copy, new_copy, dT):
        kernel_U, kernel_V = stencils  # unpack

        for i in range(timeSteps):
            if should_debug(i):
                print('STEP: {:d}'.format(i))

            old_U, old_V = np.zeros_like(old_copy.U), np.zeros_like(old_copy.V)

            new_copy.U = kernel_U(old_U, old_V, dT)
            new_copy.V = kernel_V(old_U, old_V, dT)

            old_copy.U = new_copy.U  # sync
            old_copy.V = new_copy.V

            if should_save(i):
                # todo save
                print('    saved to todo')

    return _f


def create_grid_dist(size, domain, gh, periodicity=[0, ] * 3):
    grid_node = Node()
    grid_node['dim'] = DIM
    grid_node['n props'] = 2
    grid_node['gh'] = gh
    grid_node['size'] = size
    grid_node['periodicity'] = periodicity
    grid_node['domain/low'] = domain.low.coords
    grid_node['domain/high'] = domain.high.coords

    minimon.enter()
    openfpm.create_grid(grid_node)
    minimon.leave('openfpm.create_grid')

    return Domain(
        U=np.zeros(size),
        V=np.zeros(size),
        spacing=np.float64(grid_node['spacing']),
        size=size
    )


def main():
    """ https://git.mpi-cbg.de/openfpm/openfpm_pdata/-/blob/master/example/Grid/3_gray_scott_3d/main.cpp """

    whole_domain = Box(
        low=Point(np.float64([0, 0, 0])),
        high=Point(np.float64([2.5, 2.5, 2.5]))
    )
    grid_size = np.int64([128, ] * 3)
    periodicity = np.int64([PERIODIC, ] * 3)
    g = 1  # Ghost in grid unit
    deltaT, du, dv = 1, 2e-5, 1e-5
    timeSteps = 2  # todo debug only 5000

    old_copy = create_grid_dist(
        grid_size,
        whole_domain,
        g,
        periodicity=periodicity
    )
    new_copy = deepcopy(old_copy)
    spacing = deepcopy(old_copy.spacing)

    init_domain(old_copy, new_copy, whole_domain)

    uFactor = deltaT * du / (spacing[x] * spacing[x])
    vFactor = deltaT * dv / (spacing[x] * spacing[x])
    stencils = get_stencils(uFactor, vFactor, deltaT, format='numba')
    loop = make_loop(timeSteps, debug_every=300, save_every=500)

    minimon.enter()
    loop(stencils, old_copy, new_copy, deltaT)
    minimon.leave('tot_sim')
    minimon.print_stats()


if __name__ == '__main__':
    openfpm.openfpm_init()

    main()

    openfpm.openfpm_finalize()
