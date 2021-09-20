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


def _bench_cpp(minimon):
    openfpm.openfpm_init()

    grid_node = Node()

    # todo better standardized keys
    grid_node['dim'] = 3
    grid_node['n props'] = 1
    grid_node['gh'] = 4
    grid_node['size'] = 5
    grid_node['periodicity'] = 0  # todo for all dims, MUST BE in {0, 1}
    grid_node['p1'] = 0.0  # todo for all dims
    grid_node['p2'] = 1.0  # todo for all dims

    openfpm.create_grid(grid_node)
    print(grid_node)

    openfpm.openfpm_finalize()


def _bench_py(minimon):
    print('hello world from py!')


def do_benchmark():
    minimon = MiniMon()

    minimon.enter()
    _bench_cpp(minimon)
    minimon.leave('cpp')

    minimon.enter()
    _bench_py(minimon)
    minimon.leave('python')

    minimon.print_stats()


if __name__ == '__main__':
    do_benchmark()
