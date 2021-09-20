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


def _bench_cpp(minimon, dim=3, n_props=1, gh=4, size=5, periodicity=0, p1=0.0, p2=1.0):
    """
    - periodicity for all dims, MUST BE in {0, 1}
    - p1 for all dims
    - p2 for all dims
    """

    openfpm.openfpm_init()

    grid_node = Node()

    minimon.enter()
    # todo better standardized keys
    grid_node['dim'] = dim
    grid_node['n props'] = n_props
    grid_node['gh'] = gh
    grid_node['size'] = size
    grid_node['periodicity'] = periodicity
    grid_node['p1'] = p1
    grid_node['p2'] = p2
    minimon.leave('populating Node')

    minimon.enter()
    openfpm.create_grid(grid_node)
    minimon.leave('openfpm.create_grid')

    openfpm.openfpm_finalize()
    print(grid_node)


def _bench_py(minimon, dim=3, n_props=1, gh=4, size=5, periodicity=0, p1=0.0, p2=1.0):
    Point = namedtuple('Point', 'coords')
    Box = namedtuple('Box', 'low high')
    GridInfo = namedtuple('GridInfo', 'GDbox Dbox origin')
    LocalGridsInfo = [
        GridInfo(
            GDbox=Box(
                low=Point([0, ] * dim),
                high=Point([gh * 2 - 1 + size, ] * dim)
            ),
            Dbox=Box(
                low=Point([gh, ] * dim),
                high=Point([gh - 1 + size, ] * dim)
            ),
            origin=Point([-gh, ] * dim)
        )
        for _ in range(n_props)
    ]

    grid_node = Node()
    grid_node['dim'] = dim
    grid_node['n props'] = n_props
    grid_node['gh'] = gh
    grid_node['size'] = size
    grid_node['periodicity'] = periodicity
    grid_node['p1'] = p1
    grid_node['p2'] = p2
    grid_node['patches'] = [
        Node() for _ in range(n_props)
    ]
    for i in range(n_props):
        grid_node['patches']['GDBoxLow'] = LocalGridsInfo[i].GDbox.low
        grid_node['patches']['GDBoxHigh'] = LocalGridsInfo[i].GDbox.high,
        grid_node['patches']['DBoxLow'] = LocalGridsInfo[i].Dbox.low,
        grid_node['patches']['DBoxHigh'] = LocalGridsInfo[i].Dbox.high,
        grid_node['patches']['origin'] = LocalGridsInfo[i].origin

    print(grid_node)


def do_benchmark():
    minimon = MiniMon()

    args = {
        'dim': 3,
        'n_props': 1,
        'gh': 4,
        'size': 5,
        'periodicity': 0,
        'p1': 0.0,
        'p2': 1.0
    }

    minimon.enter()
    _bench_cpp(minimon, **args)
    minimon.leave('cpp')

    minimon.enter()
    _bench_py(minimon, **args)
    minimon.leave('python')

    minimon.print_stats()


if __name__ == '__main__':
    do_benchmark()
