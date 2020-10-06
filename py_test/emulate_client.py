from openfpm import openfpm

openfpm.connect("localhost")

n_rank = openfpm.ranks()
output, error = openfpm.run_command("print(\"Hello world\")")

assert(len(error) == 0)
assert(len(output) == (len("Hello world") + 1)*n_rank)
test_out = "Hello world\n" * n_rank
assert(output == test_out)

g_env = openfpm.grids()
assert(g_env != 0)

for g in g_env:
    assert(g["name"] == "grid_dist_id")
    assert(g["sizes"] == [13,13,13])

g = g_env[0]
sl_z = openfpm.slice(g,z=2)
sl_y = openfpm.slice(g,y=5)
sl_x = openfpm.slice(g,x=8)

#openfpm.continue(steps = 50)

#sl_z = openfpm.slice(z=2)
#sl_y = openfpm.slice(y=5)
#sl_x = openfpm.slice(x=8)

#openfpm.set([2,3,4],[8,8,8],0.0)

#lap_code = "from scipy import ndimage\
#\
#stencil = numpy.array([[0,0,0][0,1,0],[0,0,0]],[[0,1,0],[1,-6,1],[0,1,0]],[[0,0,0],[0, 1, 0],[0,0,0]]])\
#\
#for patch in unit_test_grid.prop(0):\
#patchk = ndimage.convolve(patch, stencil)/(x[1]-x[0])**2"

#output, error = openfpm.run_command(lap_code)

#print(output)
#print(error)
