/**
 * For license see https://git.mpi-cbg.de/openfpm/openfpm_pdata. Especially we don't want GitHub Copilot parsing it.
*/


// todo does not work with numpy #define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "Grid/grid_dist_id.hpp"

#include <numpy/arrayobject.h>

// use  proper strdup
#ifdef CONDUIT_PLATFORM_WINDOWS
    #define _conduit_strdup _strdup
#else
    #define _conduit_strdup strdup
#endif

#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_blueprint.hpp"

// conduit python module capi header
#include "conduit_python.hpp"

constexpr int N_g = 10;

template<unsigned int dim, typename St, typename T, typename Memory = HeapMemory, typename Decomposition = CartDecomposition<dim,St>>
using sgrid_dist_soa = grid_dist_id<dim,St,T,Decomposition,Memory,sgrid_soa<dim,T,Memory>>;

template<unsigned int dim, typename St, typename T, typename Memory = HeapMemory, typename Decomposition = CartDecomposition<dim,St>>
using grid_dist_soa = grid_dist_id<dim,St,T,Decomposition,Memory,grid_soa<dim,T,Memory>>;

int c_one = 0;
grid_dist_soa<3,double,aggregate<double>> * g_one_3d[N_g];

int c_two = 0;
grid_dist_soa<3,double,aggregate<double,double>> * g_two_3d[N_g];

int c_three = 0;
grid_dist_soa<3,double,aggregate<double,double,double>> * g_three_3d[N_g];

std::string getNodePathAddress(const int i) {
    return "patches/" + std::to_string(i);  // standardize address
}

void openfpm_init_wrapped()
{
    // todo at check runtime
    openfpm_init(nullptr, nullptr);
}

static PyObject* openfpm_init_wrapper(PyObject *self, PyObject *args)
{
    // todo constrain # args
    openfpm_init_wrapped();
    Py_RETURN_NONE;
}

void openfpm_finalize_wrapped()
{
    openfpm_finalize();
}

static PyObject* openfpm_finalize_wrapper(PyObject *self, PyObject *args)
{
    // todo constrain # args
    openfpm_finalize_wrapped();
    Py_RETURN_NONE;
}

conduit::Node create_grid(const npy_int64 dim,
                         const npy_int64 n_prop,
                         const npy_int64 gh,
                         size_t (& sz)[3],
                         npy_int64 per[3],
                         double domainP1[3],
                         double domainP2[3])
{

    // Wrap n patch in numpy array
    conduit::Node node;

    // Parsing node, take dimensionality and ghost
    if (dim == 3 && n_prop == 1)
    {
        Ghost<3,long int> g(gh);
        Box<3,double> domain(domainP1, domainP2);
        periodicity<3> bc = {per[0], per[1], per[2]};
        g_one_3d[c_one] = new grid_dist_soa<3, double, aggregate<double>>(
            sz, domain, g, bc
        );

        c_one++;

        // Populate conduit node
        auto & gdb_ext = g_one_3d[c_one]->getLocalGridsInfo();
        node["name"] = "grid_3d";

        for (int i = 0 ; i < gdb_ext.size() ; i++)  // todo assuming 3D grid
        {
            node[getNodePathAddress(i) + "/GDBoxLow"] = {
                gdb_ext.get(i).GDbox.getLow(0),
                gdb_ext.get(i).GDbox.getLow(1),
                gdb_ext.get(i).GDbox.getLow(2)
            };
            node[getNodePathAddress(i) + "/GDBoxHigh"] = {
                gdb_ext.get(i).GDbox.getHigh(0),
                gdb_ext.get(i).GDbox.getHigh(1),
                gdb_ext.get(i).GDbox.getHigh(2)
            };
            node[getNodePathAddress(i) + "/DBoxLow"] = {
                gdb_ext.get(i).Dbox.getLow(0),
                gdb_ext.get(i).Dbox.getLow(1),
                gdb_ext.get(i).Dbox.getLow(2)
            };
            node[getNodePathAddress(i) + "/DBoxHigh"] = {
                gdb_ext.get(i).Dbox.getHigh(0),
                gdb_ext.get(i).Dbox.getHigh(1),
                gdb_ext.get(i).Dbox.getHigh(2)
            };
            node[getNodePathAddress(i) + "/origin"] = {
                gdb_ext.get(i).origin.get(0),
                gdb_ext.get(i).origin.get(1),
                gdb_ext.get(i).origin.get(2)
            };
            
            for (int j = 0 ; j < n_prop ; j++)
            {
                if (j == 0)
                {
                    node.set_path(
                        getNodePathAddress(i) + "data",
                        (unsigned char *) g_one_3d[c_one]->get_loc_grid(i).template getPointer<0>()
                    );
                }
                // todo run with `ddd` node["wow"] = { 1.0,2.0,3.0,4.0};
        
                /*else if (j == 1)
                {
                    node[getNodePathAddress(i) + "data"] = (char*) g_one_3d[c_one]->get_loc_grid(i).template getPointer<1>();
                }
                else if (j == 2)
                {
                    node[getNodePathAddress(i) + "data"] = (char*) g_one_3d[c_one]->get_loc_grid(i).template getPointer<2>();
                }*/
            }
        }
    }

    return node;
}

static PyObject* create_grid_wrapper(PyObject *self, PyObject *args)
{
    npy_int64 dim, n_prop, gh, p[3];
    size_t sz[3];
    npy_float64 p1[3], p2[3];

    if (
        !PyArg_ParseTuple(
            args,
            "llllllllldddddd",
            &dim, &n_prop, &gh,
            &sz[0], &sz[1], sz[2],
            &p[0], &p[1], &p[2],
            &p1[0], &p1[1], &p1[2],
            &p2[0], &p2[1], &p2[2]
        )
    )
    {
        Py_RETURN_NONE;
    }

    auto node = create_grid(dim, n_prop, gh, sz, p, p1, p2);
    return PyConduit_Node_Python_Wrap(&node, 1);  // python owns => true
}

static PyObject* delete_grid(long int dim, long int ng)
{
    delete g_one_3d;
    Py_RETURN_NONE;
}

static PyObject* delete_grid_wrapper(PyObject *self, PyObject *args)
{
    npy_int64 dim, ng;

    if (!PyArg_ParseTuple(args, "ll", &dim, &ng)) {
        Py_RETURN_NONE;
    }

    return delete_grid(dim,ng);
}

static PyObject* f_wrapper(PyObject *self, PyObject *args)
{
    PyObject *node;

    if(!PyArg_ParseTuple(args, "O", &node)) {
        return NULL;
    }

    if(!PyConduit_Node_Check(node)) {
        PyErr_SetString(PyExc_TypeError, "'argument must be a conduit.Node instance");

        return NULL;
    }

    conduit::Node& n = *PyConduit_Node_Get_Node_Ptr(node);

    npy_float64 x = n["inp"].value();
    n.set_path("out", x * 2);

    Py_RETURN_NONE;  // return PyConduit_Node_Python_Wrap(&n, 0);
}

static PyMethodDef methods[] = {
    {"create_grid", create_grid_wrapper, METH_VARARGS, ""},
    {"delete_grid", delete_grid_wrapper, METH_VARARGS, ""},
    {"openfpm_init", openfpm_init_wrapper, METH_VARARGS, ""},
    {"openfpm_finalize", openfpm_finalize_wrapper, METH_VARARGS, ""},
    {"f", f_wrapper, METH_VARARGS, ""},
    {NULL, NULL, METH_VARARGS, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "openfpm",  // name of extension module
    "",
    -1,  // size of per-interpreter state of the module, or -1 if the module keeps state in global variables.
    methods
};

PyMODINIT_FUNC PyInit_openfpm(void) {
    import_array();
    import_conduit();
    return PyModule_Create(&module);
}
