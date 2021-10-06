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
grid_dist_id<3,double,aggregate<double>> * g_one_3d[N_g];

// int c_two = 0;
// grid_dist_id<3,double,aggregate<double,double>> * g_two_3d[N_g];

// int c_three = 0;
// grid_dist_id<3,double,aggregate<double,double,double>> * g_three_3d[N_g];

std::string getNodePathAddress(const int i) {
    return "patches/" + std::to_string(i);  // standardize address
}

void openfpm_init_wrapped()
{
    // todo at check runtime
    openfpm_init(nullptr, nullptr);
    std::cout << "openfpm is initialized ? hoping so ... " << is_openfpm_init() << std::endl;  // debug only
    std::cout << "my rank is " << create_vcluster().getProcessUnitID() << std::endl;
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
    std::cout << "openfpm is initialized ? should NOT be ... " << is_openfpm_init() << std::endl;  // debug only
}

static PyObject* openfpm_finalize_wrapper(PyObject *self, PyObject *args)
{
    // todo constrain # args
    openfpm_finalize_wrapped();
    Py_RETURN_NONE;
}

static PyObject* create_grid_wrapper(PyObject *self, PyObject *args)
{
    // parse input as conduit::Node

    PyObject *node;

    if(!PyArg_ParseTuple(args, "O", &node)) {
        return NULL;
    }

    if(!PyConduit_Node_Check(node)) {
        PyErr_SetString(PyExc_TypeError, "'argument must be a conduit.Node instance");

        return NULL;
    }

    conduit::Node& n = *PyConduit_Node_Get_Node_Ptr(node);

    // parse input values
    // todo more elegant

    npy_int64 dim = n["dim"].value();
    npy_int64 n_prop = n["n props"].value();
    npy_int64 gh = n["gh"].value();
    conduit::int64_array p = n["periodicity"].as_int64_array();
    periodicity<3> bc = {p[0], p[1], p[2]};  // todo more elegant?

    conduit::int64_array size = n["size"].as_int64_array();  // todo should be unsigned
    size_t sz[3] = {
        (size_t) size[0], (size_t) size[1], (size_t) size[2]
    };  // todo more elegant casting

    conduit::float64_array _p1 = n["domain/low"].as_float64_array();
    npy_float64 p1[3] = {_p1[0], _p1[1], _p1[2]};

    conduit::float64_array _p2 = n["domain/high"].as_float64_array();
    npy_float64 p2[3] = {_p2[0], _p2[1], _p2[2]};

    Ghost<3,long int> g(gh);
    Box<3,double> domain(p1, p2);

    g_one_3d[c_one] = new grid_dist_id<3, double, aggregate<double>>(
        sz, domain, g, bc
    );

    // populate conduit node
    if (dim == 3 && n_prop == 1)  // todo debug only
    {
        auto & gdb_ext = g_one_3d[c_one]->getLocalGridsInfo();
        for (int i = 0 ; i < gdb_ext.size(); i++) {
            const float_t tmp0[3] = {
                gdb_ext.get(i).GDbox.getLow(0),
                gdb_ext.get(i).GDbox.getLow(1),
                gdb_ext.get(i).GDbox.getLow(2)
            };
            n.set_path_float32_ptr(getNodePathAddress(i) + "/GDBoxLow", tmp0, 3);

            const float_t tmp1[3] = {gdb_ext.get(i).GDbox.getHigh(0),
                    gdb_ext.get(i).GDbox.getHigh(1),
                    gdb_ext.get(i).GDbox.getHigh(2)};
            n.set_path_float32_ptr(getNodePathAddress(i) + "/GDBoxHigh", tmp1, 3);

            const float_t tmp2[3] = {gdb_ext.get(i).Dbox.getLow(0),
                    gdb_ext.get(i).Dbox.getLow(1),
                    gdb_ext.get(i).Dbox.getLow(2)};
            n.set_path_float32_ptr(getNodePathAddress(i) + "/DBoxLow", tmp2, 3);

            const float_t tmp3[3] = {gdb_ext.get(i).Dbox.getHigh(0),
                    gdb_ext.get(i).Dbox.getHigh(1),
                    gdb_ext.get(i).Dbox.getHigh(2)};
            n.set_path_float32_ptr(getNodePathAddress(i) + "/DBoxHigh", tmp3, 3);

            const float_t tmp4[3] = {gdb_ext.get(i).origin.get(0),
                    gdb_ext.get(i).origin.get(1),
                    gdb_ext.get(i).origin.get(2)};
            n.set_path_float32_ptr(getNodePathAddress(i) + "/origin", tmp4, 3);

            for (int j = 0 ; j < n_prop ; j++) {
                if (j == 0) {
                    // todo node.set_path(
                    //    getNodePathAddress(i) + "data",
                    //    (unsigned char *)g_one_3d[c_one]->get_loc_grid(i).template getPointer<0>());
                }
                // todo run with `ddd` node["wow"] = { 1.0,2.0,3.0,4.0};
        
                /* todo else if (j == 1)
                {
                    node[getNodePathAddress(i) + "data"] = (char*) g_one_3d[c_one]->get_loc_grid(i).template getPointer<1>();
                }
                todo else if (j == 2)
                {
                    node[getNodePathAddress(i) + "data"] = (char*) g_one_3d[c_one]->get_loc_grid(i).template getPointer<2>();
                }*/
            }
        }
        c_one += 1;
    }

    Py_RETURN_NONE;  // return PyConduit_Node_Python_Wrap(&n, 0);
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
