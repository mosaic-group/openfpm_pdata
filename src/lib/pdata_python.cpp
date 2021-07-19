// todo does not work with numpy #define PY_SSIZE_T_CLEAN
#include "Grid/grid_dist_id.hpp"
#include <Python.h>

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

void init_openfpm_python()
{
    openfpm_init();
}

void finalize_openfpm_python()
{
    openfpm_finalize();
}

constexpr int N_g = 10;

int c_one = 0;
grid_dist_soa<3,double,aggregate<double>> * g_one_3d[N_g];

int c_two = 0;
grid_dist_soa<3,double,aggregate<double,double>> * g_two_3d[N_g];

int c_three = 0;
grid_dist_soa<3,double,aggregate<double,double,double>> * g_three_3d[N_g];

PyObject * create_grid(long int dim, 
                         long int n_prop,
                         long int gh, 
                         size_t (& sz)[3],
                         long int per[3],
                         double domainP1[3], 
                         double domainP2[3])
{
    // Parsing node, take dimensionality and ghost

    if (dim == 3 && n_prop == 1)
    {
        Ghost<3,long int> g(gh);
        Box<3,double> domain(domainP1,domainP2);
        periodicity<3> bc = {per[0],per[1],per[2]};

        g_one_3d[c_one] = new grid_dist_soa<3,double,aggregate<double>>(sz,domain,g,bc);

        c_one++;

        // Wrap n patch in numpy array

        // 
        conduit::Node node;

        // Populate conduit node

        auto & gdb_ext = g_one_3d.getLocalGrids();

        node["name"] = "grid_3d";

        for (int i = 0 ; i < gdb_ext.size() ; i++)
        {

            node["patches/" + std::to_string(i) + "/GDBoxLow"] = {gdb_ext.get(i).getLow(0),...};
            node["patches/" + std::to_string(i) + "/GDBoxHigh"] = {5,5,5};
            node["patches/" + std::to_string(i) + "/DBoxLow"] = 
            node["patches/" + std::to_string(i) + "/origin"] = {gdb_ext.get(i).get(0),.....};

            for (int j = 0 ; j < n_prop ; j++)
            {
                if (j == 0)
                {
                    node["patches/" +  std::to_string(i)  + "data"] = g_one.get_loc_grid(i).template getPointer<0>();
                }
                else if (j == 1)
                {
                    node["patches/" +  std::to_string(i)  + "data"] = g_one.get_loc_grid(i).template getPointer<1>();
                }
                else if (j == 2)
                {
                    node["patches/" +  std::to_string(i)  + "data"] = g_one.get_loc_grid(i).template getPointer<2>();
                }
            }
        }
    }

    // Convert to python

    

    return NULL;
}

PyObject * create_grid_wrapper(PyObject *self, PyObject *args)
{
    npy_int64 dim, n_prop, gh, sz[3],p[3];
    npy_float64 p1[3],p2[3]

    if (!PyArg_ParseTuple(args, "llllllllldddddd", &dim, &nprop,&gh,&sz[0],&sz[1],sz[2]
                                                                    &p[0],&p[1],&p[2]),
                                                                    &p1[0],&p1[1],&p1[2],
                                                                    &p2[0],&p2[1],&p2[2]); 
    {
        return NULL;
    }

    return create_grid(dim,gh);
}



void delete_grid(long int dim, long int ng)
{
    delete g_one;
}

void delete_grid_wrapper(PyObject *self, PyObject *args)
{
    npy_int64 dim, ng;

    if (!PyArg_ParseTuple(args, "ll", &dim, &ng)) {
        return NULL;
    }

    return delete_grid(dim,ng);
}

static PyMethodDef methods[] = {
    {"create_grid", create_grid_wrapper, METH_VARARGS, ""},
    {"delete_grid", delete_grid_wrapper, METH_VARARGS, ""},
    {"openfpm_init", openfpm_init_python, METH_VARARGS, ""},
    {"openfpm_finalize", openfpm_finalize_python, METH_VARARGS, ""},
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
