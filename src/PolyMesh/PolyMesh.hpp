#ifndef __POLYMESH_HPP__
#define __POLYMESH_HPP__

#include <Vector/vector_dist.hpp>
#include "voro++.hh"
#include "util/polymesh_base_functions.hpp"
#include "Space/Shape/Plane.hpp"
#include "hash_map/hopscotch_map.h"
#include "hash_map/hopscotch_set.h"
#include "PolyMesh_gradients.hpp"

template<unsigned int dim, typename T, typename TypeVolume, typename TypeSurface, typename TypeEdge, typename TypeVertex>
class PolyMesh
{
    // Transform volumes, we add connectivity
	// Connectivity (int,int)
	// The first integer is the number of faces
	// The second integer indicate where in Surfaces the information are stored
    typedef typename AggregateAppend<int,TypeVolume>::type TypeVolume_p;
    typedef typename AggregateAppend<int,TypeVolume_p>::type TypeVolume_pp;
    typedef typename AggregateAppend<int,TypeVolume_pp>::type TypeVolume_ppp;
    typedef typename AggregateAppend<int,TypeVolume_ppp>::type TypeVolume_pppp;
    typedef typename AggregateAppend<int,TypeVolume_pppp>::type TypeVolume_ppppp;

    // Transform surfaces
    typedef typename AggregateAppend<int,TypeSurface>::type TypeSurface_p;
    typedef typename AggregateAppend<int,TypeSurface_p>::type TypeSurface_pp;
    typedef typename AggregateAppend<int,TypeSurface_pp>::type TypeSurface_ppp;
    typedef typename AggregateAppend<int,TypeSurface_ppp>::type TypeSurface_pppp;
    typedef typename AggregateAppend<int[2],TypeSurface_pppp>::type TypeSurface_ppppp;

    // Transform edges
    typedef typename AggregateAppend<int,TypeEdge>::type TypeEdge_p;
    typedef typename AggregateAppend<int,TypeEdge_p>::type TypeEdge_pp;
    typedef typename AggregateAppend<int,TypeEdge_pp>::type TypeEdge_ppp;
    typedef typename AggregateAppend<int,TypeEdge_ppp>::type TypeEdge_pppp;

    // Transform Vertex
    typedef typename AggregateAppend<int,TypeVertex>::type TypeVertex_p;
    typedef typename AggregateAppend<int,TypeVertex_p>::type TypeVertex_pp;
    typedef typename AggregateAppend<int[4],TypeVertex_pp>::type TypeVertex_ppp;

    static constexpr int LowConIndex = 0;
    static constexpr int PositionFix = 1;
    static constexpr int StartPosFix = 2;

    // We use the volume_surfaces_connectivity as an example, in the case we operate on the
    // surfaces array. If we re-arrange surfaces (make them unique or reorder the indexes), the
    // connectivity become invalid. The way to reconstruct is before reorder, take the old position of
    // ,
    // take start position of 
    // 
    // int: index in the lower connectivity
    // T[3]: Position fix: Connectivity is searched by position when we use the relink.
    //       The position fix is filled in FillOldPositionConnectivity and used in relink
    // int: Start pos fix: The start pos stored by the surfaces is still valid, we have to infact
    //      consider that the edge connectivity remain fixed ergo the start pos in the surface is valid
    //      (StartPos in surfaces contain the position in the edge connectivity of the surface)
    openfpm::vector<aggregate<int,T[3],int>> volume_surfaces_connectivity;
    openfpm::vector<aggregate<int,T[3],int>> surface_edges_connectivity;
    openfpm::vector<aggregate<int,T[3],int>> edge_vertices_connectivity;


    // Delauney tetraedrons
    // 4 indices of the tetraedron
    // 1 index of the vetex
    openfpm::vector<aggregate<int[4],int>> Tet_del;

    //! List of cells
    vector_dist<dim,T,TypeVolume_ppppp> Volumes;

    //! List of cells surfaces
    vector_dist<dim,T,TypeSurface_ppppp> Surfaces;

    //! List of edges
    vector_dist<dim,T,TypeEdge_pppp> Edges;

    //! List of vertices
    vector_dist<dim,T,TypeVertex_ppp> Vertices;

    //! Tollerance (in distance) for which two edges/vertices are considered the same
    T same_element_tollerance = 1e-7;

    bool tet_found(Point<dim,T> & xv,int p, int q, int s)
    {
        Point<dim,T> xp = Volumes.getPos(p);
        Point<dim,T> xq = Volumes.getPos(q);
        Point<dim,T> xs = Volumes.getPos(s);

        xp -= xv;
        xq -= xv;
        xs -= xv;

        double tot = (xp[1]*xq[2] - xp[2]*xq[1])*xs[0] + (xp[0]*xq[2] - xp[2]*xq[0])*xs[1] + (xp[0]*xq[1] - xp[1]*xq[0])*xs[2];

        return fabs(tot) > same_element_tollerance;
    }

    /*! \brief Calculate the delaunay tetraedrons
     *
     * \tparam TetDelStart Tetraedra of the delauny start point
     * \tparam TetDelNum number of the delauny tetraedra for the volume V
     * 
     */
    template <unsigned int TetDelStart, unsigned int TetDelNum, typename vector_deltet_type>
    void calc_delaunay_tetraedron(vector_deltet_type & Tet_del)
    {
        double max_radius = 0.0;

        auto it = Volumes.getDomainIterator();

        while (it.isNext())
        {
            auto p = it.get().getKey();

            Point<dim,T> xp = Volumes.getPos(p);

            Volumes.template getProp<TetDelStart>(p) = Tet_del.size();


            ForAllVolumeVertices(*this,p,[&](openfpm::vector<int> & indices){

                for (int i = 0 ; i < indices.size() ; i++)
                {
                    auto v = indices.get(i);

                    Point<dim,T> xv = Vertices.getPos(v);

                    auto radius = xp.distance(xv);

                    if (radius > max_radius)
                    {max_radius = radius;}
                }
            });

            Volumes.template getProp<TetDelNum>(p) = Tet_del.size() - Volumes.template getProp<TetDelStart>(p);

            ++it;
        }

        Tet_del.clear();

        Volumes.template ghost_get<>();
        auto cl_v = Volumes.getCellList(max_radius*1.1);

        auto it2 = Volumes.getDomainIterator();

        while (it2.isNext())
        {
            auto p = it2.get().getKey();

            Point<dim,T> xp = Volumes.getPos(p);

            Volumes.template getProp<TetDelStart>(p) = Tet_del.size();


            ForAllVolumeVertices(*this,p,[&](openfpm::vector<int> & indices){

                for (int i = 0 ; i < indices.size() ; i++)
                {
                    int id_tot = 0;
                    int ids[16];

                    auto v = indices.get(i);

                    Point<dim,T> xv = Vertices.getPos(v);

                    auto radius = xp.distance(xv);

                    auto it_nn = cl_v.getNNIterator(cl_v.getCell(xv));

                    while (it_nn.isNext())
                    {
                        auto q = it_nn.get();

                        Point<dim,T> xq = Volumes.getPos(q);

                        if (fabs(radius - xq.distance(xv)) < same_element_tollerance)
                        {
                            // found element tetra
                            if (id_tot < 16)
                            {ids[id_tot] = q;}
                            id_tot++;
                        }

                        ++it_nn;
                    }

                    // We should have 4 elements.

                    if (id_tot >= 3)
                    {
                        Tet_del.add();
                        Tet_del.last().template get<DelTelIndices>()[0] = p;
                        Tet_del.last().template get<DelTelIndices>()[1] = ids[0];
                        Tet_del.last().template get<DelTelIndices>()[2] = ids[1];

                        // The last point must not be coplanar

                        for (int j = 2 ; j < id_tot ; j++)
                        {
                            Tet_del.last().template get<DelTelIndices>()[3] = ids[j];

                            if (tet_found(xp,ids[0],ids[1],ids[j]) == true)
                            {break;}
                        }

                        Tet_del.last().template get<DelTelVertId>() = v;
                    }
                    else
                    {
                        std::cout << __FILE__ << ":" << __LINE__ << " was not able to construct delaunay" << std::endl;
                    }
                }
            });

            Volumes.template getProp<TetDelNum>(p) = Tet_del.size() - Volumes.template getProp<TetDelStart>(p);

            ++it2;
        }
    }

    template<unsigned int NeleDwProp, unsigned int StartProp, unsigned int unique_id, typename Ele_type, typename Ele_type_up, typename Connectivity_type>
    void FixIndexes(Ele_type_up & Elements_up,
                    Ele_type & Elements,
                    Connectivity_type & conn)
    {
        auto it_e2 = Elements_up.getDomainIterator();

        while (it_e2.isNext())
        {
            auto pk = it_e2.get();
            auto p = pk.getKey();

            auto start_conn = Elements_up.template getProp<StartProp>(pk);
            auto n_ele_dw = Elements_up.template getProp<NeleDwProp>(pk);

            for (int k = 0 ; k < n_ele_dw ; k++)
            {
                auto ele_dw = conn.template get<0>(start_conn + k);
                auto ele_dw_unique = Elements.template getProp<unique_id>(ele_dw);
                conn.template get<0>(start_conn + k) = ele_dw_unique;
            }

            ++it_e2;
        }
    }

    /* \brief 
     *
     * 
     */
    template<unsigned int prp, unsigned int GID, typename Ele_type, typename cl_type>
    void unify_elements(Ele_type & Elements, cl_type & cl_e,
                        typename std::remove_reference<decltype(Elements.getPosVector())>::type & positions_unique,
                        typename std::remove_reference<decltype(Elements.getPropVector())>::type & properties_unique)
    {
        auto it = Elements.getDomainAndGhostIterator();

        typedef typename boost::mpl::at<typename Ele_type::value_type::type, boost::mpl::int_<prp>>::type index_type;
        constexpr int max_id = std::numeric_limits<index_type>::max();

        while (it.isNext())
        {
            auto p = it.get();

            Elements.template getProp<prp>(p) = max_id;

            ++it;
        }

        positions_unique.clear();
        properties_unique.clear();
        auto it_e = Elements.getDomainIterator();

        while (it_e.isNext())
        {
            auto p = it_e.get().getKey();
            Point<dim,T> xp = Elements.getPos(p);

            auto it_nn = cl_e.getNNIterator(cl_e.getCell(xp));

            while (it_nn.isNext())
            {
                auto q = it_nn.get();

                if (p == q)   {++it_nn; continue;}

                Point<dim,T> xq = Elements.getPos(q);

                if (xp.distance(xq) < same_element_tollerance)
                {
                    auto p_gid = Elements.template getProp<GID>(p);
                    auto q_gid = Elements.template getProp<GID>(q);

                    auto smallest = (p_gid < q_gid)?p_gid:q_gid;
                    smallest = (smallest < Elements.template getProp<prp>(q))?smallest:Elements.template getProp<prp>(q);
                    smallest = (smallest < Elements.template getProp<prp>(p))?smallest:Elements.template getProp<prp>(p);

                    Elements.template getProp<prp>(q) = smallest;
                    Elements.template getProp<prp>(p) = smallest;


                }

                ++it_nn;
            }

            ++it_e;
        }

        // Create unique

        auto it_e2 = Elements.getDomainIterator();

        while (it_e2.isNext())
        {
            auto pk = it_e2.get();
            auto p = pk.getKey();

            if (Elements.template getProp<GID>(p) == Elements.template getProp<prp>(p) || 
                Elements.template getProp<prp>(p) == max_id)
            {
                positions_unique.add();
                properties_unique.add();
                int last = positions_unique.size()-1;

                positions_unique.get(last) = Elements.getPosVector().get(p);
                properties_unique.get(last) = Elements.getPropVector().get(p);

                Elements.template getProp<prp>(p) = last;
            }
            else
            {
                int link = Elements.template getProp<prp>(p);
                int last = Elements.template getProp<prp>(link);

                Elements.template getProp<prp>(p) = last;
            }

            ++it_e2;
        }
    }


    template<typename Ele_type, typename cl_type>
    bool check_uniqueness(Ele_type & Elements, cl_type & cl_e)
    {
        auto it = Elements.getDomainIterator();

        auto it_e = Elements.getDomainIterator();

        while (it_e.isNext())
        {
            int n_unique = 0;

            auto p = it_e.get().getKey();
            Point<dim,T> xp = Elements.getPos(p);

            auto it_nn = cl_e.getNNIterator(cl_e.getCell(xp));

            while (it_nn.isNext())
            {
                auto q = it_nn.get();

                if (p == q)   {++it_nn; continue;}

                Point<dim,T> xq = Elements.getPos(q);

                if (xp.distance(xq) < same_element_tollerance)
                {
                    n_unique++;
                }

                ++it_nn;
            }

            if (n_unique != 0)
            {
                std::cout << __FILE__ << ":" << __LINE__ << " error elements are not unique" << std::endl;
                return false;
            }

            ++it_e;
        }

        return true;
    }

    template<unsigned int StartProp, typename Element_conn_type, typename Element_type>
    void FillOldPositionConnectivity_element(Element_conn_type & element_conn, Element_type & elements)
    {
        for (int i = 0 ; i < element_conn.size() ; i++)
        {
            auto e = element_conn.template get<0>(i);

            for (int j = 0 ; j < dim ; j++)
            {
                element_conn.template get<1>(i)[j] = elements.getPos(e)[j];
                element_conn.template get<2>(i) = elements.template getProp<StartProp>(e);
            }
        }
    }


    template<unsigned int StartProp,typename Elements_conn_type, typename Elements_type, typename CellList_type>
    void relink_elements(Elements_conn_type & elements_conn, Elements_type & elements, CellList_type & cl_s)
    {
        for (int i = 0 ; i < elements_conn.size() ; i++)
        {
            Point<dim,T> xp_old = elements_conn.template get<1>(i);

            // relink
            auto NN_s = cl_s.getNNIterator(cl_s.getCell(xp_old));

#ifdef SE_CLASS1

            int num_surf = 0;

#endif

            while (NN_s.isNext())
            {
                auto q = NN_s.get();

                Point<dim,T> xs_nn = elements.getPos(q);

                if (xs_nn.distance(xp_old) < same_element_tollerance)
                {
                    // Found surface to connect

                    elements_conn.template get<0>(i) = q;
                    elements.template getProp<StartProp>(q) = elements_conn.template get<2>(i);

#ifdef SE_CLASS1
                    num_surf++;
#else
                    break;
#endif
                }

                ++NN_s;

            }

#ifdef SE_CLASS1
            if (num_surf == 0)
            {
                std::cout << __FILE__ << ":" << __LINE__ << " error found " << num_surf << "  elements in connectivity" << std::endl;
            }
#endif

        }
    }



    template<unsigned int GID,typename elements_type>
    void CreateGlobalIds(elements_type & elements)
    {
        auto n_local_ele = elements.size_local();

        auto & v_cl = create_vcluster();

        v_cl.scan_exclusive(n_local_ele);
        v_cl.execute();

        for (int i = 0 ; i < elements.size_local() ; i++)
        {
            elements.template getProp<GID>(i) = i+n_local_ele;
        }
    }

public:

    static constexpr int dims = dim;
    typedef T stype; 

    static constexpr int VolumeNumberOfSurfaces = TypeVolume::max_prop;
    static constexpr int SurfacesStart = TypeVolume::max_prop + 1;
    static constexpr int VolumeNumberOfDelTet = TypeVolume::max_prop + 2;
    static constexpr int DelTetStart = TypeVolume::max_prop + 3;
    static constexpr int VolumeGID = TypeVolume::max_prop + 4;
    static constexpr int SurfaceGID = TypeSurface::max_prop;
    static constexpr int SurfaceNumberOfEdges = TypeSurface::max_prop + 1;
    static constexpr int EdgesStart = TypeSurface::max_prop + 2;
    static constexpr int SurfacesMarking = TypeSurface::max_prop + 3;
    static constexpr int VerticesStart = TypeEdge::max_prop;
    static constexpr int EdgeMarking = TypeEdge::max_prop + 1;
    static constexpr int EdgeNumberOfVertices = TypeVertex::max_prop + 2;
    static constexpr int EdgeGID = TypeVertex::max_prop + 3;
    static constexpr int VertexMarking = TypeVertex::max_prop;
    static constexpr int VertexGID = TypeVertex::max_prop + 1;

    static constexpr int DelTelIndices = 0;
    static constexpr int DelTelVertId = 1;

    PolyMesh(Box<dim,T> & box, size_t (& bc)[dim], Ghost<dim,double> & ghost)
    :Volumes(0,box,bc,ghost),Surfaces(0,box,bc,ghost),Edges(0,box,bc,ghost),Vertices(0,box,bc,ghost)
    {}

    template<unsigned int ... prp>
    void ghost_get_volumes()
    {
        Volumes.template ghost_get<VolumeGID,prp ...>();
    }


    auto getDomainVolumeIterator() -> decltype(Volumes.getDomainIterator())
    {
        return Volumes.getDomainIterator();
    }

    template<typename lamb_type>
    void ForAllSurfaceEdges(int surf_index , lamb_type lamb)
    {
        auto n_edges = Surfaces.template getProp<SurfaceNumberOfEdges>(surf_index);
        auto edges_start = Surfaces.template getProp<EdgesStart>(surf_index);

        for (int j = 0 ; j < n_edges ; j++)
        {
            auto e = surface_edges_connectivity.template get<0>(j+edges_start);

            lamb(e);
        }
    }

    template<typename lamb_type>
    void ForAllVolumeSurfaces(size_t vol_index , lamb_type lamb)
    {
        auto n_faces = Volumes.template getProp<VolumeNumberOfSurfaces>(vol_index);
        auto faces_start = Volumes.template getProp<SurfacesStart>(vol_index);

        for (int j = 0 ; j < n_faces ; j++)
        {
            auto s = volume_surfaces_connectivity.template get<0>(j+faces_start);

            lamb(s,j+faces_start);
        }
    }

    /*! \brief Number of local Volumes
     *
     * \return the number of volumes
     * 
     */
    size_t numberOfLocalVolumes()
    {
        return Volumes.size_local();
    }    

    /*! \brief Add a Volume
     *
     * \param p Point where to add the volume
     * 
     */
    void addVolume(const Point<dim,T> & p)
    {
        Volumes.add();
        for (int i = 0 ; i < dim ; i++)
        {Volumes.getLastPos()[i] = p[i];}
    }

    /*! \brief Create the Voronoi tesellation
     *
     * 
     */
    void createVoronoi()
    {
        if (dim != 3)
        {std::cout << __FILE__ << ":" << __LINE__ << " Error voronoi tesellation works only in 3D" << std::endl;}

        unsigned int i,j;
        int id,nx,ny,nz;
        double x,y,z;

        // Create the Voronoi with Voro++ or Quentin Code

        voro::voronoicell_neighbor c;
        std::vector<int> neigh,f_vert;
        std::vector<double> v;
 
        auto domain = Volumes.getDecomposition().getProcessorBounds();

        Ghost<3,double> g = Volumes.getDecomposition().getGhost();

        for (int i = 0 ; i < dim ; i++)
        {
            if (Volumes.getDecomposition().periodicity(i) == NON_PERIODIC)
            {
                g.setLow(i,0.0);
                g.setHigh(i,0.0);
            }
        }

        domain.enlarge(g);

        // Create a pre-container class to import the input file and guess the
        // best computational grid size to use.
        voro::pre_container pcon(domain.getLow(0),domain.getHigh(0),
                                 domain.getLow(1),domain.getHigh(1),
                                 domain.getLow(2),domain.getHigh(2),false,false,false);

        // Randomly add particles into the container
        i = 0;
        auto it = Volumes.getDomainAndGhostIterator();
        while (it.isNext())
        {
            auto c = it.get();

            pcon.put(i,Volumes.getPos(c)[0],
                      Volumes.getPos(c)[1],
                      Volumes.getPos(c)[2]);

            ++i;
            ++it;
        }

        pcon.guess_optimal(nx,ny,nz);
 
        // Set up the container class and import the particles from the
        // pre-container
        voro::container con(domain.getLow(0),domain.getHigh(0),
                                 domain.getLow(1),domain.getHigh(1),
                                 domain.getLow(2),domain.getHigh(2),nx,ny,nz,false,false,false,8);
        pcon.setup(con);

        Surfaces.clear();
        Edges.clear();
        Vertices.clear();
        volume_surfaces_connectivity.clear();
        edge_vertices_connectivity.clear();
        surface_edges_connectivity.clear();

        // Open the output files
        // FILE *fp4=voro::safe_fopen("polygons4_v.pov","w"),
        //      *fp5=voro::safe_fopen("polygons5_v.pov","w"),
        //      *fp6=voro::safe_fopen("polygons6_v.pov","w");

        // Loop over all particles in the container and compute each Voronoi
        // cell

        voro::c_loop_all cl(con);
        if(cl.start())
        {
            do
            {
                cl.pos(x,y,z);

                Point<dim,T> p({x,y,z});
                if (!Volumes.getDecomposition().getProcessorBounds().isInside(p))   {continue;}

                if(con.compute_cell(c,cl)) 
                {
                    id=cl.pid();
        
                    // Gather information about the computed Voronoi cell
                    c.neighbors(neigh);
                    c.face_vertices(f_vert);
                    c.vertices(x,y,z,v);

                    // OpenFPM does not allow points outside the domain + ghost
                    // (it delete the point)
                    for (int i = 0 ; i < v.size() ; i++)
                    {
                        for (int j = 0 ; j < 3 ; j++)
                        {
                            if (v[3*i + j] < domain.getLow(j))
                            {v[3*i + j] = domain.getLow(j);}
                            if (v[3*i + j] > domain.getHigh(j))
                            {v[3*i + j] = domain.getHigh(j);}
                        }
                    }
                    

                    Volumes.template getProp<VolumeNumberOfSurfaces>(id) = neigh.size();
                    Volumes.template getProp<SurfacesStart>(id) = volume_surfaces_connectivity.size();

                    // Loop over all faces of the Voronoi cell
                    for( i = 0,j = 0; i < neigh.size() ; i++) 
                    {
                        Point<dim,T> previous;
                        volume_surfaces_connectivity.add();
                        volume_surfaces_connectivity.last().template get<0>() = Surfaces.size_local();

                        int k,l,n=f_vert[j];

                        Surfaces.add();
                        int last_surf = Surfaces.size_local() - 1;
                        Surfaces.template getLastProp<SurfaceNumberOfEdges>() = n;
                        Surfaces.template getLastProp<EdgesStart>() = surface_edges_connectivity.size();

                        Point<dim,T> center = {0.0,0.0,0.0};
                        Point<dim,T> first;
                        
                        int begining_vertex = Vertices.size_local();

                        for(k = 0; k < n; k++) 
                        {
                            l=3*f_vert[j+k+1];
                            center[0] += v[l];
                            center[1] += v[l+1];
                            center[2] += v[l+2];

                            Vertices.add();
                            Vertices.getLastPos()[0] = v[l];
                            Vertices.getLastPos()[1] = v[l+1];
                            Vertices.getLastPos()[2] = v[l+2];

                            if (k != 0)
                            {
                                surface_edges_connectivity.add();
                                surface_edges_connectivity.last().template get<0>() = Edges.size_local();

                                Edges.add();
                                Edges.getLastPos()[0] = (previous[0] + v[l]) / 2.0;
                                Edges.getLastPos()[1] = (previous[1] + v[l+1]) / 2.0;
                                Edges.getLastPos()[2] = (previous[2] + v[l+2]) / 2.0;

                                Edges.template getLastProp<EdgeNumberOfVertices>() = 2;

                                Edges.template getLastProp<VerticesStart>() = edge_vertices_connectivity.size();

                                edge_vertices_connectivity.add();
                                edge_vertices_connectivity.last().template get<0>() = Vertices.size_local() - 2;
                            
                                edge_vertices_connectivity.add();
                                edge_vertices_connectivity.last().template get<0>() = Vertices.size_local() - 1;
                            }
                            else
                            {
                                first[0] = v[l];
                                first[1] = v[l+1];
                                first[2] = v[l+2];
                            }

                            previous = Point<dim,T>({v[l],v[l+1],v[l+2]});
                        }
                        l=3*f_vert[j+n];

                        surface_edges_connectivity.add();
                        surface_edges_connectivity.last().template get<0>() = Edges.size_local();

                        Edges.add();
                        Edges.getLastPos()[0] = (first[0] + v[l]) / 2.0;
                        Edges.getLastPos()[1] = (first[1] + v[l+1]) / 2.0;
                        Edges.getLastPos()[2] = (first[2] + v[l+2]) / 2.0;

                        Edges.template getLastProp<VerticesStart>() = edge_vertices_connectivity.size();
                        Edges.template getLastProp<EdgeNumberOfVertices>() = 2;

                        edge_vertices_connectivity.add();
                        edge_vertices_connectivity.last().template get<0>() = begining_vertex;
                            
                        edge_vertices_connectivity.add();
                        edge_vertices_connectivity.last().template get<0>() = Vertices.size_local() - 1;

                        Surfaces.getLastPos()[0] = center[0] / n;
                        Surfaces.getLastPos()[1] = center[1] / n;
                        Surfaces.getLastPos()[2] = center[2] / n;

                        // Skip to the next entry in the face vertex list
                        j+=f_vert[j]+1;
                    }
                }

            } while (cl.inc());
        }

        // Create global IDs
        CreateGlobalIds<VolumeGID>(Volumes);
        CreateGlobalIds<SurfaceGID>(Surfaces);
        CreateGlobalIds<EdgeGID>(Edges);
        CreateGlobalIds<VertexGID>(Vertices);

        // Make unique

        int n_volumes = Volumes.size_local();
        double r_cut = pow(domain.getVolume() / n_volumes, 1.0 / dim) / 8;
        auto cl_vol = Volumes.getCellList(r_cut);

        typename std::remove_reference<decltype(Surfaces.getPosVector())>::type positions_unique;
        typename std::remove_reference<decltype(Surfaces.getPropVector())>::type properties_surfaces_unique;


        //////// SURFACE FIX ///////////////////////////////

        FillOldPositionConnectivity_element<EdgesStart>(volume_surfaces_connectivity,Surfaces);

        Surfaces.map();
        Surfaces.template ghost_get<SurfaceGID>();
        auto cl_s = Surfaces.getCellList(r_cut);

        unify_elements<SurfacesMarking,SurfaceGID>(Surfaces,cl_s,positions_unique,properties_surfaces_unique);

        FixIndexes<VolumeNumberOfSurfaces,SurfacesStart,SurfacesMarking>(Volumes,Surfaces,volume_surfaces_connectivity);

        int sz = positions_unique.size();
        Surfaces.getPosVector().swap(positions_unique);
        Surfaces.getPropVector().swap(properties_surfaces_unique);
        Surfaces.resize(sz);

        // Now we relink surfaces

        Surfaces.template ghost_get<SurfaceNumberOfEdges>();
        cl_s = Surfaces.getCellList(r_cut);

        relink_elements<EdgesStart>(volume_surfaces_connectivity,Surfaces,cl_s);

        FillOldPositionConnectivity_element<VerticesStart>(surface_edges_connectivity,Edges);
        typename std::remove_reference<decltype(Edges.getPropVector())>::type properties_edge_unique;


        // Make edge unique

        Edges.map();
        Edges.template ghost_get<EdgeGID>();
        auto cl_e = Edges.getCellList(r_cut);

        unify_elements<EdgeMarking,EdgeGID>(Edges,cl_e,positions_unique,properties_edge_unique);

        FixIndexes<SurfaceNumberOfEdges,EdgesStart,EdgeMarking>(Surfaces,Edges,surface_edges_connectivity);

        sz = positions_unique.size();
        Edges.getPosVector().swap(positions_unique);
        Edges.getPropVector().swap(properties_edge_unique);
        Edges.resize(sz);

        Edges.template ghost_get<EdgeNumberOfVertices>();
        cl_e = Edges.getCellList(r_cut);

        relink_elements<VerticesStart>(surface_edges_connectivity,Edges,cl_e);

        // Make vertex unique

        typename std::remove_reference<decltype(Vertices.getPropVector())>::type properties_vertices_unique;
        FillOldPositionConnectivity_element<VertexMarking>(edge_vertices_connectivity,Vertices);

        Vertices.map();
        Vertices.template ghost_get<VertexGID>();
        auto cl_v = Vertices.getCellList(r_cut);

        unify_elements<VertexMarking,VertexGID>(Vertices,cl_v,positions_unique,properties_vertices_unique);

        FixIndexes<EdgeNumberOfVertices,VerticesStart,VertexMarking>(Edges,Vertices,edge_vertices_connectivity);

        sz = positions_unique.size();
        Vertices.getPosVector().swap(positions_unique);
        Vertices.getPropVector().swap(properties_vertices_unique);
        Vertices.resize(sz);

        Vertices.template ghost_get<>();
        cl_v = Vertices.getCellList(r_cut);

        relink_elements<VertexMarking>(edge_vertices_connectivity,Vertices,cl_v);

        // Calculate Delauny

        calc_delaunay_tetraedron<DelTetStart,VolumeNumberOfDelTet>(Tet_del);

        return;
    }

    /*! \brief Get tge Delauney tetraedron vector
     *
     * 
     */
    openfpm::vector<aggregate<int[4],int>> & getDelTetraedron()
    {
        return Tet_del;
    }

    /*! \brief toString
     *
     * convert the content of the structure into a human readable string
     * 
     */
    std::string toString()
    {
        std::stringstream content;
        auto it = Volumes.getDomainIterator();

        content << "VOLUMES: " << std::endl;

        while (it.isNext())
        {
            auto pk = it.get();
            auto p = pk.getKey();

            content << "    " << p << "  " << Point<dim,T>(Volumes.getPos(p)).toString() << std::endl;

            auto n_surf = Volumes.template getProp<VolumeNumberOfSurfaces>(p);
            auto start_conn_vol = Volumes.template getProp<SurfacesStart>(p);

            for (int j = 0 ; j < n_surf ; j++)
            {
                auto surf_id = volume_surfaces_connectivity.template get<0>(start_conn_vol + j);

                content << "    " << "    " << " Surface ID:" << surf_id << "  " << Point<dim,T>(Surfaces.getPos(surf_id)).toString() << std::endl;


                auto n_edge =  Surfaces.template getProp<SurfaceNumberOfEdges>(surf_id);
                auto start_conn_surf = Surfaces.template getProp<EdgesStart>(surf_id);

                for (int k = 0 ; k < n_edge ; k++)
                {
                    auto edge_id = surface_edges_connectivity.template get<0>(start_conn_surf + k);

                    content << "    " << "    " << "    " << " Edge id: " << edge_id << "  " << Point<dim,T>(Edges.getPos(edge_id)).toString() << std::endl;

                    auto vstart = Edges.template getProp<VerticesStart>(edge_id);

                    auto vid0 = edge_vertices_connectivity.template get<0>(vstart);
                    auto vid1 = edge_vertices_connectivity.template get<0>(vstart+1);

                    content << "    " << "    " << "    " <<  "    Vertex id 0: " << vid0 << "  " << Point<dim,T>(Vertices.getPos(vid0)).toString() << std::endl;
                    content << "    " << "    " << "    " <<  "    Vertex id 1: " << vid1 << "  " << Point<dim,T>(Vertices.getPos(vid1)).toString() << std::endl;
                }
            }

            ++it;
        }

        return content.str();
    }

    /*! \brief Return the domain where the polymesh is defined
     *
     * \return the domain where the polymesh is defined
     * 
     */
    auto getDomain() -> decltype(Volumes.getDecomposition().getDomain())
    {
        return Volumes.getDecomposition().getDomain();
    }

    /*! \brief Return the volume iterator 
     * 
     * \return the volume iterator
     * 
     */
    auto getVolumeDomainIterator() -> decltype(Volumes.getDomainIterator())
    {
        return Volumes.getDomainIterator();
    }

    /*! \brief Return the vector of the volumes
     *
     * \return the vector of the volumes
     * 
     */
    auto getVolumesDist() -> decltype(Volumes) &
    {
        return Volumes;
    }

    /*! \brief Return the volume position 
     *
     * \param i volume element
     * 
     * \return the volume position
     * 
     */
    auto getVolumePos(size_t i) -> decltype(Volumes.getPos(0))
    {
        return Volumes.getPos(i);
    }

    /*! \brief Return the volume property
     *
     * \param i volume element
     * 
     * \return the volume property
     * 
     */
    template<unsigned int prp>
    auto getVolumeProp(size_t i) -> decltype(Volumes.template getProp<prp>(0))
    {
        return Volumes.template getProp<prp>(i);
    }

    /*! \brief Return the volume global-id
     *
     * \param i volume element
     * 
     * \return the volume GID
     * 
     */
    auto getVolumeGID(size_t i) -> decltype(Volumes.template getProp<VolumeGID>(0))
    {
        return Volumes.template getProp<VolumeGID>(i);
    }

    /*! \brief Return the Vertex index of an edge
     *
     * \param e edge
     * \param v vertex id
     * 
     * \return the vertex index
     * 
     */
    auto getEdgeVertexID(size_t e, int v)
    {
        #ifdef SE_CLASS1

        if (v > 2)
        {std::cout << __FILE__ << ":" << __LINE__ << " error every edge has maximum two vertices" << std::endl;}

        #endif

        auto vstart = Edges.template getProp<VerticesStart>(e);

        return edge_vertices_connectivity.template get<0>(v+vstart);
    }

    auto getEdgeVertexPos(size_t e, int v) -> decltype(Vertices.template getPos(0))
    {
        #ifdef SE_CLASS1

        if (v > 2)
        {std::cout << __FILE__ << ":" << __LINE__ << " error every edge has maximum two vertices" << std::endl;}

        #endif

        auto vstart = Edges.template getProp<VerticesStart>(e);

        return Vertices.template getPos(edge_vertices_connectivity.template get<0>(v+vstart));
    }

    /*! \brief Return the surface position 
     *
     * \param i surface element
     * 
     * \return the surface position
     * 
     */
    auto getSurfacePos(size_t i) -> decltype(Surfaces.getPos(0))
    {
        return Surfaces.getPos(i);
    }

    /*! \brief Return the volumes 
     *
     * \return the volumes
     * 
     */
    auto getVolumes() -> decltype(Volumes.getPropVector())
    {
        return Volumes.getPropVector();
    }

    /*! \brief Return the surfaces
     *
     * \return the surfaces
     * 
     */
    auto getSurfaces() -> decltype(Surfaces.getPropVector())
    {
        return Surfaces.getPropVector();
    }

    /*! \brief Return the surface connectivity
     *
     *  \return the surface to edge connectivity
     * 
     */
    openfpm::vector<aggregate<int,T[3],int>> & getVolumesSurfacesConn()
    {
        return volume_surfaces_connectivity;
    }

    /*! \brief Return the edge connectivity
     *
     *  \return the edge to vertices connectivity
     * 
     */
    openfpm::vector<aggregate<int,T[3],int>> & getEdgesVerticesConn()
    {
        return edge_vertices_connectivity;
    }

    /*! \brief Return the
     *
     *  return the Volume to surface connectivity
     * 
     */
    openfpm::vector<aggregate<int,T[3],int>> & getSurfacesEdgesConn()
    {
        return surface_edges_connectivity;
    }

    /*! Return the Vertex position
     *
     * \param v index of the position
     * 
     */
    template<typename key_type>
    auto getVertexPos(key_type & v) -> decltype(Vertices.getPos(v))
    {
        return Vertices.getPos(v);
    }

    /*! \brief Return the edges
     *
     * \return the edges
     * 
     */
    auto getEdges() -> decltype(Edges.getPropVector())
    {
        return Edges.getPropVector();
    }

    /*! \brief Return the volumes 
     *
     * \return the volumes
     * 
     */
    auto getSurfacesPosition() -> decltype(Surfaces.getPosVector())
    {
        return Surfaces.getPosVector();
    }

    /*! \brief Return the volumes 
     *
     * \return the volumes
     * 
     */
    auto getVolumesPosition() -> decltype(Volumes.getPosVector())
    {
        return Volumes.getPosVector();
    }

    /*! \brief Return the vertices 
     *
     * \return the vertices
     * 
     */
    auto getVerticesPosition() -> decltype(Vertices.getPosVector())
    {
        return Vertices.getPosVector();
    }

    /*! \brief Return the volume of a cell
     *
     */
    T getVolume(size_t vol_id)
    {
        T vt = 0;
        
        ForAllVolumeTetraedron(*this,vol_id, [&](Tetraedron<dim,T> & tet) {

            vt += tet.volume();

        });

        return vt;
    }


    /*! \brief check consistency of the connection
     * 
     */
    bool check_consistent()
    {
        bool consistent = true;

        {
        auto it = Volumes.getDomainIterator();

        while (it.isNext())
        {
            auto pk = it.get();
            auto p = pk.getKey();

            auto n_surf = Volumes.template getProp<VolumeNumberOfSurfaces>(p);
            auto start_conn_vol = Volumes.template getProp<SurfacesStart>(p);

            for (int j = 0 ; j < n_surf ; j++)
            {
                if (start_conn_vol + j >= volume_surfaces_connectivity.size())
                {
                    std::cout << __FILE__ << ":" << __LINE__ << " error a volume points to an invalid id of surface connectivity" << std::endl;

                    return false;
                }

                auto surf_id = volume_surfaces_connectivity.template get<0>(start_conn_vol + j);

                if (surf_id >= Surfaces.size_local_with_ghost())
                {
                    std::cout << __FILE__ << ":" << __LINE__ << " error, volume surface connectivity point to an inexistent surface" << std::endl;

                    return false;
                }

                auto n_edge =  Surfaces.template getProp<SurfaceNumberOfEdges>(surf_id);
                auto start_conn_surf = Surfaces.template getProp<EdgesStart>(surf_id);

                for (int k = 0 ; k < n_edge ; k++)
                {
                    if (start_conn_surf + k >= surface_edges_connectivity.size())
                    {
                        std::cout << __FILE__ << ":" << __LINE__ << " error a surface points to an invalid id of edge connectivity" << std::endl;

                        return false;
                    }

                    auto edge_id = surface_edges_connectivity.template get<0>(start_conn_surf + k);

                    if (edge_id >= Edges.size_local_with_ghost())
                    {
                        std::cout << __FILE__ << ":" << __LINE__ << " error, surface edge connectivity point to an inexistent edge" << std::endl;

                        return false;
                    }

                    auto vstart = Edges.template getProp<VerticesStart>(edge_id);

                    if (vstart+1 >= edge_vertices_connectivity.size())
                    {
                        std::cout << __FILE__ << ":" << __LINE__ << " error a surface points to an invalid id of edge connectivity" << std::endl;

                        return false;
                    }

                    auto vid0 = edge_vertices_connectivity.template get<0>(vstart);
                    auto vid1 = edge_vertices_connectivity.template get<0>(vstart+1);
                }
            }

            ++it;
        }
        }

        // first we check that all surfaces are surfaces

        auto it = getDomainVolumeIterator();

        std::unordered_map<size_t,size_t> all_verts;
        openfpm::vector<int> indices;

        while(it.isNext())
        {
            auto v = it.get().getKey();

            ForAllVolumeSurfaces(v,[&](unsigned int s, unsigned int cs){
            
                auto n_e = Surfaces.template getProp<SurfaceNumberOfEdges>(s);
                auto start_e = Surfaces.template getProp<EdgesStart>(s);

                if (n_e < 3)
                {
                    std::cout << __FILE__ << ":" << __LINE__ << " error, found a face with lower than 3 vertices" << std::endl;

                    consistent = false;

                    return;
                }

                all_verts.clear();

                for (int j = 0 ; j < n_e ; j++)
                {
                    auto e = surface_edges_connectivity.template get<0>(start_e + j);
                    auto nv = Edges.template getProp<EdgeNumberOfVertices>(e);
                    auto start_v = Edges.template getProp<VerticesStart>(e);

                    if (nv != 2)
                    {consistent = false;}

                    auto v = edge_vertices_connectivity.template get<0>(start_v);
                    all_verts[v]++;
                    v = edge_vertices_connectivity.template get<0>(start_v+1);
                    all_verts[v]++;
                }

                indices.clear();

                // Check edge create a closed loop
                for (std::pair<size_t, size_t> element : all_verts)
                {
                    if (element.second != 2)
                    {
                        consistent = false;
                        return;
                    }

                    indices.add(element.first);
                }

                if (indices.size() < 3)
                {
                    std::cout << __FILE__ << ":" << __LINE__ << " error, found a face with lower than 3 vertices" << std::endl;

                    consistent = false;

                    return;
                }

                Point<dim,T> p1 = Vertices.template getPos(indices.template get<0>(0));
                Point<dim,T> p2 = Vertices.template getPos(indices.template get<0>(1));
                int p3_id = 2;

                Plane<dim,T> pln;

                bool unique = false;

                while (unique == false && p3_id < indices.size())
                {
                    Point<dim,T> p3 = Vertices.template getPos(indices.template get<0>(p3_id));

                    unique = pln.construct(p1,p2,p3);

                    p3_id++;
                }

                // We check all point lie on the plane

                for (int i = 0 ; i < indices.size() ; i++)
                {
                    Point<dim,T> p = Vertices.template getPos(indices.template get<0>(i));

                    auto res = pln.residual(p);

                    if (res >= 1e-5)
                    {
                        std::cout << __FILE__ << ":" << __LINE__ << " error, a face is not on a plane" << std::endl;
                        consistent = false;

                        return;
                    }
                }

            });

            ++it;
        }

        if (consistent == false)
        {return false;}

        double vol = 0.0;

        {
        auto it = getDomainVolumeIterator();

        while(it.isNext())
        {
            auto v = it.get().getKey();

            auto single_vol = getVolume(v);

            vol += single_vol;

            ++it;
        }
        }

        auto TotVol = getDomain().getVolume();

        if (fabs(TotVol - vol) / TotVol > 0.001)
        {
            std::cout << __FILE__ << ":" << __LINE__ << " error tesellation does not cover the entire space" << std::endl;

            return false;
        }

        // Check unique elements

        auto domain = Volumes.getDecomposition().getProcessorBounds();
        domain.enlarge(Volumes.getDecomposition().getGhost());

        int n_volumes = Volumes.size_local();
        double r_cut = pow(domain.getVolume() / n_volumes, 1.0 / dim) / 8;

        auto cl_vol = Volumes.getCellList(r_cut);
        auto cl_s = Surfaces.getCellList(r_cut);
        auto cl_e = Edges.getCellList(r_cut);
        auto cl_v = Vertices.getCellList(r_cut);

        if (check_uniqueness<>(Volumes,cl_vol) == false) {return false;}
        if (check_uniqueness<>(Surfaces,cl_s) == false) {return false;}
        if (check_uniqueness<>(Edges,cl_e) == false) {return false;}
        if (check_uniqueness<>(Vertices,cl_v) == false) {return false;}

        return true;
    }

    mutable tsl::hopscotch_map<int, int> map;

    tsl::hopscotch_map<int, int> & getTempUnorderedMap()
    {
        return map;
    }

    mutable openfpm::vector<int> indices;

    openfpm::vector<int> & getTempVectorInteger()
    {
        return indices;
    }

    /*! \brief write on VTK
     *
     */
    void write(std::string file)
    {
        // VTK writer

        typedef typename std::remove_reference<decltype(Volumes.getPosVector())>::type vector_of_volumes_position;
        typedef typename std::remove_reference<decltype(Volumes.getPropVector())>::type vector_of_volumes;
        typedef typename std::remove_reference<decltype(Surfaces.getPropVector())>::type vector_of_surfaces;
        typedef typename std::remove_reference<decltype(Edges.getPropVector())>::type vector_of_edges;
        typedef typename std::remove_reference<decltype(Vertices.getPropVector())>::type vector_of_vertices;

        auto & Volumes_prop = Volumes.getPropVector();
        auto & Volumes_pos = Volumes.getPosVector();
        auto & Surfaces_prop = Surfaces.getPropVector();
        auto & Surfaces_pos = Surfaces.getPosVector();
        auto & Edges_prop = Edges.getPropVector();
        auto & Vertices_prop = Vertices.getPropVector();
        auto & Vertices_pos = Vertices.getPosVector();

		// Create a writer and write
		VTKWriter<boost::mpl::vector<vector_of_volumes,
									 vector_of_surfaces,
									 vector_of_edges,
									 vector_of_vertices,
									 vector_of_volumes_position
									 >,VTK_POLYMESH> vtk_ply;

        openfpm::vector<std::string> null;

        vtk_ply.add(Volumes_prop,Volumes_pos,
                    Surfaces_prop,Surfaces_pos,
                    Edges_prop,
                    Vertices_prop,Vertices_pos,
                    volume_surfaces_connectivity,surface_edges_connectivity,edge_vertices_connectivity,
                    Volumes.size_local(),Surfaces.size_local(),Edges.size_local(),Vertices.size_local());

        vtk_ply.write(file + ".vtu",null,null,null,null,"vtk output","",file_type::ASCII);
    }

};


#endif