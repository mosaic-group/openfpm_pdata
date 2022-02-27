#ifndef __POLYMESH_HPP__
#define __POLYMESH_HPP__

#include <Vector/vector_dist.hpp>
#include "voro++.hh"

template<unsigned int dim, typename T, typename TypeVolume, typename TypeSurface, typename TypeEdge, typename TypeVertex>
class PolyMesh
{
    static constexpr int SurfaceVolume = TypeSurface::max_prop;
    static constexpr int SurfaceNumberOfEdges = TypeSurface::max_prop + 1;
    static constexpr int EdgesStart = TypeSurface::max_prop + 2;
    static constexpr int EdgeVertices = TypeEdge::max_prop;
    static constexpr int EdgeSurfaces = TypeEdge::max_prop + 1;
    static constexpr int EdgeMarking = TypeEdge::max_prop + 2;
    static constexpr int VertexMarking = TypeVertex::max_prop;

    // Transform volumes
    typedef typename AggregateAppend<int,TypeVolume>::type TypeVolume_plus;
    typedef typename AggregateAppend<int,TypeVolume_plus>::type TypeVolume_plusplus;

    // Transform surfaces
    typedef typename AggregateAppend<int[2],TypeSurface>::type TypeSurface_plus;
    typedef typename AggregateAppend<int,TypeSurface_plus>::type TypeSurface_plusplus;
    typedef typename AggregateAppend<int,TypeSurface_plusplus>::type TypeSurface_plusplusplus;

    // Transform edges
    typedef typename AggregateAppend<int[2],TypeEdge>::type TypeEdge_plus;
    typedef typename AggregateAppend<int[2],TypeEdge_plus>::type TypeEdge_plusplus;
    typedef typename AggregateAppend<int,TypeEdge_plusplus>::type TypeEdge_plusplusplus;

    // Transform Vertex
    typedef typename AggregateAppend<int,TypeVertex>::type TypeVertex_plus;


    //! List of cells
    vector_dist<dim,T,TypeVolume_plusplus> Volumes;

    //! List of cells surfaces
    vector_dist<dim,T,TypeSurface_plusplusplus> Surfaces;

    //! List of edges
    vector_dist<dim,T,TypeEdge_plusplusplus> Edges;

    //! List of vertices
    vector_dist<dim,T,TypeVertex_plus> Vertices;

    //! Tollerance (in distance) for which two edges/vertices are considered the same
    T same_element_tollerance = 1e-5;

    template<unsigned int not_unique_id, unsigned int unique_id, typename Ele_type, typename Ele_type_up>
    void FixIndexes(Ele_type_up & Elements_up,
                    Ele_type & Elements)
    {
        auto it_e2 = Elements_up.getDomainIterator();

        while (it_e2.isNext())
        {
            auto pk = it_e2.get();
            auto p = pk.getKey();

            int element_down_id_not_unique = Elements_up.template getProp<not_unique_id>(pk);
            int unique_elements_id = Elements.template getProp<unique_id>(element_down_id_not_unique);

            Elements_up.template getProp<not_unique_id>(pk) = unique_elements_id;

            ++it_e2;
        }
    }

    template<unsigned int not_unique_id, unsigned int unique_id, typename Ele_type, typename Ele_type_up>
    void FixIndexes_i2(Ele_type_up & Elements_up,
                    Ele_type & Elements)
    {
        auto it_e2 = Elements_up.getDomainIterator();

        while (it_e2.isNext())
        {
            auto pk = it_e2.get();
            auto p = pk.getKey();

            int element_down_id_not_unique_0 = Elements_up.template getProp<not_unique_id>(pk)[0];
            int element_down_id_not_unique_1 = Elements_up.template getProp<not_unique_id>(pk)[1];
            int unique_elements_id_0 = Elements.template getProp<unique_id>(element_down_id_not_unique_0);
            int unique_elements_id_1 = Elements.template getProp<unique_id>(element_down_id_not_unique_1);

            Elements_up.template getProp<not_unique_id>(pk)[0] = unique_elements_id_0;
            Elements_up.template getProp<not_unique_id>(pk)[1] = unique_elements_id_1;

            ++it_e2;
        }
    }

    /* \brief 
     *
     * 
     */
    template<unsigned int prp, typename Ele_type, typename cl_type>
    void unify_elements(Ele_type & Elements, cl_type & cl_e,
                        typename std::remove_reference<decltype(Elements.getPosVector())>::type & positions_unique,
                        typename std::remove_reference<decltype(Elements.getPropVector())>::type & properties_unique)
    {
        auto it_e = Elements.getDomainIterator();

        while (it_e.isNext())
        {
            auto p = it_e.get();
            Point<dim,T> xp = Elements.getPos(p);

            auto it_nn = cl_e.getNNIterator(cl_e.getCell(xp));

            while (it_nn.isNext())
            {
                auto q = it_nn.get();

                if (p.getKey() == q)   {++it_nn; continue;}

                Point<dim,T> xq = Elements.getPos(p);

                if (xp.distance(xq) < same_element_tollerance)
                {
                    Elements.template getProp<prp>(q) = (p.getKey() < q)?p.getKey():q;
                    Elements.template getProp<prp>(p) = (p.getKey() < q)?p.getKey():q;
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

            if (p == Elements.template getProp<prp>(p))
            {
                positions_unique.add();
                properties_unique.add();
                int last = positions_unique.size()-1;

                positions_unique.get(last) = Elements.getPosVector().get(p);
                properties_unique.get(last) = Elements.getPropVector().get(p);

                Elements.template getProp<prp>(p) = last;
            }

            ++it_e2;
        }
    }

public:

    PolyMesh(Box<dim,T> & box, size_t (& bc)[dim], Ghost<dim,double> & ghost)
    :Volumes(0,box,bc,ghost),Surfaces(0,box,bc,ghost),Edges(0,box,bc,ghost),Vertices(0,box,bc,ghost)
    {}

    void addCell(const Point<dim,T> & p)
    {
        Volumes.add();
        for (int i = 0 ; i < dim ; i++)
        {Volumes.getLastPos()[i] = p[i];}
    }

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
 
        auto & domain = Volumes.getDecomposition().getProcessorBounds();

        // Create a pre-container class to import the input file and guess the
        // best computational grid size to use.
        voro::pre_container pcon(domain.getLow(0),domain.getHigh(0),
                                 domain.getLow(1),domain.getHigh(1),
                                 domain.getLow(2),domain.getHigh(2),true,false,false);

        // Randomly add particles into the container
        i = 0;
        auto it = Volumes.getDomainIterator();
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
                                 domain.getLow(2),domain.getHigh(2),nx,ny,nz,true,false,false,8);
        pcon.setup(con);

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
                if(con.compute_cell(c,cl)) 
                {
                    cl.pos(x,y,z);
                    id=cl.pid();
        
                    // Gather information about the computed Voronoi cell
                    c.neighbors(neigh);
                    c.face_vertices(f_vert);
                    c.vertices(x,y,z,v);

                    // Loop over all faces of the Voronoi cell
                    for( i = 0,j = 0; i < neigh.size();i++) 
                    {
                        // Draw all quadrilaterals, pentagons, and hexagons.
                        // Skip if the neighbor information is smaller than
                        // this particle's ID, to avoid double counting. This
                        // also removes faces that touch the walls, since the
                        // neighbor information is set to negative numbers for
                        // these cases.
                        if(neigh[i]>id) 
                        {
                            Surfaces.add();

                            Point<dim,T> center = {0.0,0.0,0.0};
                            Point<dim,T> previous;
                            int k,l,n=f_vert[j];

                            for(k = 0; k < n; k++) 
                            {
                                l=3*f_vert[j+k+1];
                                center[0] += v[l];
                                center[1] += v[l+1];
                                center[2] += v[l+2];

                                previous = Point<dim,T>({v[l],v[l+1],v[l+2]});

                                Vertices.add();
                                Vertices.getLastPos()[0] = v[l];
                                Vertices.getLastPos()[1] = v[l+1];
                                Vertices.getLastPos()[2] = v[l+2];

                                if (k != 0)
                                {
                                    Edges.add();
                                    Edges.getLastPos()[0] = (previous[0] + v[l]) / 2.0;
                                    Edges.getLastPos()[1] = (previous[1] + v[l+1]) / 2.0;
                                    Edges.getLastPos()[2] = (previous[2] + v[l+2]) / 2.0;
                                }
                            }

                            Surfaces.getLastPos()[0] = center[0];
                            Surfaces.getLastPos()[1] = center[1];
                            Surfaces.getLastPos()[2] = center[2];
                        }
                        // Skip to the next entry in the face vertex list
                        j+=f_vert[j]+1;
                    }
                }
            } while (cl.inc());
        }

        // Make unique

        int n_volumes = Volumes.size_local();
        double r_cut = pow(domain.getVolume() / n_volumes, 1.0 / dim);

        typename std::remove_reference<decltype(Edges.getPosVector())>::type positions_unique;
        typename std::remove_reference<decltype(Edges.getPropVector())>::type properties_edge_unique;

        // Make edge unique

        auto cl_e = Edges.getCellList(r_cut);
        unify_elements<EdgeMarking>(Edges,cl_e,positions_unique,properties_edge_unique);

        FixIndexes<EdgesStart,EdgeMarking>(Surfaces,Edges);

        Edges.getPosVector().swap(positions_unique);
        Edges.getPropVector().swap(properties_edge_unique);

        // Done

        // Make vertex unique

        positions_unique.clear();

        typename std::remove_reference<decltype(Vertices.getPropVector())>::type properties_vertices_unique;

        auto cl_v = Vertices.getCellList(r_cut);
        unify_elements<VertexMarking>(Vertices,cl_v,positions_unique,properties_vertices_unique);

        FixIndexes_i2<EdgeVertices,VertexMarking>(Edges,Vertices);

        Vertices.getPosVector().swap(positions_unique);
        Vertices.getPropVector().swap(properties_vertices_unique);
    }

    /*! \brief write on VTK
     *
     */
    void write(std::string file)
    {
        // VTK writer

        
    }

};


#endif