//
// Created by aryaman on 6/2/21.
//

#ifndef OPENFPM_PDATA_GRID_TO_INVIS_HPP
#define OPENFPM_PDATA_GRID_TO_INVIS_HPP

#include "InVis.hpp"

namespace grid_to_InVis {

    template<bool is_vector>
    struct vis_code
    {
        template<unsigned int prop, typename vector_type, typename key_type>
        static float execute(vector_type * vd, key_type & key, vector_visualize vect_vis)
        {
            return (float) vd->template get<prop>(key);
        };
    };

    template<>
    struct vis_code<true>
    {
        template<unsigned int prop, typename vector_type, typename key_type>
        static inline float execute(vector_type * vd, key_type & key, vector_visualize vect_vis)
        {

            if(vect_vis == vector_visualize::x_)
            {return (float) vd->template get<prop>(key)[0];}
            else if(vect_vis == vector_visualize::y_)
            {return (float) vd->template get<prop>(key)[1];}
            else if(vect_vis == vector_visualize::z_)
            {return (float) vd->template get<prop>(key)[2];}
            else //magnitude
            {
                return (float) sqrt (vd->template get<prop>(key)[0] * vd->template get<prop>(key)[0] +
                                     vd->template get<prop>(key)[1] * vd->template get<prop>(key)[1] +
                                     vd->template get<prop>(key)[2] * vd->template get<prop>(key)[2]);
            }
        };
    };

    void setup_shared_mem(size_t shm_rank)
    {
        std::cout<<"Shm rank is "<<shm_rank<<std::endl;
        if (global_option == init_options::in_situ_visualization)
        {
            std::fstream fs;
            if(shm_rank == 0)
            {
                //if file used for initializing shm location indicating data type does not exist, create it
                fs.open(datatype_path);
                if(fs.fail())
                {
                    fs.open(datatype_path, std::ios::out);
                }
                fs.close();

                // specify to visualization that grid data needs to be rendered
                std::cout<<"Specifying data type"<<std::endl;
                dtype_flag = create_shmanager().create(datatype_path, 0);
                int * ptr = (int *)create_shmanager().alloc(dtype_flag, sizeof(int));
                *ptr = 1; // 1: Grid Data, 2: Particle Data
            }

            //for this process, if the file path that will be used for shared memory initialization of grid
            // locations does not exist, create it
            fs.open(grid_shm_path + std::to_string(shm_rank));
            if(fs.fail())
            {
                fs.open(grid_shm_path + std::to_string(shm_rank), std::ios::out);
            }
            fs.close();
        }

    }

    template<unsigned int dim, unsigned int prop, typename St, typename grid_dist_id_type, typename vcluster_type>
    void visualize(grid_dist_id_type &grid, scale_vis scale, St low, St high, vector_visualize vect_vis, vcluster_type &v_cl)
    {
        typedef typename grid_dist_id_type::d_grid::value_type T;
        auto &gdb_ext = grid.getLocalGridsInfo();
        if (global_option != init_options::in_situ_visualization)
        {

            std::cerr<<__FILE__<<":"<<__LINE__<<" ERROR: In the 'visualize' method. You need to call openfpm_init with 'init_options::in_situ_visualization'. Example: 'openfpm_init(&argc,&argv, init_options::in_situ_visualization);'"<<std::endl;
            return;
        }

        if(dim != 3)
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" ERROR: In the 'visualize' method. Currently only three-dimensional grids can be visualized."<<std::endl;
            return;
        }

        typedef typename boost::mpl::at<typename T::type,boost::mpl::int_<prop>>::type prop_type;
        int rank_type = std::rank<prop_type>::value;
        int vec_dim = 0;


        if(rank_type == 1) //The property is a vector
        {
            vec_dim = std::extent<prop_type, 0>::value;
        }

        bool resize_happened = false; // TODO: get information on whether resize has taken place

        if(resize_happened)
        {
            for(size_t i = Vis_ptrs.size()-1; i>=0; i--)
            {
                create_shmanager().destroy(hgrids.get(i));
            }
            Vis_ptrs.clear();
            create_shmanager().destroy(hgdb);
            create_shmanager().destroy(dtype_flag);
        }

//        grid_dist_id<3, St, aggregate<unsigned short>,Decomposition>* Vis_new = (grid_dist_id<3, St, aggregate<unsigned short>,Decomposition> *) Vis_new_void;
//        unsigned short * Vis_new = (unsigned short *) Vis_new_void;
//        openfpm::vector<unsigned short *> Vis_ptrs = Vis_ptrs_void;
//        long * Vis_header = Vis_header_void;
        if(Vis_ptrs.size() == 0) // if no vis grids have been created yet
        {
            //create the files required for System V shared memory
            setup_shared_mem(v_cl.shmRank());

//            Vis_new = new grid_dist_id<3, St, aggregate<unsigned short>,Decomposition>(getDecomposition(),ginfo.getSize(),Ghost<dim,St>(0.0));
//            handle_shmem h = create_shmanager().create(grid_shm_path + std::to_string(v_cl.shmRank()), 0);

//            Vis_new = new grid_dist_id<3, St, aggregate<unsigned short>,Decomposition>(getDecomposition(),ginfo.getSize(),Ghost<dim,St>(0.0));
            for(size_t i = 0; i < gdb_ext.size(); i++)
            {
                // for each grid
                hgrids.add(create_shmanager().create(grid_shm_path + std::to_string(v_cl.shmRank()), (int)i+1));

                Vis_ptrs.add(create_shmanager().alloc(hgrids.get(hgrids.size()-1), gdb_ext.get(i).DBox.getVolumeKey()));
            }
//            Vis_new = (unsigned short *)create_shmanager().alloc(h, total_size);
            //TODO: delete Vis_new_void in destructor of grid_dist_id
//            Vis_new_void = (void *) Vis_new;
//            Vis_new->to_shared_mem();
//            Vis_ptrs_void = Vis_ptrs;
        }
        if(Vis_header == nullptr)
        {
            hgdb = create_shmanager().create(grid_shm_path + std::to_string(v_cl.shmRank()), 0);

            // From each box, we need origin (3), grid low and high (6) and domain low and high (6)
            Vis_header = (long *) create_shmanager().alloc(hgdb, gdb_ext.size() * 15 * sizeof(long));
            size_t cnt = 0;
            for(int i = 0; i<gdb_ext.size(); i++)
            {
                Vis_header[cnt] = gdb_ext.get(i).GDbox.getLow(0);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).GDbox.getLow(1);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).GDbox.getLow(2);
                cnt++;

                Vis_header[cnt] = gdb_ext.get(i).GDbox.getHigh(0);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).GDbox.getHigh(1);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).GDbox.getHigh(2);
                cnt++;

                Vis_header[cnt] = gdb_ext.get(i).DBox.getLow(0);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).DBox.getLow(1);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).DBox.getLow(2);
                cnt++;

                Vis_header[cnt] = gdb_ext.get(i).DBox.getHigh(0);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).DBox.getHigh(1);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).DBox.getHigh(2);
                cnt++;

                Vis_header[cnt] = gdb_ext.get(i).origin.get(0);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).origin.get(1);
                cnt++;
                Vis_header[cnt] = gdb_ext.get(i).origin.get(2);
                cnt++;
            }
        }

        static grid_key_dx<3> star_stencil_3D[1] = {{0,0,0}};

        if(scale == fixed)
        {
            if(low >= high)
            {
                std::cerr<<__FILE__<<":"<<__LINE__<<" ERROR: In the 'visualize' method. You have chosen to visualize your grid with fixed scaling. Please provide a value of 'high' that is greater than 'low'"<<std::endl;
                return;
            }
        }
        else if(scale == global)
        {
            low = std::numeric_limits<St>::max();
            high = std::numeric_limits<St>::min();

            auto it1 = grid->getDomainIteratorStencil(star_stencil_3D);;

            while(it1.isNext())
            {
                auto Cp = it1.template getStencil<0>();

                float cur;
                cur = vis_code<std::rank<prop_type>::value == 1>::template execute<prop>(grid,Cp,vect_vis);

                if(cur > high) {high = cur;}
                if(cur < low) {low = cur;}

                ++it1;
            }

            v_cl.max(high);
            v_cl.min(low);
            v_cl.execute();
        }

//        // calculate the magnitude of velocity
//        auto it2 = this->getDomainIteratorStencil(star_stencil_3D);
////        auto it2_vis = Vis_new->getDomainIteratorStencil(star_stencil_3D);
//
//        size_t cnt = 0;
//
//        while (it2.isNext())
//        {
//            auto Cp = it2.template getStencil<0>();
//
////            auto Cp_vis = it2_vis.template getStencil<0>();
//
//            double cur = vis_code<std::rank<prop_type>::value == 1>::template execute<prop>(this,Cp,vect_vis);
//
//            double scaled = (cur / (high - low)) * 65535;
//            // copy
////            Vis_new->template get<0>(Cp_vis) = (unsigned short)(scaled);
//            Vis_new[cnt] = (unsigned short)(scaled);
//
//            ++it2;
//            ++cnt;
////            ++it2_vis;
//        }

        for(size_t i = 0; i < grid.getN_loc_grid(); i++)
        {
            grid_key_dx <dim> center[1];
            center[0].zero();

            auto &src = grid.get_loc_grid(i);
            auto it = grid->get_loc_grid_iterator_stencil(i, center);
            size_t cnt = 0;

            while(it.isNext())
            {
                auto Cp = it.template getStencil<0>();

                float cur = vis_code<std::rank<prop_type>::value == 1>::template execute<prop>(src,Cp,vect_vis);

                float scaled = (cur / (high - low)) * 65535;

                Vis_ptrs.get(i)[cnt] = (unsigned short)(scaled);

                ++it;
                ++cnt;
            }
        }
    }
};


#endif //OPENFPM_PDATA_GRID_TO_INVIS_HPP
