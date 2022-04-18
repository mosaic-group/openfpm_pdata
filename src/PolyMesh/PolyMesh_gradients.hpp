#ifndef POLYMESH_GRADIENTS_HPP
#define POLYMESH_GRADIENTS_HPP

template<unsigned int dim, typename T>
T det(T (& mat)[dim][dim])
{
    Eigen::Map<Eigen::Matrix<T,dim,dim,Eigen::RowMajor> > m((T *)mat);

    return m.determinant();
}

template<unsigned int dim>
Point<dim,int> index_perm(int v_c)
{
    //makes a sequence of indices from 0 to 3 but not including index v_c
    //which must be between 0 and 3
    Point<dim,int> x({0,0,0});
    int i;
    int ind=0;
    for(i = 0 ; i < 4 ; i++)
    {
        if(i != v_c)
        {
            x[ind]=i;
            ind++;
        }
    }
    return x;
}


template<unsigned int dim, typename T>
void compute_D(T (& D)[dim][4],T (& dDm)[dim][dim][dim],int v_c,Point<dim,T> (& V)[4], T (& r2)[4],int (& One_array)[4])
{
    //computes the minor determinants for D[3][4], first dimension is for Dx,Dy,Dz and for each of these there are 4 3x3 determinants to compute
    //r2[4] is the column containing the norm squared of each vector in V, it is provided as input because it depends on which rows have been modified due to offsets
    //Also computes in dDm[i][j][k] the derivative of D[i][j+1] with respect to dimension k (D[i][0] has a zero derivative).
    //We also need the array One_array (1 or 0) telling us whihch row has been modified
    int i,j;
    //we need some indices array again
    //this is for the variations between D[0],D[1] and D[2]
    int y0[3]; y0[0]=1; y0[1]=0; y0[2]=0;
    int y1[3]; y1[0]=2; y1[1]=2; y1[2]=1;
    Point<dim,int> x=index_perm<dim>(v_c); //containing the indices of the other vertice than v_c
    // v_c= 0 -> x=(1,2,3)
    //      1 -> x=(0,2,3)
    //      2 -> x=(0,1,3)
    //      3 -> x=(0,1,2)
    T m[3][3];//matrix for minor determinants
    //compute D[i][j]
    int ro;
    //struct MatDoub m2=Mat_create(4, 4); //matrix for testing consistency by computing Dx,Dy,Dz directly
    for(i=0;i<3;i++)//This is to compute D[i][j]
    {
        //make the whole matrix D[i] (TEST)
        /*
         for(ro=0;ro<4;ro++)
         {
         m2.M[ro][0]=r2[ro];
         m2.M[ro][1]=V[ro].coo[y0[i]];
         m2.M[ro][2]=V[ro].coo[y1[i]];
         m2.M[ro][3]=1;
         }
         affiche_mat(&m2);
         det(&m2);
         */
        
        
        //D[i][0] is with V[][y0] V[][y1] 1
        //fill the matrix
        for(ro=0;ro<3;ro++)
        {
            m[ro][0]=V[x[ro]][y0[i]];
            m[ro][1]=V[x[ro]][y1[i]];
            m[ro][2]=One_array[x[ro]];
        }
        //affiche_mat(&m);
        D[i][0]=det(m);
        //D[i][1] is with r2 V[y1] 1
        for(ro=0;ro<3;ro++)
        {
            m[ro][0]=r2[x[ro]];
            m[ro][1]=V[x[ro]][y1[i]];
            m[ro][2]=One_array[x[ro]];
        }
        //affiche_mat(&m);
        D[i][1]=det(m);
        //D[i][2] is r2 V[y0] 1
        for(ro=0;ro<3;ro++)
        {
            m[ro][0]=r2[x[ro]];
            m[ro][1]=V[x[ro]][y0[i]];
            m[ro][2]=One_array[x[ro]];
        }
        //affiche_mat(&m);
        D[i][2]=det(m);
        //D[i][3] i r2 V[y0] V[y1]
        for(ro=0;ro<3;ro++)
        {
            m[ro][0]=r2[x[ro]];
            m[ro][1]=V[x[ro]][y0[i]];
            m[ro][2]=V[x[ro]][y1[i]];
        }
        //affiche_mat(&m);
        D[i][3]=det(m);
    }
    
    //now we compute the dDm[i][j][k]
    int sgn,k;
    for(k=0;k<3;k++)
    {
        for(i=0;i<3;i++)
        {
            //D[i][0] has zero derivative
            //D[i][1] might have a derivative if some rows have been modified due to offsets
            //Same for D[i][2], the only difference being to replace y1[i] by y0[i]
            //For D[i][3], it's like D[i][2] but replacing One_array by V.coo[y1[i]]
            dDm[i][0][k]=0;
            dDm[i][1][k]=0;
            dDm[i][2][k]=0;
            sgn=1; //sign for expanding determinant
            for(j=0;j<3;j++)
            {
                if(One_array[x[j]]==0)//detecting modified row
                {
                    //adding the corresponding term
                    dDm[i][0][k]=dDm[i][0][k]+sgn*2*V[x[j]][k]*(V[x[y0[j]]][y1[i]]*One_array[x[y1[j]]]-V[x[y1[j]]][y1[i]]*One_array[x[y0[j]]]);
                    dDm[i][1][k]=dDm[i][1][k]+sgn*2*V[x[j]][k]*(V[x[y0[j]]][y0[i]]*One_array[x[y1[j]]]-V[x[y1[j]]][y0[i]]*One_array[x[y0[j]]]);
                    dDm[i][2][k]=dDm[i][2][k]+sgn*2*V[x[j]][k]*(V[x[y0[j]]][y0[i]]*V[x[y1[j]]][y1[i]]-V[x[y0[j]]][y1[i]]*V[x[y1[j]]][y0[i]]);
                }
                sgn=-sgn;
            }
        }
    }
}

template<unsigned int dim, typename T>
void V_tmp(int nb_vert,int v_c[4],Point<dim,T> (& V)[4],int One_array[4],Point<dim,T> (& Vtmp)[4])
{
    //make a temporary copy of the tetrahedron coordinates V, with some alterations
    //in the indices indicated by v_c[1:nb_vert-1], we replace the point by only its difference with the point in v_c[0]
    //We also put one or zero in One_array, 1 if the row is untouched, 0 if there was a substraction
    //we need a memory allocation so dont forget to free Vtmp out of the function!
    
    int i,to_modif,j;
    int ind=0;
    for(i=0;i<4;i++)
    {
        to_modif=0;
        if(ind<nb_vert)
        {
            if(v_c[ind]==i)
            {
                //vertex to modify (IF ind>=1!!)
                if(ind>=1)
                {
                    to_modif=1;
                }
                ind++;
            }
        }
        if(to_modif==1)//this is the case when we need to make the difference with v_c[0]
        {
            One_array[i]=0;
            for(j=0;j<3;j++)
            {
                //we only keep the offset part
                Vtmp[i][j]=V[i][j]-V[v_c[0]][j];
                //printf("%e ",Vtmp[i].coo[j]);
            }
        }
        else
        {
            One_array[i]=1;
            //else we copy-paste the vector V
            for(j=0;j<3;j++)
            {
                Vtmp[i][j]=V[i][j];
                //printf("%e ",Vtmp[i].coo[j]);
            }
        }
        //printf("\n");
    }
    //printf("\n");
}


template<unsigned int dim, typename T>
void compute_a(T (& a)[4],int v_c,Point<dim,T> (& V)[4],int (& One_array)[4])
{
    //compute the minor determinants of "a" knowing the row v_c that contains coordinates that vary, and the coordinates of the tetrahedron V
    //Also we introduce a row One_array[4] of 1 or 0 because we may have done operations (linear combinations) on the rows because of offsets
    int ro,co;
    T m[3][3];//matrix for minor determinants
    Point<dim,int> x=index_perm<dim>(v_c); //containing the index of the other vertice than v_c
    // v_c= 0 -> x=(1,2,3)
    //         1 -> x=(0,2,3)
    //         2 -> x=(0,1,3)
    //         3 -> x=(0,1,2)
    
    //test of consistency with the determinant of "a" entirely
    /*
     struct MatDoub m2=Mat_create(4, 4);
     for(ro=0;ro<4;ro++)
     {
     for(co=0;co<3;co++)
     {
     m2.M[ro][co]=V[ro].coo[co];
     }
     m2.M[ro][3]=1;
     }
     affiche_mat(&m2);
     det(&m2);
     free_mat(&m2);
     */
    
    int i;
    int y0[3]; y0[0]=1; y0[1]=0; y0[2]=0; //this is to eliminate certain columns to extract the minor
    int y1[3]; y1[0]=2; y1[1]=2; y1[2]=1;
    for(i=0;i<3;i++)
    {
        //constructing a[i] for i=0,1,2
        //fill the matrix
        for(ro=0;ro<3;ro++)
        {
            m[ro][0]=V[x[ro]][y0[i]];
            m[ro][1]=V[x[ro]][y1[i]];
            m[ro][2]=One_array[x[ro]];
        }
        //get the determinant
        //affiche_mat(&m);
        a[i]=det(m);
    }
    //then a[3] which is not including 1
    for(ro=0;ro<3;ro++)
    {
        for(co=0;co<3;co++)
        {
            m[ro][co]=V[x[ro]][co];
        }
    }
    //affiche_mat(&m);
    a[3]=det(m);
}

template<typename PolyMesh_type, typename Matrix_type>
void grad_voronoi_vertex_cellmove(PolyMesh_type & poly, int t,int c, Point<PolyMesh_type::dims,typename PolyMesh_type::stype> & circum, Matrix_type & M)
{
    constexpr int dims = PolyMesh_type::dims;
    typedef typename PolyMesh_type::stype T;

    auto & DelTet = poly.getDelTetraedron();
    int t_start = poly.template getVolumeProp<PolyMesh_type::DelTetStart>(c);

    //computes the gradient G of the voronoi vertex which is at the center of tetrahedron t
    //G is a 3x3 matrix, G[i][j] is the derivative of coordinate i of the vertex with respect to coordinate j of cell center c.
    //in circum we put the coordinates of the circumcenter (we have to compute it anyway)
    
    //first we compute some numerical quantities involved in the circumcenter that are required for the gradient
    //we need to compute the coordinates of every points
    Point<dims,T> V[4]; //stored in this array
    Point<dims,T> Vtmp[4];
    int i,j,c_ind;
    int v_c[4]={-1,-1,-1,-1};
    int nb_vert = 0;
    //if the cell does not belong in the tetrahedron this vector will remain equal to (-1,-1,-1,-1) and we can easily return (0,0,0)
    //if more than one vertex is connected to cell c, then it means that more than one vertex move at the same time and this requires specific computations
    for(i=0;i<4;i++)
    {
        c_ind=DelTet.template get<0>(t_start + t)[i];

        if(poly.getVolumeGID(c_ind) == poly.getVolumeGID(c))
        {
            v_c[nb_vert]=i;
            nb_vert++;
        }

        V[i] = poly.getVolumePos(c_ind);
    }
    
    if(nb_vert==0)//here gradient is zero because the cell center does not belong to the tetrahedron
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                M[i][j]=0;
            }
        }
        //we must exit the function in this case
        return;
    }

    //the circumcenter has coordinates (D[0]/(2*a),D[1]/(2*a),D[2]/(2*a))
    //where D[] and "a" are specific determinants of matrices formed from the coordinates of the tetrahedron
    //They are 4x4 determinants but we will develop them in 4 3x3 determinants because anyway this is what we need to do to compute the gradient
    // a is given by:
    // | V[0].coo[0] V[0].coo[1] V[0].coo[2] 1 |
    // | V[1].coo[0] ...... .....   ...... . 1 |
    // |  .................................... |
    // | V[3].coo[0] ......      V[3].coo[2] 1 |
    
    //we expand the determinant along the row v_c
    //we get a= sum_{i=0}^{i=2} (-1)^(v_c+i)*V[v_c].coo[i]*a[i] + (-1)^(xc+3)*a[3]
    //where a[i] is the minor determinant (3x3) of "a" taken at row v_c and column i
    double a[4];
    //Careful, if n_vert>1, then we must temporarily modify certain rows of a, the trick is to take all rows v_c[1:n_vert-1] and remove from them the row v_c[0]
    //because it does not change the value of the determinant, but it is practical because then the only "variable" part of the determinant is the row v_c[0] and we can use the same method to compute the gradient
    //by doing this operation some 1 become zeroes in the last column, this must also be sent to the function compute a
    int One_array[4]; //array storing ones or zeroes for computing "a"
    if(nb_vert==1)
    {
        for(i=0;i<4;i++)
        {
            One_array[i]=1;
        }
        compute_a(a, v_c[0], V,One_array);
    }
    else
    {
        V_tmp(nb_vert,v_c,V,One_array,Vtmp);
        compute_a(a, v_c[0], Vtmp,One_array);
    }
    //compute the value of "a" using a[]
    double a_f=0;
    int sgn;
    if(v_c[0]%2 == 0){sgn=1;}else{sgn=-1;}
    for(i=0;i<3;i++)
    {
        a_f=a_f+sgn*V[v_c[0]][i]*a[i];
        sgn=-sgn;
    }
    a_f=a_f+sgn*a[3];
    //the D[] array is given by
    // D[i] for i=0,1,2 is:
    // D[i]= (-1)^(i+2) * | r2[0] V[0].coo[y0[i]] V[0].coo[y1[i]] 1 |
    //                    | r2[1].................................  |
    //                    | r2[2].................................. |
    //                    | r2[3].................V[3].coo[y1[i]] 1 |
    //with r2[j] is the norm squared of V[j]
    //y0[]=[1,0,0]
    //y1[]=[2,2,1]
    //again we develop along the row v_c
    //D[i]=(-1)^(i+2)*( (-1)^(v_c)*r2[v_c]*D[i][0] + (-1)^(v_c+1)*V[v_c].coo[y0[i]]*D[i][1] + (-1)^(v_c+2)*V[v_c].coo[y1[i]]*D[i][2] +(-1)^(v_c+3)*D[i][3]
    //if nb_vert>1 then we need to modify certain rows of D by substracting the row v_c from them
    //this has already been done by Vtmp and One_array but we also need to do it for r2
    double r2[4];
    //first we compute the norm squared of each vector
    for(i=0;i<4;i++)
    {
        r2[i] = norm2(V[i]);
    }
    if(nb_vert>=2)//then if needed we perform the substractions
    {
        for(i=0;i<4;i++)
        {
            if(One_array[i]==0)//this row has a substraction
            {
                r2[i]=r2[i]-r2[v_c[0]];
            }
        }
    }
    double D[3][4]; //this is to store the minor determinants of D[i]
    double dDm[3][3][3];//And this to store the derivative of D[i][1:3][k] (with respect to dim k) which might be non zero!
    if(nb_vert==1)
    {
        compute_D(D,dDm, v_c[0], V,r2,One_array);
    }
    else
    {
        //we use Vtmp in this case because some rows are modified
        compute_D(D,dDm, v_c[0], Vtmp,r2,One_array);
    }
    double D_f[3]; //the value of D[i]
    double tmp;
    int sgni=1; //starting (-1)^(i+2) is 1
    int sgnvc;
    int y0[3]; y0[0]=1; y0[1]=0; y0[2]=0;
    int y1[3]; y1[0]=2; y1[1]=2; y1[2]=1;
    for(i=0;i<3;i++)
    {
        if(v_c[0]%2==0){sgnvc=1;}else{sgnvc=-1;}//(-1)^(v_c)
        D_f[i]=0;
        tmp=sgnvc*r2[v_c[0]]*D[i][0];
        sgnvc=-sgnvc;
        tmp=tmp+sgnvc*V[v_c[0]][y0[i]]*D[i][1];
        sgnvc=-sgnvc;
        tmp=tmp+sgnvc*V[v_c[0]][y1[i]]*D[i][2];
        sgnvc=-sgnvc;
        tmp=tmp+sgnvc*D[i][3];
        
        D_f[i]=sgni*tmp;
        sgni=-sgni;
    }
    
    //lets return the circumcenter
    for(i=0;i<3;i++)
    {
        circum[i]=D_f[i]/(2*a_f);
    }
    
    //now the gradient is just going to be given by the following formula (derivative of circumcenter cen with respect to coordinate k of cell c)
    //d(cen.coo[i])/d(k) = (1/(2*a))*d(D[i])/d(k)-(D[i]/(2*a^2))*d(a)/d(k)
    //start first by d(a)/d(k), d(k) actually means d(V[v_c].coo[k])
    //d(a)/d(k) = (-1)^(v_c+k)*a[k]
    Point<dims,T> da;
    int k;
    for(k=0;k<3;k++)
    {
        if((v_c[0]+k)%2 == 0){sgn=1;}else{sgn=-1;}
        da[k]=sgn*a[k];
    }
    //then d(D[i])/d(k) = (-1)^(i+2)*( (-1)^(v_c)*2*V[v_c].coo[k]*D[i][0] + if(y0[i]==k){ (-1)^(v_c+1)*D[i][1] } + if(y1[i]==k){ (-1)^(v_c+2)*D[i][2]} )
    //if nb_vert>1 there will be additional terms because the D[i][1:3] can have a non zero derivative in this case
    //what we do is that we expand the 3x3 determinant along the first column (r2) because this first column contains the only variable terms (susceptible to contribute to the derivative). If a row x has been modified the column contains (at x)  sum_i(2*V[v_c].coo[i]*Off[x].coo[i]+(Off[x].coo[i])^2)
    //the V[v_c] terms show that there will be a non-zero derivative
    double dD[3][3]; //dD[i][k] = d(D[i])/d(k)
    for(k=0;k<3;k++)
    {
        sgni=1; //first even i
        for(i=0;i<3;i++)
        {
            if(v_c[0]%2 == 0){sgn=1;}else{sgn=-1;} //sign of vc
            tmp=sgn*2*V[v_c[0]][k]*D[i][0];
            sgn=-sgn;
            if(y0[i]==k)
            {
                tmp=tmp+sgn*D[i][1];
            }
            tmp=tmp+sgn*V[v_c[0]][y0[i]]*dDm[i][0][k]; //derivative of D[i][1], no need that y0[i]=k because if there is a modified column then there is always a term to put here
            sgn=-sgn;
            if(y1[i]==k)
            {
                tmp=tmp+sgn*D[i][2];
            }
            tmp=tmp+sgn*V[v_c[0]][y1[i]]*dDm[i][1][k];
            sgn=-sgn;
            //also potential derivative of D[i][3]
            tmp=tmp+sgn*dDm[i][2][k];
            dD[i][k]=sgni*tmp;
            sgni=-sgni;
        }
    }
    
    //finally compute the gradient G[i][k] the derivative of coo i of vertex with respect to coo k of cell center
    for(k=0;k<3;k++)
    {
        for(i=0;i<3;i++)
        {
            M[i][k]=dD[i][k]/(2*a_f)-D_f[i]*da[k]/(2*a_f*a_f);
        }
    }
}


/*! \brief Calculate the Gv quantity on each vertex
 *
 * Fills the allocated Gv with Matrices so that Gv[i]->M[j][k] corresponds to the derivative of coordinate 
 * j of vertex i of the cell c_v (v as "volume"), when the cell c_m (m as "moving") is moved along direction k
 * 
 */
template<typename PolyMesh_type, typename Gv_type>
void grad_voronoi_cell(PolyMesh_type & Poly, unsigned int c, Gv_type & Gv)
{
    Gv.clear();

    constexpr int DelTetStart = PolyMesh_type::DelTetStart;
    constexpr int DelTetNumber = PolyMesh_type::VolumeNumberOfDelTet;

    auto & DelTet = Poly.getDelTetraedron();

    auto it = Poly.getVolumeDomainIterator();

    auto Nid = Poly.template getVolumeProp<DelTetNumber>(c);

    Gv.resize(Nid);

    for (int t = 0 ; t < Nid ; t++)
    {
        Point<PolyMesh_type::dims,typename PolyMesh_type::stype> circum;
        grad_voronoi_vertex_cellmove(Poly,t,c,circum,Gv.template get<0>(t));
    }
}


template<unsigned int dim, typename T>
void derivate_V(T (& dV)[4][dim],T (& A)[4][4], bool (& Mask)[4][dim])
{
    T Ac44[4][4];
    T Ac43[3][3];
    //computes the derivative of the volume of a tetrahedron with respect to a movement of each of its points in each direction
    //The tetrahedron is given in the shape of a matrix A(4x4) like this:
    //A = | x0 y0 z0 1 |
    //    | x1 y1 z1 1 |
    //    | x2 y2 z2 1 |
    //    | x3 y3 z3 1 |
    //The determinant of this matrix is the 6*volume*sgn (sgn=+1 or -1 depending on row order)
    //we return the derivative with respect to xi,yi,zi in a matrix dV (4x3)
    //dV = sgn/6 * | dA/dx0 dA/dy0 dA/dz0 |
    //             | dA/dx1 dA/dy1 dA/dz1 |
    //             | dA/dx2 dA/dy2 dA/dz2 |
    //             | dA/dx3 dA/dy3 dA/dz3 |
    
    //We now use a Mask to speed up the computation. wo->Mask is a 4x3 matrix containing 0 or 1
    //If Mask->M[i][j] is zero it means that we don't need to compute dV->M[i][j] because the coordinate [i][j] does not change, so it will not be used. We simply put dV->M[i][j] to zero in this case. This leads to less determinant computations!
    
    //the matrices wo->Ac44 (4x4) and wo->Ac43 (4x3) are just temporary matrices used for internal computations, we put them as input in order to avoid making an unecessary large amount of memory allocations
    
    //first copy of A needed to know the sign of the determinant
    int r,c;
    for(r=0;r<4;r++)
    {
        for(c=0;c<4;c++)
        {
            Ac44[r][c]=A[r][c];
        }
    }
    //compute determinant to know the sign
    double detA=det(Ac44);
    //printf("det: %e\n",detA/6);
    int sgn;
    if(detA<0)
    {
        sgn=-1;
    }
    else
    {
        sgn=1;
    }
    //now the derivative of A (row r, column c) is given by (-1)^(r+c)*M[r][c]
    //where M is the determinant of the minor matrix M (row r and column c removed from A)
    //we need an array to generate the indices to put in the minor
    Point<dim,int> y[4];
    y[0][0]=1; y[0][1]=2; y[0][2]=3;
    y[1][0]=0; y[1][1]=2; y[1][2]=3;
    y[2][0]=0; y[2][1]=1; y[2][2]=3;
    y[3][0]=0; y[3][1]=1; y[3][2]=2;
    
    //Ac43 will store the minor matrix
    int r2,c2,sgn2;
    for(r=0;r<4;r++)
    {
        for(c=0;c<3;c++)
        {
            //we do something only if it's needed according to the mask
            if(Mask[r][c] != 0)
            {
                //creating the minor matrix
                for(r2=0;r2<3;r2++)
                {
                    for(c2=0;c2<3;c2++)
                    {
                        Ac43[r2][c2]=A[y[r][r2]][y[c][c2]];
                    }
                }
                //computing determinant of it
                //detA=det(Ac43);
                detA=det(Ac43); //faster determinant
                //sign to put
                if((r+c)%2 == 0){sgn2=1;}else{sgn2=-1;}
                dV[r][c]=sgn*sgn2*detA/6.0;
            }
            else //else we put 0
            {
                dV[r][c]=0;
            }
        }
    }
}

template<typename PolyMesh_type, typename Gv_type, typename Gb_type>
Point<PolyMesh_type::dims,typename PolyMesh_type::stype> derivate_V_cell(PolyMesh_type & poly,
                            int c_v, 
                            int c_m, 
                            Gv_type & Gv, 
                            Gb_type & Gb, 
                            bool (& Mask)[4][PolyMesh_type::dims])
{
    constexpr int dim = PolyMesh_type::dims;
    typedef typename PolyMesh_type::stype T;

    constexpr int DelTetStart = PolyMesh_type::DelTetStart;
    constexpr int DelTetNumber = PolyMesh_type::VolumeNumberOfDelTet;

    constexpr int DelTelVertId = PolyMesh_type::DelTelVertId;
    auto & DelTet = poly.getDelTetraedron();

    //Returns a vector V so that V[i] is the derivative of the volume of the cell c_v with respect to the coordinate i of the center of the cell c_m
    //wo->Ac44 (4x4) and wo->Ac43 (4x3) are just temporary matrices allocated outside of the function to minimize the number of malloc calls, wo->Mask is used by derivate_V
    //Gv contains the derivative of each voronoi vertex of cell c_v with respect to cell c_m (and should actually correspond to either wo->Gv or wo->Gv2)
    //Gb contains the derivative of the center of each face of the cell c_v with respect to c_m (and should actually correspond to either wo->Gb or wo->Gb2)
    
    //sum the variation of the volume of each tetrahedron (center of cell, center of face, vertex 0, vertex1) that make up the bulk of the cell
    Point<dim,T> g_vol({0,0,0}); //init gradient
    Point<dim,T> g_vol_f; //temporary gradient for each face
    int igv2,igv3;
    int f = 0;
    T A[4][4];
    T dV[4][3];
    Point<dim,T> o({0,0,0});
    Point<dim,T> v;

//    for(f=0;f<s->Cell[c_v].nf;f++) //loop on faces
//    {

    poly.ForAllVolumeSurfaces(c_v,[&](int f_ind, int conn){

        //init face gradient to zero
        g_vol_f = 0.0;
        //f_ind=s->Cell[c_v].F[f];
        //start constructing the tetrahedron in matrix form
        //At the same time we create the mask used for the derivate_V function, kicking out terms that don't need to be computed
        //cell center
        for(int x=0;x<3;x++)
        {
            A[0][x]=poly.getVolumePos(c_v)[x];
        }
        //The cell c_m is moving so if c_m != c_v we can skip these terms
        if(c_m != c_v)
        {
            for(int x=0;x<3;x++)
            {
                Mask[0][x]=0;
            }
        }
        else
        {
            for(int x=0;x<3;x++)
            {
                Mask[0][x]=1;
            }
        }

        for(int x=0;x<3;x++)
        {
            A[1][x]=poly.getSurfacePos(f_ind)[x];
            //we can skip dV.M[1][x] depending on Gb[f].M[x][:]
            if(Gb.template get<0>(f)[x][0] == 0 && Gb.template get<0>(f)[x][1] == 0 && Gb.template get<0>(f)[x][2] == 0)
            {
                Mask[1][x]=0;
            }
            else{Mask[1][x]=1;}
        }
        A[0][3]=1; //dont forget the column of one
        A[1][3]=1;
        A[2][3]=1;
        A[3][3]=1;

        poly.ForAllSurfaceEdges(f_ind,[&](int e){


            for(int x=0;x<3;x++)
            {
                A[2][x]=poly.getEdgeVertexPos(e,0)[x];
            }
            //let's look in Gv to which index it corresponds
            auto start_t = poly.template getVolumeProp<DelTetStart>(c_v);
            auto nt = poly.template getVolumeProp<DelTetNumber>(c_v);
            for(igv2 = 0 ;igv2 < nt  ;igv2++)
            {
                if( DelTet.template get<DelTelVertId>(start_t + igv2) == poly.getEdgeVertexID(e,0) )
                {
                    break;
                }
            }

            //now we skip dV.M[2][x] depending on Gv[igv2].M[x][:]
            for(int x=0;x<3;x++)
            {
                if(Gv.template get<0>(igv2)[x][0]==0 && 
                   Gv.template get<0>(igv2)[x][1]==0 && 
                   Gv.template get<0>(igv2)[x][2]==0)
                {
                    Mask[2][x]=0;
                }
                else{Mask[2][x]=1;}
            }
            for(int x=0;x<3;x++)
            {
                A[3][x]=poly.getEdgeVertexPos(e,1)[x];
            }
            //let's look in Gv to which index it corresponds
            for(igv3 = 0 ;igv3 < nt ; igv3++)
            {
                if( DelTet.template get<DelTelVertId>(start_t + igv3) == poly.getEdgeVertexID(e,1) )
                {
                    break;
                }
            }
            //now we skip dV.M[3][x] depending on Gv[igv3].M[x][:]
            for(int x=0;x<3;x++)
            {
                if(Gv.template get<0>(igv3)[x][0]==0 && Gv.template get<0>(igv3)[x][1]==0 && Gv.template get<0>(igv3)[x][2]==0)
                {
                    Mask[3][x]=0;
                }
                else{Mask[3][x]=1;}
            }
            //computes the gradient of the tetrahedron.
            //using Ac44 and Ac43 as temporary matrices(avoid memory allocation)
            derivate_V(dV, A, Mask);
            //then from this we compute the volume change by multiplying the matrix dV with all the variations of the points of the tetrahedron
            //the point 0 is the cell center, it only moves along the dimension k IF c_v==c_m, the gradient is 1 so we just put dV[0][k]
            if(c_v==c_m)
            {
                for(int k=0;k<3;k++)
                {
                    g_vol_f[k]=g_vol_f[k]+dV[0][k];
                }
            }
            //the point 1 is the face_center which moves by the gradient Gb[f]
            for(int x=0;x<3;x++)
            {
                for(int k=0;k<3;k++)
                {
                    g_vol_f[k]=g_vol_f[k]+dV[1][x]*(Gb.template get<0>(f)[x][k]);
                }
            }
            //point 2
            for(int x=0;x<3;x++)
            {
                for(int k=0;k<3;k++)
                {
                    g_vol_f[k]=g_vol_f[k]+dV[2][x]*(Gv.template get<0>(igv2)[x][k]);
                }
            }
            //point 3
            for(int x=0;x<3;x++)
            {
                for(int k=0;k<3;k++)
                {
                    g_vol_f[k]=g_vol_f[k]+dV[3][x]*(Gv.template get<0>(igv3)[x][k]);
                }
            }
            
        });
        //we dont forget to multiply by self_ad
        //but careful, only do this to the gradient of the current face
        for(int k=0;k<3;k++)
        {
            g_vol_f[k]=g_vol_f[k];
            //then add it to the global gradient
            g_vol[k]=g_vol[k]+g_vol_f[k];
        }
        ++f;
    });

    return g_vol;
}

template<typename PolyMesh_type, typename Gv_type, typename Gb_type>
void grad_barycenter(PolyMesh_type & poly, int c, Gv_type & Gv, Gb_type & Gb)
{
    constexpr int DelTetStart = PolyMesh_type::DelTetStart;
    constexpr int DelTetNumber = PolyMesh_type::VolumeNumberOfDelTet;
    constexpr int DelTelVertId = PolyMesh_type::DelTelVertId;

    constexpr int VolumeNumberOfSurfaces = PolyMesh_type::VolumeNumberOfSurfaces;

    auto & DelTet = poly.getDelTetraedron();
    auto & Volumes = poly.getVolumes();

    //Fills Gb with the gradient on the barycenters of all the faces of cell c, provided an array of gradients of individual vertices Gv (array obtained from function grad_voronoi_cell). It is important that the cell c is the cell that has been used in grad_voronoi_cell, otherwise there will be memory errors
    //Gb is ordered in the same way as the faces of the cell are ordered
    // Gb[i]->M[j][k] = derivative of coo j of barycenter i with respect to coo k of cell c

    int f = 0;
    //loop on faces
    int f_ind;
    int e,k,x,ind_v,igv0,igv1;

    Gb.resize(Volumes.template getProp<VolumeNumberOfSurfaces>(c));

    poly.ForAllVolumeSurfaces(c,[&](int f_ind, int conn){

        //init gradient to zero
        for(k=0;k<3;k++)
        {
            for(e=0;e<3;e++)
            {
                Gb.template get<0>(f)[e][k]=0;
            }
        }

        int ne = 0;

        //we loop on edges which means that we count each vertex twice
        poly.ForAllSurfaceEdges(f_ind,[&](int e){

            auto ind_v = poly.getEdgeVertexID(e,0);

            auto tet_start = poly.template getVolumeProp<DelTetStart>(c);
            auto nt = poly.template getVolumeProp<DelTetNumber>(c);

            for(igv0 = tet_start; igv0 < nt + tet_start ; igv0++)
            {
                if(DelTet.template get<DelTelVertId>(igv0)==ind_v)
                {
                    break;
                }
            }

            //vertex 1
            ind_v = poly.getEdgeVertexID(e,1);
            for(igv1 = tet_start; igv1 < nt + tet_start ; igv1++)
            {
                if(DelTet.template get<DelTelVertId>(igv1) == ind_v)
                {
                    break;
                }
            }

            for(x=0;x<3;x++)//coo loop
            {
                for(k=0;k<3;k++)
                {
                    Gb.template get<0>(f)[x][k] = Gb.template get<0>(f)[x][k]+Gv.template get<0>(igv0 - tet_start)[x][k]; //here memory error can happen if wronf cell has been used
                    Gb.template get<0>(f)[x][k] = Gb.template get<0>(f)[x][k]+Gv.template get<0>(igv1 - tet_start)[x][k];
                }
            }
            ne++;
        });

        //each vertice has been counted twice so there are 2*ne weights
        for(x=0;x<3;x++)
        {
            for(k=0;k<3;k++)
            {
                Gb.template get<0>(f)[x][k] = Gb.template get<0>(f)[x][k]/(2*ne);
            }
        }
        f++;
    });
}


#endif