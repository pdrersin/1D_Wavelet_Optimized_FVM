#ifndef Init_Domain_h_inluded
#define Init_Domain_h_inluded
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file Init_Domain.h
/// \brief Header for initializing domain
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void create_domain(unordered_map<key,Cell*,hasher_key> &Cellvect)
/// \brief This creates the initial domain
/// \brief Also calls add_init_children(unordered_map<key,Cell*,hasher_key> &Cellvect);
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void create_domain(unordered_map<int,Cell*> &Cellvect)
{
    //Creating root of this particular domain
    Cell *nc = new Cell;
    nc->i = 1;
    nc->level = 0;//Root is kept level 0
    nc->cell_length = 1.0;
    nc->leaf = 0;

    nc->xx = 0.5;

    Cellvect[1]=nc;

    //Adding other instances
    add_init_children(Cellvect);

    //Checking leaves
    check_leaves(Cellvect);

    //Initial physical values to finest cells
    leaf_phys_vals(Cellvect);

    //Link Neighbors
    link_level_neighbors(Cellvect);

}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void add_init_children(unordered_map<int,Cell*> &Cellvect);
/// \brief This adds children to root
/// \brief Note that the first element of map must be i=1,j=1
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void add_init_children(unordered_map<int,Cell*> &Cellvect)
{
    int level = 1;

    int init_max_level=max_level;

    while(level<=init_max_level)
    {
        int upper_limit = int(pow(2,level+1))-1;
        int lower_limit = int(pow(2,level));

        for (int i=lower_limit;i<=upper_limit;i++)
        {
                Cell *nc = new Cell;
                nc->i = i;
                nc->level = level;


                nc->cell_length = 1.0/pow(2,level);


                nc->leaf = 0;
                nc->xx = 0.0;

                nc->parent = Cellvect[i/2];

                if (i%2==0)//left child
                {
                    nc->parent->left_child = nc;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;
                }
                else if (i%2!=0)//right child
                {
                    nc->parent->right_child = nc;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;
                }
                Cellvect[i]=nc;
        }
    level++;
    }

}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void leaf_phys_vals(unordered_map<int,Cell*> &Cellvect);
/// \brief Initial condition for physical quantities at leaves
/// \brief Note that only leaves have actual physical quantities
/// \brief This subroutine also specifies location of scalar phi
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void leaf_phys_vals(unordered_map<int,Cell*> &Cellvect)
{
    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        if ((*citer_one).second->leaf == 1)
        {
            double x = (*citer_one).second->xx;

            if (n_eq==1)
            {
                //double pi = 4.0*atan(1.0);
                //(*citer_one).second->q[0] = sin(pi*x*2.0);
                //(*citer_one).second->q[0] = exp(-50.0*(x-0.5)*(x-0.5));
                if (x>=0.5 && x<0.6)
                {
                    (*citer_one).second->q[0] = 1.0;
                }
                else
                {
                    (*citer_one).second->q[0] = 0.0;
                }
            }
            else if (n_eq==2)
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = height_1;
                    (*citer_one).second->q[1] = 0.0;
                }
                else
                {
                    (*citer_one).second->q[0] = height_2;
                    (*citer_one).second->q[1] = 0.0;
                }
            }
            else if (n_eq==3)//SOD Shock tube
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = 1.0;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 1.0/(gmm-1.0);
                }
                else
                {
                    (*citer_one).second->q[0] = 0.125;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.1/(gmm-1.0);
                }
            }
            else if (n_eq==7)
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = 1.0;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.0;
                    (*citer_one).second->q[3] = 0.0;
                    (*citer_one).second->q[4] = 1.0;
                    (*citer_one).second->q[5] = 0.0;
                    (*citer_one).second->q[6] = 1.0/(gmm-1.0) + 0.5*(bx*bx + 1.0);
                }
                else
                {
                    (*citer_one).second->q[0] = 0.125;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.0;
                    (*citer_one).second->q[3] = 0.0;
                    (*citer_one).second->q[4] = -1.0;
                    (*citer_one).second->q[5] = 0.0;
                    (*citer_one).second->q[6] = 0.1/(gmm-1.0) + 0.5*(bx*bx + 1.0);//Refer http://www.csun.edu/~jb715473/examples/mhd1d.htm
                }
            }
            else if (n_eq==14)
            {
                if (x<=0.5)
                {
                    double detF = 0.98*(1.0);
                    double r = rho0/detF;

                    (*citer_one).second->q[0] = r;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = r*0.5;
                    (*citer_one).second->q[3] = r;

                    double f[3][3] = {0.0};
                    f[0][0] = 0.98;
                    f[0][1] = 0.0;
                    f[0][2] = 0.0;
                    f[1][0] = 0.02;
                    f[1][1] = 1.0;
                    f[1][2] = 0.1;
                    f[2][0] = 0.0;
                    f[2][1] = 0.0;
                    f[2][2] = 1.0;

                    double s = 1.0e-3;

                    double energy = nle_en_calc(f,s);

                    (*citer_one).second->q[4] = r*(energy + 0.5*(0.5*0.5 + 1.0));
                    (*citer_one).second->q[5] = r*f[0][0];
                    (*citer_one).second->q[6] = r*f[0][1];
                    (*citer_one).second->q[7] = r*f[0][2];
                    (*citer_one).second->q[8] = r*f[1][0];
                    (*citer_one).second->q[9] = r*f[1][1];
                    (*citer_one).second->q[10] = r*f[1][2];
                    (*citer_one).second->q[11] = r*f[2][0];
                    (*citer_one).second->q[12] = r*f[2][1];
                    (*citer_one).second->q[13] = r*f[2][2];

                }
                else
                {
                    double detF = 1.0;
                    double r = rho0/detF;

                    (*citer_one).second->q[0] = r;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.0;
                    (*citer_one).second->q[3] = 0.0;

                    double f[3][3] = {0.0};

                    f[0][0] = 1.0;
                    f[0][1] = 0.0;
                    f[0][2] = 0.0;
                    f[1][0] = 0.0;
                    f[1][1] = 1.0;
                    f[1][2] = 0.1;
                    f[2][0] = 0.0;
                    f[2][1] = 0.0;
                    f[2][2] = 1.0;

                    double s = 0.0;

                    double energy = nle_en_calc(f,s);

                    (*citer_one).second->q[4] = r*(energy);
                    (*citer_one).second->q[5] = r*f[0][0];
                    (*citer_one).second->q[6] = r*f[0][1];
                    (*citer_one).second->q[7] = r*f[0][2];
                    (*citer_one).second->q[8] = r*f[1][0];
                    (*citer_one).second->q[9] = r*f[1][1];
                    (*citer_one).second->q[10] = r*f[1][2];
                    (*citer_one).second->q[11] = r*f[2][0];
                    (*citer_one).second->q[12] = r*f[2][1];
                    (*citer_one).second->q[13] = r*f[2][2];
                }
            }
//            else if (n_eq==4)//Elasto-Plastic Problem
//            {
//                if (x<=0.1)
//                {
//                    (*citer_one).second->q[0] = 2785.0;
//                    (*citer_one).second->q[1] = (*citer_one).second->q[0]*800.0;
//
//                    double p = 1.0e-6;
//
//                    double eta = 2785.0/ep_rho_0;
//                    double mg_f = (eta-1.0)*(eta-ep_gamma_0*(eta-1.0)/2.0)/((eta-ep_s*(eta-1.0))*(eta-ep_s*(eta-1.0)));
//
//                    double int_en = (p/ep_rho_0 - ep_a_0*ep_a_0*mg_f)/ep_gamma_0;
//
//                    double t_en = int_en + 0.5*(*citer_one).second->q[1]*(*citer_one).second->q[1]/(*citer_one).second->q[0];
//
//                    (*citer_one).second->q[2] = (*citer_one).second->q[0]*t_en;
//
//                    (*citer_one).second->q[3] = 2.0/3.0*ep_y_0;
//
//
//                }
//                else
//                {
//                    (*citer_one).second->q[0] = 2785.0;
//                    (*citer_one).second->q[1] = (*citer_one).second->q[0]*0.0;
//
//                    double p = 1.0e-6;
//
//                    double eta = 2785.0/ep_rho_0;
//                    double mg_f = (eta-1.0)*(eta-ep_gamma_0*(eta-1.0)/2.0)/((eta-ep_s*(eta-1.0))*(eta-ep_s*(eta-1.0)));
//
//                    double int_en = (p/ep_rho_0 - ep_a_0*ep_a_0*mg_f)/ep_gamma_0;
//
//                    double t_en = int_en + 0.5*(*citer_one).second->q[1]*(*citer_one).second->q[1]/(*citer_one).second->q[0];
//
//                    (*citer_one).second->q[2] = (*citer_one).second->q[0]*t_en;
//                    (*citer_one).second->q[3] = 2.0/3.0*ep_y_0;
//                }
//            }
        }
    }
}







#endif
