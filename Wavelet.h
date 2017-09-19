#ifndef Wavelet_h_inluded
#define Wavelet_h_inluded

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file Wavelet.h
/// \brief Header for wavelet coefficient calculation
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void encode(unordered_map<int,Cell*> &);
/// \brief This is a control panel for calculation of wavelet coefficients
/// \brief Calls two functions for projection and detail calculation respectively
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void encode(unordered_map<int,Cell*> &Cellvect)
{
    projection(Cellvect);

    details_calculation(Cellvect);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void projection(unordered_map<int,Cell*> &Cellvect);
/// \brief Projection from leaves towards root
/// \brief Calculates average cell value of variable at parent
/// \brief Uses equation 2.28 in Tenaud tutorial
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void projection(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level-1;

    while(level>=0)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->leaf==0)
            {
                double current_cell_area = (*citer_one).second->cell_length;
                double child_cell_area = (*citer_one).second->left_child->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double left_q = (*citer_one).second->left_child->q[j];
                    double right_q = (*citer_one).second->right_child->q[j];

                    double local_q = child_cell_area/current_cell_area*(left_q+right_q);
                    (*citer_one).second->q[j] = local_q;
                }
            }
        }
        levelvector.clear();
        level--;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void details_calculation(unordered_map<int,Cell*> &Cellvect);
/// \brief Uses projection to reinterpolate at finer levels
/// \brief Stores wavelet coefficients as difference of true vals and interpolated ones
/// \brief Uses step 6 in Algorithm 1 Tenaud tutorial
/// \brief Note that there are no details calculated for level 0 (or root since it has no parents)
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void details_calculation(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;

    while(level>=1)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {

                int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double parent_val = (*citer_one).second->parent->q[j];

                    double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                    double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                    double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                    double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                    }
                    else if (i%2!=0)//right child
                    {
                        interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));
                    }

                    (*citer_one).second->det[j] =  fabs((*citer_one).second->q[j] - interp_val);
                }


        }
        level--;
        levelvector.clear();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void preserve_graded_structure(unordered_map<int,Cell*> &Cellvect);
/// \brief This subroutine is used after marking cells for deletion or retention
/// \brief Any cell that must be kept has its parent stencil retained
/// \brief This follows steps 13-25: Algorithm 3 Tenaud tutorial
/// \brief Assumes tree is already graded at the end of previous refinement
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void preserve_graded_structure(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;

    while(level>1)//Parent of level 1 is root which is always included
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if (  (*citer_one).second->keep_flag == 1 )
            {
                int i = (*citer_one).second->i;
                int level_ref = (*citer_one).second->level;

                //Keeping parent
                keep_cell(i/2,level_ref-1,Cellvect);

                //Keeping left_level neighbor of parent
                keep_cell(i/2-1,level_ref-1,Cellvect);

                //Keeping right_level neighbor of parent
                keep_cell(i/2+1,level_ref-1,Cellvect);

                //Keeping left_left_level neighbor of parent
                keep_cell(i/2-2,level_ref-1,Cellvect);

                //Keeping right_right_level neighbor of parent
                keep_cell(i/2+2,level_ref-1,Cellvect);
            }
        }
        level--;
        levelvector.clear();//This makes sure Cell* pointers etc are not corrupted
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void hartens_predictive_thresholding(unordered_map<int,Cell*> &Cellvect);
/// \brief Used to determine if a cell must be kept in the tree structure or not
/// \brief In addition, adds virtual cells and decodes them to ensure stencil retention for fluxes
//////////////////////////////////////////////////////////////////////////////////////////////////
void hartens_predictive_thresholding(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level - 1;

    for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)//Step 3 in Tenaud tutorial Algorithm 5
    {
        (*citer_one).second->keep_flag = 0;
        (*citer_one).second->virt = 0;
        (*citer_one).second->leaf = 0;
        (*citer_one).second->new_cell = 0;
    }

    Cellvect[1]->keep_flag = 1;//Root must stay

//    double global_max[n_eq] = {0.0};
//    find_max_vals(Cellvect,global_max);

    while(level>0)
    {

        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        if (n_eq==1)
        {
            for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
            {
                //double local_det = (*citer_one).second->det[0]/global_max[0];
                double local_det = (*citer_one).second->det[0];
                double local_threshold = wc_threshold/pow(2,max_level-level);

                if (local_det>=local_threshold)
                {
                    int iref = (*citer_one).second->i;
                    int level_ref = (*citer_one).second->level;

                    //Keeping left_child
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    //Keeping left_child->left_level
                    keep_cell(2*iref-1,level_ref+1,Cellvect);
                   //Keeping right_child
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                    //Keeping right_child->right_level
                    keep_cell(2*iref+2,level_ref+1,Cellvect);

                    if (local_det>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    {
                            for (int i=4*iref-2;i<=4*iref+5;i++)
                            {
                                keep_cell(i,level_ref+2,Cellvect);
                            }
                    }
                }
            }
        }
        else
        {
            double max_detail = find_max_detail(levelvector);

            for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
            {

                double local_det = 0.0;

                for (int j=0;j<=n_eq-1;j++)
                {
                    local_det = local_det+fabs((*citer_one).second->det[j]);
                    //local_det = max(local_det,fabs((*citer_one).second->det[j]/global_max[j]));
                }

                double local_threshold = wc_threshold/pow(2,max_level-level);

                if (local_det/max_detail>local_threshold)
                //if (local_det>=local_threshold)
                {

                    int iref = (*citer_one).second->i;
                    int level_ref = (*citer_one).second->level;

                    //Keeping left_child
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    //Keeping left_child->left_level
                    keep_cell(2*iref-1,level_ref+1,Cellvect);
                   //Keeping right_child
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                    //Keeping right_child->right_level
                    keep_cell(2*iref+2,level_ref+1,Cellvect);

                    if (local_det/max_detail>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    //if (local_det>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    {
                            for (int i=4*iref-2;i<=4*iref+5;i++)
                            {
                                keep_cell(i,level_ref+2,Cellvect);
                            }
                    }
                }
            }
        }

        level--;
        levelvector.clear();
    }
    preserve_graded_structure(Cellvect);
    conservative_refinement(Cellvect);
    link_level_neighbors(Cellvect);
    decode_new_cells(Cellvect);
    delete_unnecessary_cells(Cellvect);
    check_leaves(Cellvect);
    add_virtual_cells(Cellvect);
    preserve_graded_structure(Cellvect);
    conservative_refinement(Cellvect);
    link_level_neighbors(Cellvect);
    decode_new_cells(Cellvect);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double find_max_detail(unordered_map<int,Cell*> &Cellvect);
/// \brief Finds the maximum value of the wavelet coefficient at a particular level of refinement
/// \brief Needed for Harten's thresholding
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
double find_max_detail(unordered_map<int,Cell*> &Cellvect)
{
    double max_detail=0.0;

    if (n_eq!=1)
    {
        for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            double local_det = 0.0;
            for (int j=0;j<=n_eq-1;j++)
            {
                local_det = local_det+fabs((*citer_one).second->det[j]);
            }

            max_detail = max(local_det,max_detail);
        }
    }
    else
    {
        for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            double l1sum = (*citer_one).second->det[0];
            max_detail = max(fabs(l1sum),max_detail);
        }
    }

    return max_detail;

}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void find_max_vals(unordered_map<int,Cell*> &Cellvect,double (&global_max)[n_eq]);
/// \brief Finds the maximum value of each conserved variable globally
/// \brief Needed for Harten's thresholding
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void find_max_vals(unordered_map<int,Cell*> &Cellvect,double (&global_max)[n_eq])
{
    for (int j=0;j<=n_eq-1;j++)
    {
        global_max[j] = 0.0;

        for (auto citer_one = Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            global_max[j] = max(global_max[j],fabs((*citer_one).second->q[j]));
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void conservative_refinement(unordered_map<int,Cell*> &Cellvect);
/// \brief Used to ensure that Hartens thresholding always has children of children etc
/// \brief If any cell has atleast one child - it should have all 2 children.
/// \brief Used within hartens_predictive_thresholding function
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void conservative_refinement(unordered_map<int,Cell*> &Cellvect)
{
    int level = 0;

    while(level<=max_level-1)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            int iref = (*citer_one).second->i;
            int level_ref = (*citer_one).second->level;

            if (find_cell_exist(2*iref,Cellvect))//left
            {
                if (find_cell(2*iref,Cellvect)->keep_flag == 1)
                {
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                }
            }

            if (find_cell_exist(2*iref+1,Cellvect))//right
            {
                if (find_cell(2*iref+1,Cellvect)->keep_flag == 1)
                {
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                }
            }
        }
    level++;
    levelvector.clear();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void add_virtual_cells(unordered_map<int,Cell*> &Cellvect);
/// \brief Refer page 27 Tenaud Tutorial
/// \brief All newly created cells must have virt=1
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void add_virtual_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = 0;

    while(level<=max_level)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
                if ((*citer_one).second->leaf==1)
                {
                int iref = (*citer_one).second->i;
                int lev = (*citer_one).second->level;
                int ileft = periodize(iref-1,lev);
                int iright = periodize(iref+1,lev);


                if (level<max_level)
                {
                    if (!find_cell_exist(iright,Cellvect))//Right Direction
                    {
                        keep_cell_virt(iref+1,level,Cellvect);
                        keep_cell_virt(iref+2,level,Cellvect);
                        keep_cell_virt(iref+3,level,Cellvect);
                        keep_cell_virt(iref+4,level,Cellvect);
                    }
                    else if (find_cell(iright,Cellvect)->leaf == 0)
                    {
                        keep_cell_virt(2*iref,level+1,Cellvect);
                        keep_cell_virt(2*iref+1,level+1,Cellvect);
                        keep_cell_virt(2*iref-1,level+1,Cellvect);
                        keep_cell_virt(2*iref-2,level+1,Cellvect);
                    }

                    if (!find_cell_exist(ileft,Cellvect))//Left Direction
                    {
                        keep_cell_virt(iref-1,level,Cellvect);
                        keep_cell_virt(iref-2,level,Cellvect);
                        keep_cell_virt(iref-3,level,Cellvect);
                        keep_cell_virt(iref-4,level,Cellvect);
                    }
                    else if (find_cell(ileft,Cellvect)->leaf == 0)
                    {
                        keep_cell_virt(2*iref,level+1,Cellvect);
                        keep_cell_virt(2*iref+1,level+1,Cellvect);
                        keep_cell_virt(2*iref+2,level+1,Cellvect);
                        keep_cell_virt(2*iref+3,level+1,Cellvect);
                    }
                }
                else
                {
                    if (!find_cell_exist(iright,Cellvect))//Right Direction
                    {
                        keep_cell_virt(iref+1,level,Cellvect);
                        keep_cell_virt(iref+2,level,Cellvect);
                        keep_cell_virt(iref+3,level,Cellvect);
                        keep_cell_virt(iref+4,level,Cellvect);
                    }

                    if (!find_cell_exist(ileft,Cellvect))//Left Direction
                    {
                        keep_cell_virt(iref-1,level,Cellvect);
                        keep_cell_virt(iref-2,level,Cellvect);
                        keep_cell_virt(iref-3,level,Cellvect);
                        keep_cell_virt(iref-4,level,Cellvect);
                    }
                }
                }
        }
        level++;
        levelvector.clear();
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void delete_unnecessary_cells(unordered_map<int,Cell*> &Cellvect);
/// \brief Deletes all cells with keep_flag zero. Simple.
/// \brief Must be used after all thresholding, virtual cells, conservative refinement etc
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void delete_unnecessary_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;
    while(level>0)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->keep_flag == 0)
            {
                int iref = (*citer_one).second->i;
                delete_cell(iref,Cellvect);
            }
        }
        level--;
        levelvector.clear();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void decode_new_cells(unordered_map<int,Cell*> &Cellvect);
/// \brief New cells created are given approximate field values
/// \brief Linkages must be complete (and tree must be graded prior to this)
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void decode_new_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = 1;

    while(level<=max_level)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
			continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->new_cell == 1)
            {
                int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

                //double parent_cell_length = (*citer_one).second->parent->cell_length;
                //double cell_length = (*citer_one).second->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        double parent_val = (*citer_one).second->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];
                    }
                    else if (i%2!=0)//right child
                    {
                        //interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));

                        double parent_val = (*citer_one).second->left_level->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->left_level->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->left_level->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->left_level->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->left_level->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];

                        interp_val = 2.0*parent_val - interp_val;

                    }

                    (*citer_one).second->q[j] = interp_val;
                }

                (*citer_one).second->new_cell = 0;
            }
        }
        level++;
        levelvector.clear();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void substage_decode(unordered_map<int,Cell*> &Cellvect);
/// \brief Virtual cells are updated for the substage solution of RK3
/// \brief Linkages must be complete (and tree must be graded prior to this)
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void substage_decode(unordered_map<int,Cell*> &Cellvect)
{

    encode(Cellvect);

    unordered_map<int,Cell*>leafvector;
    leafvector = leaf_vector(Cellvect);

    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        if ((*citer_one).second->virt == 1)
        {
            int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

//                double parent_cell_length = (*citer_one).second->parent->cell_length;
//                double cell_length = (*citer_one).second->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        double parent_val = (*citer_one).second->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];
                    }
                    else if (i%2!=0)//right child
                    {
                        //interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));

                        double parent_val = (*citer_one).second->left_level->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->left_level->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->left_level->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->left_level->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->left_level->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];

                        interp_val = 2.0*parent_val - interp_val;

                    }

                    (*citer_one).second->q[j] = interp_val;// + (*citer_one).second->det[j];
                }
        }
    }
}


#endif // Wavelet_h_inluded
