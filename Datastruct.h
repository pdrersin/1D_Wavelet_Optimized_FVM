#ifndef Datastruct_h_inluded
#define Datastruct_h_inluded
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file Datastruct.h
/// \brief Header for data structure related functions
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void check_leaves(unordered_map<int,Cell*> &Cellvect);
/// \brief Marks leaves
/// \brief Just sweeps through entire unordered map
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void check_leaves(unordered_map<int,Cell*> &Cellvect)
{
    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        int iref = (*citer_one).second->i;
        (*citer_one).second->leaf = 0;

        if (!find_cell_exist(2*iref+1,Cellvect) && !find_cell_exist(2*iref,Cellvect))//No Children
        {
            (*citer_one).second->leaf = 1;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn bool find_cell_exist(int i, unordered_map<int,Cell*> &Cellvect);
/// \brief Logical function to find if cell exists
/// \brief Note that the first element of vector must be i=1
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
bool find_cell_exist(int i, unordered_map<int,Cell*> &Cellvect)
{
    //Note that first element of Cellvect has to be the root.
    int child_tracker[max_level];
    int level;

    if (i==1)//Root always remains in vector
    {
       return true;
    }


    for (level=0;level<max_level;level++)//Increment after body of loop executes
    {
        if (i%2==0)//Left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//Right Child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
    }

    if (i!=1)
    {
        return false;
    }

    Cell* loc = Cellvect[1];

    do{
        if (child_tracker[level]==1)// left child
        {

            if (loc->left_child == nullptr)//Cell does not exist
            {
                return false;
            }
            else
            {
                loc = loc->left_child;
            }
        }
        else if (child_tracker[level]==2)//Right child
        {
            if (loc->right_child == nullptr)//Cell does not exist
            {
                return false;
            }
            else
            {
                loc = loc->right_child;
            }
        }
        level--;

    }while(level>=0);

    return true;//If it has made it here that means the cell exists

}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn Cell* find_cell(int, int, unordered_map<int,Cell*> &Cellvect);
/// \brief Returns cell* to particular i value
/// \brief ALWAYS use once cell is confirmed to exist
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
Cell* find_cell(int i, unordered_map<int,Cell*> &Cellvect)//Use after verifying that cell exists
{
    //Note that first element of Cellvect has to be the root.
    int child_tracker[max_level];
    int level;

    if (i==1)//Root always remains in vector
    {
        Cell* loc = Cellvect[1];
        return loc;
    }


    for (level=0;level<max_level;level++)//Increment after body of loop executes
    {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)// right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
    }

    Cell* loc=Cellvect[1];

    do{

        if (child_tracker[level]==1)//left child
        {
                loc = loc->left_child;
        }
        else if (child_tracker[level]==2)//right child
        {
                loc = loc->right_child;
        }
        level--;
    }while(level>=0);

    return loc;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void link_level_neighbors(unordered_map<int,Cell*> &Cellvect)
/// \brief Links to neighboring level if cell exists
/// \brief Required for interpolations
/// \brief Assumes periodicity(1) or Open(2)
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void link_level_neighbors(unordered_map<int,Cell*> &Cellvect)
{

    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        int iright = (*citer_one).second->i + 1;
        int ileft = (*citer_one).second->i - 1;

        //Linking right cells within domain
        if (find_cell_exist(iright,Cellvect))
        {
            (*citer_one).second->right_level = find_cell(iright,Cellvect);
        }

        //Linking left cells within domain
        if (find_cell_exist(ileft,Cellvect))
        {
            (*citer_one).second->left_level = find_cell(ileft,Cellvect);
        }

        //Linking right direction boundary condition
        if (iright == pow(2,(*citer_one).second->level+1))
        {
            if (bc_type==1)//Periodic
            {
                int per_iright = pow(2,(*citer_one).second->level);
                if (find_cell_exist(per_iright,Cellvect))
                {
                    (*citer_one).second->right_level = find_cell(per_iright,Cellvect);
                }
            }
            else//Open
            {
                (*citer_one).second->right_level = (*citer_one).second;
            }

        }

        //Linking left direction boundary condition
        if (ileft == pow(2,(*citer_one).second->level)-1)
        {
            if (bc_type==1)
            {
                int per_ileft = pow(2,(*citer_one).second->level+1)-1;
                if (find_cell_exist(per_ileft,Cellvect))
                {
                    (*citer_one).second->left_level = find_cell(per_ileft,Cellvect);
                }
            }
            else
            {
                (*citer_one).second->left_level = (*citer_one).second;
            }

        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void keep_cell(int i, unordered_map<int,Cell*> &Cellvect);
/// \brief Either changes existent cell keep_flag to 1 or makes new cell
/// \brief Note that this does not decode
/// \brief Must set new_cell=1 in this function
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void keep_cell(int i, int level, unordered_map<int,Cell*> &Cellvect)
{


    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }


    //Periodicity is taken care of now

    if (find_cell_exist(i,Cellvect))
    {
        Cell* loc = find_cell(i,Cellvect);
        loc->keep_flag = 1;
    }
    else
    {
        int child_tracker[max_level];
        int level;

        for (level=0;level<max_level;level++)//Charting route of cell from root
        {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//Right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
        }

        Cell* loc = Cellvect[1];
        int insert_level = 1;

        do{

            if (child_tracker[level]==1)//Lower left child
            {

                if (loc->left_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i;
                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;;

                    loc->left_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i]=nc;//Insertion into unique key

                    loc = nc;

                    i = 2*i;

                }
                else
                {
                    loc = loc->left_child;
                    i = 2*i;
                }


            }
            else if (child_tracker[level]==2)//Right child
            {
                if (loc->right_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;
                    nc->i = 2*i+1;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;

                    loc->right_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i+1]=nc;//Insertion into unique key

                    loc = nc;
                    i = 2*i+1;

                }
                else
                {
                    loc = loc->right_child;
                    i = 2*i+1;
                }
            }
            level--;
            insert_level++;
        }while(level>=0);
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn unordered_map<Cell*> level_vector(int level, unordered_map<int,Cell*> &Cellvect);
/// \brief Returns a vector of memory addresses to elements of the same level
/// \brief Useful for level operations which are incremented or decremented
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
unordered_map<int,Cell*> level_vector(int level, unordered_map<int,Cell*> &Cellvect)
{
    unordered_map<int,Cell*>levelvector;

    if (level==0)
    {
        Cell* loc = Cellvect[1];
        levelvector[1] = loc;
        return levelvector;
    }

    int low_limit = int(pow(2,level));
    int high_limit = int(pow(2,level+1))-1;

    for (int i=low_limit;i<=high_limit;i++)
    {
            if (find_cell_exist(i,Cellvect))
            {
            Cell* loc = find_cell(i,Cellvect);
            levelvector[i] = loc;
            }
    }
    return levelvector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn unordered_map<Cell*> leaf_vector(unordered_map<int,Cell*> &Cellvect);
/// \brief Returns a vector of memory addresses to leaves
/// \brief Useful for final stage time advancement
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
unordered_map<int,Cell*> leaf_vector(unordered_map<int,Cell*> &Cellvect)
{
    unordered_map<int,Cell*>leafvector;
    for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        if((*citer_one).second->leaf==1)
        {
            int i = (*citer_one).second->i;
            Cell* loc = find_cell(i,Cellvect);
            leafvector[i] = loc;
        }
    }
    return leafvector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void keep_cell_virt(int i, unordered_map<int,Cell*> &Cellvect);
/// \brief Makes new virtual cell
/// \brief Note that this does not decode
/// \brief Must set new_cell=1 and virt=1 in this function
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void keep_cell_virt(int i, int level, unordered_map<int,Cell*> &Cellvect)
{

    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }

    //Periodicity is taken care of now

    if (find_cell_exist(i,Cellvect))
    {
        Cell* loc = find_cell(i,Cellvect);
        loc->keep_flag = 1;
    }
    else
    {
        int child_tracker[max_level];
        int level;

        for (level=0;level<max_level;level++)//Charting route of cell from root
        {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
        }

        Cell* loc = Cellvect[1];
        int insert_level = 1;
        int virt_tracker = 0;

        do{
            if (child_tracker[level]==1)//left child
            {

                if (loc->left_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;;

                    loc->left_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i] = nc;

                    loc = nc;
                    i = 2*i;
                    virt_tracker++;

                }
                else
                {
                    loc = loc->left_child;
                    i = 2*i;
                }
            }
            else if (child_tracker[level]==2)//Upper left child
            {
                if (loc->right_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i+1;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;;

                    loc->right_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;
                    Cellvect[2*i+1] = nc;

                    loc = nc;
                    i = 2*i+1;
                    virt_tracker++;
                }
                else
                {
                    loc = loc->right_child;
                    i = 2*i+1;
                }
            }

            level--;
            insert_level++;

        }while(level>=0);

        if (virt_tracker!=0)
        {
        find_cell(i,Cellvect)->virt = 1;
        find_cell(i,Cellvect)->leaf = 1;
        find_cell(i,Cellvect)->keep_flag = 1;
        }

        find_cell(i,Cellvect)->keep_flag = 1;

    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void delete_cell(int i, unordered_map<key,Cell*> &Cellvect);
/// \brief Severs pointers for parent and neighbors
/// \brief Does not sever pointers children as if cell has children it must be in graded structure
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void delete_cell(int i, unordered_map<int,Cell*> &Cellvect)
{
    Cell* loc = find_cell(i,Cellvect);
    //Severing pointer connection with parent

    if (i%2==0)//left child
    {
        loc->parent->left_child = nullptr;
    }
    else if (i%2!=0)//right child
    {
        loc->parent->right_child = nullptr;
    }

    //Delinking neighbor connections
    if (loc->left_level != nullptr)
    {
        loc->left_level->right_level = nullptr;
    }

    if (loc->right_level != nullptr)
    {
        loc->right_level->left_level = nullptr;
    }

    Cellvect.erase(i);

    delete(loc);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn void periodize(int i, int level);
/// \brief Periodizes a particular i value for Periodic BC
/// \brief Connects boundary node to itself for Open BC
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
int periodize(int i, int level)
{
    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }

    return i;
}

#endif // Datastruct_h_inluded
