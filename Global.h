#ifndef Global_h_inluded
#define Global_h_inluded
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file Global.h
/// \brief Global header file containing Cell class
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

//Headers for global use
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <fstream>      //myfile.open()
#include <iomanip>      // std::setw
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This is the object for a finite volume cell.
///        Each cell has a unique i,j value to locate it in memory
/// \return
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
class Cell
{
public:


    double q[n_eq] = {0.0};
    int i;

    //MRA based storage
    double det[n_eq]={0.0};      //Detail - init at zero

    int level;
    double xx=0.0;
    int keep_flag = 1;
    double cell_length=0.0;
    int leaf = 0;//if 1->leaf
    int new_cell = 0;//if 1- this is a newly formed cell

    Cell* parent=nullptr;//Parent

    Cell* left_child=nullptr;
    Cell* right_child=nullptr;

    Cell* left_level=nullptr;
    Cell* right_level=nullptr;

    int virt = 0;

    //Face fluxes

    double f_p[n_eq] = {0.0};//
    double f_m[n_eq] = {0.0};//

     //Runge Kutta storage
    double q_rhs[n_eq] = {0.0};
    double q_temp[n_eq] = {0.0};
    double q_1[n_eq] = {0.0};
    double q_2[n_eq] = {0.0};

    //Source Term
    double q_source[n_eq] = {0.0};

    //Perona Malik Storage
    double filt_q[n_eq]={0.0};

    //Constructor
    Cell(){};
    //Deconstructor
    ~Cell(){};

};


#endif // Global_h_inluded
