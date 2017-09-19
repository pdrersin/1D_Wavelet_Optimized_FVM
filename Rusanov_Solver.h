#ifndef Rusanov_h_inluded
#define Rusanov_h_inluded


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void rusanov_solver(Cell* &cell)
/// \brief Reconstructs fluxes according to Rusanov Approach
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rusanov_solver(Cell* &cellval)
{
    double f_right[n_eq] = {0.0};
    double f_left[n_eq] = {0.0};
    double q_right[n_eq] = {0.0};
    double q_left[n_eq] = {0.0};

    double q[n_eq] = {0.0};
    double qp1[n_eq] = {0.0};
    double qp2[n_eq] = {0.0};
    double qp3[n_eq] = {0.0};
    double qm1[n_eq] = {0.0};
    double qm2[n_eq] = {0.0};
    double qm3[n_eq] = {0.0};

    for (int j=0;j<=n_eq-1;j++)
    {
        q[j] = cellval->q[j];
        qp1[j] = cellval->right_level->q[j];
        qp2[j] = cellval->right_level->right_level->q[j];
        qp3[j] = cellval->right_level->right_level->right_level->q[j];
        qm1[j] = cellval->left_level->q[j];
        qm2[j] = cellval->left_level->left_level->q[j];
        qm3[j] = cellval->left_level->left_level->left_level->q[j];
    }

    double c=0.0;

    if (mflux==1)
    {
        //Calculating right face flux
        weno3(cellval->level,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno3(qm1,q,qp1,qp2);//Wavespeed calculation
    }
    else if (mflux==2)
    {
        //Calculating right face flux
        muscl(qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno3(qm1,q,qp1,qp2);//Wavespeed calculation
    }
    else if (mflux==3)
    {
        //Calculating right face flux
        weno5(cellval->level,qm2,qm1,q,qp1,qp2,qp3,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno5(qm2,qm1,q,qp1,qp2,qp3);//Wavespeed calculation
    }
    else if (mflux==4)
    {
        //Calculating right face flux
        weno6(cellval->level,qm2,qm1,q,qp1,qp2,qp3,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno5(qm2,qm1,q,qp1,qp2,qp3);//Wavespeed calculation
    }



    for (int j=0;j<=n_eq-1;j++)
    {
        cellval->f_p[j] = 0.5*(f_right[j]+f_left[j]) - 0.5*c*(q_right[j]-q_left[j]);//Rusanov reconstruction
    }


    if (mflux==1)
    {
        //Calculating left face flux
        weno3(cellval->level,qm2,qm1,q,qp1,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno3(qm2,qm1,q,qp1);//Wavespeed calculation
    }
    else if (mflux==2)
    {
        //Calculating left face flux
        muscl(qm2,qm1,q,qp1,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno3(qm2,qm1,q,qp1);//Wavespeed calculation
    }
    else if (mflux==3)
    {
        //Calculating left face flux
        weno5(cellval->level,qm3,qm2,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno5(qm3,qm2,qm1,q,qp1,qp2);//Wavespeed calculation
    }
    else if (mflux==4)
    {
        //Calculating left face flux
        weno6(cellval->level,qm3,qm2,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
        c = calc_wave_speed_weno5(qm3,qm2,qm1,q,qp1,qp2);//Wavespeed calculation
    }

    for (int j=0;j<=n_eq-1;j++)
    {
        cellval->f_m[j] = 0.5*(f_right[j]+f_left[j]) - 0.5*c*(q_right[j]-q_left[j]);//Rusanov reconstruction
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double calc_max_wave_speed(unordered_map<int,Cell*> &leafvector)
/// \brief Returns wavespeed for entire leafvector
/// \brief Required for CFL based timestep updation - Rusanov
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double calc_max_wave_speed(unordered_map<int,Cell*> &leafvector)
{
    double max_c = 0.0;

    for (auto citer_one = leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        if ((*citer_one).second->virt==0)
        {
            double q[n_eq]={0.0};
            for(int j=0;j<=n_eq-1;j++)
            {
                q[j] = (*citer_one).second->q[j];
            }

            if (n_eq==1)//Burgers
            {
                double c = fabs(q[0]);
                max_c = max(max_c,fabs(c));
            }
            else if (n_eq==2)//SWE Equations
            {
                double c = swe_wave_speed(q);
                max_c = max(max_c,fabs(c));
            }
            else if (n_eq==3)//Euler Equations
            {
                double c = euler1d_wave_speed(q);
                max_c = max(max_c,fabs(c));
            }
            else if (n_eq==7)
            {
                double c = mhd_wave_speed(q);
                max_c = max(max_c,fabs(c));
            }
            else if (n_eq==14)
            {
                double c = nle_wave_speed(q);
                max_c = max(max_c,fabs(c));
            }
//            else if (n_eq==4)//Euler Equations
//            {
//                double c = elastoplastic_wave_speed(q);
//                max_c = max(max_c,fabs(c));
//            }
        }
    }



    return max_c;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double calc_wave_speed_weno3(double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq])
/// \brief Returns wavespeed at interface of cell i and cell i+1
/// \brief Required for Rusanov solver
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double calc_wave_speed_weno3(double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq])
{
    if (n_eq==1) //Burgers
    {
        double qspeed = fabs(q[0]);
        double qp1speed = fabs(qp1[0]);
        double qp2speed = fabs(qp2[0]);
        double qm1speed = fabs(qm1[0]);

        return max({qspeed,qp1speed,qp2speed,qm1speed});
        //return max({qspeed,qp1speed});

    }
    if (n_eq==2)//SWE Equations
    {
        double qspeed = swe_wave_speed(q);
        double qm1speed = swe_wave_speed(qm1);
        double qp1speed = swe_wave_speed(qp1);
        double qp2speed = swe_wave_speed(qp2);

        return max({qspeed,qp1speed,qp2speed,qm1speed});

    }
    else if (n_eq==3) //Euler
    {
        double qspeed = euler1d_wave_speed(q);
        double qm1speed = euler1d_wave_speed(qm1);
        double qp1speed = euler1d_wave_speed(qp1);
        double qp2speed = euler1d_wave_speed(qp2);

        return max({qspeed,qp1speed,qp2speed,qm1speed});
    }
    else if (n_eq==7)
    {
            double qspeed = mhd_wave_speed(q);
            double qm1speed = mhd_wave_speed(qm1);
            double qp1speed = mhd_wave_speed(qp1);
            double qp2speed = mhd_wave_speed(qp2);

            return max({qspeed,qp1speed,qp2speed,qm1speed});

    }
    else if (n_eq==14)
    {
            double qspeed = nle_wave_speed(q);
            double qm1speed = nle_wave_speed(qm1);
            double qp1speed = nle_wave_speed(qp1);
            double qp2speed = nle_wave_speed(qp2);

            return max({qspeed,qp1speed,qp2speed,qm1speed});
    }
//    else if (n_eq==4) //Euler
//    {
//        double qspeed = elastoplastic_wave_speed(q);
//        double qm1speed = elastoplastic_wave_speed(qm1);
//        double qp1speed = elastoplastic_wave_speed(qp1);
//        double qp2speed = elastoplastic_wave_speed(qp2);
//
//        return max({qspeed,qp1speed,qp2speed,qm1speed});
//    }


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double calc_wave_speed_weno5(double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq])
/// \brief Returns wavespeed at interface of cell i and cell i+1
/// \brief Required for Rusanov solver
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double calc_wave_speed_weno5(double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq])
{
    if (n_eq==1) //Burgers
    {
        double qspeed = fabs(q[0]);
        double qp1speed = fabs(qp1[0]);
        double qp2speed = fabs(qp2[0]);
        double qp3speed = fabs(qp3[0]);
        double qm1speed = fabs(qm1[0]);
        double qm2speed = fabs(qm2[0]);

        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});
        //return max({qspeed,qp1speed});

    }
    else if (n_eq==2)//SW Equations
    {
        double qspeed = swe_wave_speed(q);
        double qm1speed = swe_wave_speed(qm1);
        double qm2speed = swe_wave_speed(qm2);
        double qp1speed = swe_wave_speed(qp1);
        double qp2speed = swe_wave_speed(qp2);
        double qp3speed = swe_wave_speed(qp3);

        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});

    }
    else if (n_eq==3)//Euler Equations
    {
        double qspeed = euler1d_wave_speed(q);
        double qm1speed = euler1d_wave_speed(qm1);
        double qm2speed = euler1d_wave_speed(qm2);
        double qp1speed = euler1d_wave_speed(qp1);
        double qp2speed = euler1d_wave_speed(qp2);
        double qp3speed = euler1d_wave_speed(qp3);

        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});

    }
    else if (n_eq==7)
    {
        double qspeed = mhd_wave_speed(q);
        double qm1speed = mhd_wave_speed(qm1);
        double qm2speed = mhd_wave_speed(qm2);
        double qp1speed = mhd_wave_speed(qp1);
        double qp2speed = mhd_wave_speed(qp2);
        double qp3speed = mhd_wave_speed(qp3);

        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});
    }
    else if (n_eq==14)
    {
        double qspeed = nle_wave_speed(q);
        double qm1speed = nle_wave_speed(qm1);
        double qm2speed = nle_wave_speed(qm2);
        double qp1speed = nle_wave_speed(qp1);
        double qp2speed = nle_wave_speed(qp2);
        double qp3speed = nle_wave_speed(qp3);

        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});
    }
//    else if (n_eq==4)//ElastoPlastic Equations
//    {
//        double qspeed = elastoplastic_wave_speed(q);
//        double qm1speed = elastoplastic_wave_speed(qm1);
//        double qm2speed = elastoplastic_wave_speed(qm1);
//        double qp1speed = elastoplastic_wave_speed(qp1);
//        double qp2speed = elastoplastic_wave_speed(qp2);
//        double qp3speed = elastoplastic_wave_speed(qp2);
//
//        return max({qspeed,qp1speed,qp2speed,qm1speed,qm2speed,qp3speed});
//
//    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double mhd_wave_speed(double (&q)[n_eq])
/// \brief Returns MHD wavespeed at q
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double mhd_wave_speed(double (&q)[n_eq])
{
    double p = (gmm-1.0)*(q[6]-0.5*q[1]*q[1]/q[0]-0.5*q[2]*q[2]/q[0]-0.5*q[3]*q[3]/q[0]-0.5*(bx*bx + q[4]*q[4] + q[5]*q[5]));

    if (p<0.0){p=0.0;}
    if (q[0]<0.0){q[0]=0.0;}

    double a = sqrt(gmm*p/q[0]);
    double ca = sqrt(bx*bx/q[0]);

    double cf = a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0];
    cf = cf + sqrt((a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0])*(a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0]) - 4.0*a*a*bx*bx/q[0]);
    cf = sqrt(0.5*cf);

    double cs = a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0];
    cs = cs - sqrt((a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0])*(a*a + (bx*bx + q[4]*q[4] + q[5]*q[5])/q[0]) - 4.0*a*a*bx*bx/q[0]);
    cs = sqrt(0.5*cs);

    double a1 = fabs(q[1]/q[0] - cf);
    double a2 = fabs(q[1]/q[0] - ca);
    double a3 = fabs(q[1]/q[0] - cs);
    double a4 = fabs(q[1]/q[0]);
    double a5 = fabs(q[1]/q[0] + cs);
    double a6 = fabs(q[1]/q[0] + ca);
    double a7 = fabs(q[1]/q[0] + cf);


    return max({a1,a2,a3,a4,a5,a6,a7});

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double euler1d_wave_speed(double (&q)[n_eq])
/// \brief Returns Euler1D wavespeed at q
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double euler1d_wave_speed(double (&q)[n_eq])
{
    double p = (gmm-1.0)*(q[2]-0.5*q[1]*q[1]/q[0]);
    double a = sqrt(gmm*p/q[0]);

    double l1 = fabs(q[1]/q[0]);
    double l2 = fabs(q[1]/q[0] - a);
    double l3 = fabs(q[1]/q[0] + a);

    return max({l1,l2,l3});

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double swe_wave_speed(double (&q)[n_eq])
/// \brief Returns Shallow Water Equations wavespeed at q
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double swe_wave_speed(double (&q)[n_eq])
{
    double a = sqrt(grav*q[0]);

    double l1 = fabs(q[1]/q[0]);
    double l2 = fabs(q[1]/q[0] - a);
    double l3 = fabs(q[1]/q[0] + a);

    return max({l1,l2,l3});

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn    double nle_wave_speed(double (&q)[n_eq])
/// \brief Returns NLE wavespeed at q - Unfinished
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double nle_wave_speed(double (&q)[n_eq])
{
    double c = 0.0;


    return c;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// \fn    double elastoplastic_wave_speed(double (&q)[n_eq])
///// \brief Returns elasto-plastic problem wavespeed at q - Unfinished
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//double elastoplastic_wave_speed(double (&q)[n_eq])
//{
//    double eta = q[0]/ep_rho_0;
//    double dfdn = (ep_gamma_0*(eta-1.0) + ep_s - eta*(1.0+ep_s))/pow((eta*(ep_s-1.0)-ep_s),3);
//
//    double int_en = q[2]/q[0] - 0.5*q[1]*q[1]/q[0];
//    double mg_f = (eta-1.0)*(eta-ep_gamma_0*(eta-1.0)/2.0)/((eta-ep_s*(eta-1.0))*(eta-ep_s*(eta-1.0)));
//    double p = ep_rho_0*(ep_a_0*ep_a_0)*mg_f + ep_rho_0*ep_gamma_0*int_en;
//
//    double a = ep_a_0*ep_a_0*dfdn + p/(q[0]*q[0])*ep_rho_0*ep_gamma_0;
//
//    double l1 = fabs(q[1]/q[0]);
//    double l2 = fabs(q[1]/q[0] - a);
//    double l3 = fabs(q[1]/q[0] + a);
//
//    return max({l1,l2,l3});
//
//}

#endif // Rusanov_h_inluded
