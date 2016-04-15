/*
=========================================================================================
3D Image Toolkit

Copyright 2013:
Marco Fedele, Luca Barbarotta, Francesco Cremonesi, Elena Faggiano.
All rights reserved.
Read details of the license in the file license.txt provided with the library.
=========================================================================================
*/

/*!
  \file  rsfe_splitbregman_imp.hxx

  \brief file with the implementation of all the methods of class \ref segm::rsfe_splitbregman
*/


#ifndef SEGMENTATION_IMP_HXX_INCLUDED
#define SEGMENTATION_IMP_HXX_INCLUDED




template <typename T>
segm::rsfe_splitbregman<T>::rsfe_splitbregman() :
    segmentation<T>(), getpotfile ("data"), verbosity (true), maxiter (100),
    showfrequency (0), dumpfrequency (0),  outputname ("result"), logfilename (""),
    onthego (true), save_current (false), end_now (false),
    sigma (3), lamda (5.e-2), lamda1 (1.e-5), lamda2 (1.e-5),
    beta (100.), dt (1), ls_steps (5), epsilon (1.e-3), a0 (0.), b0 (1.),
    auto_extract_conn_comp (false), cc_init_pixel_x (-1), cc_init_pixel_y (-1),
    cc_init_pixel_z (-1), cc_binary_output (true), cc_init_variables (true), auto_extract_conn_comp_freq (0), bc (NEU),
    ls_tol (0), edge_detector (GRAD), edge_detector_sigma (0), nu (1.),
    onestep_poisson (lapl::DirGaussSeidel), cv (conv::acyclic_fftw_convolution),
    gaussian_pixel_approach (false), tol (1.e-5),
    ix (-1), iy (-1), iz (-1), fx (-1), fy (-1), fz (-1), dimx (0), dimy (0), dimz (1),
    space_dim (3), hx (1), hy (1), hz (1), init_variables (false)
{
    this->alpha = ( this->a0 + this->b0 ) / 2.;
}




// PRIVATE MEMBER IMPLEMENTATION

template <typename T>
void segm::rsfe_splitbregman<T>::initialize_phi_with_cube ()
{
    this->phi.setdim (this->dimx, this->dimy, this->dimz);

    if (this->ix == -1)
    {
        double frac = 0.33;

        this->ix = this->dimx * frac;
        this->iy = this->dimy * frac;
        this->iz = this->dimz * frac;
        this->fx = this->dimx * (1 - frac);
        this->fy = this->dimy * (1 - frac);
        this->fz = this->dimz * (1 - frac);
    }

    // initializing phi with the external value a0
    this->phi = this->a0;


    // modifying values of internal nodes

    // 2d case
    if (this->dimz == 1)
    {
        #pragma omp parallel for
        for (int i = this->ix; i < this->fx; ++i)
            for (int j = this->iy; j < this->fy; ++j)
            {
                this->phi (i, j, 0) = this->b0;
            }
    }
    // 3d case
    else
    {

        #pragma omp parallel for
        for (int i = this->ix; i < this->fx; ++i)
            for (int j = this->iy; j < this->fy; ++j)
                for (int k = this->iz; k < this->fz; ++k)
                {
                    this->phi (i, j, k) = this->b0;
                }
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::initialize_phi_with_init (im3d::image3d<T> const& init)
{
    uint X (init.getdimx() ), Y (init.getdimy() ), Z (init.getdimz() );

    // 1st CASE: init has the same dimensions of the image to segment
    // just assign it to phi
    if ( X == this->dimx && Y == this->dimy && Z == this->dimz )
    {
        this->phi = init;
        this->phi.change_range_of_intensity (this->b0, this->a0);
    }

    // 2nd CASE: in at least one direction init is a pixel smaller
    // assign init to phi and copy results of last but one faces on the last ones
    // (this allows to pass algorithm an interpolation of a result with a lower resolution)
    else if ( (X == this->dimx || X == this->dimx - 1) &&
              (Y == this->dimy || Y == this->dimy - 1) &&
              (Z == this->dimz || Z == this->dimz - 1) )
    {
        this->phi.setdim (this->dimx, this->dimy, this->dimz);

        #pragma omp parallel for
        for (uint i = 0; i < X; ++i)
            for (uint j = 0; j < Y; ++j)
                for (uint k = 0; k < Z; ++k)
                {
                    this->phi (i, j, k) = init (i, j, k);
                }

        // copying result of last but one faces on last threes in case of different
        // dimensions between phi and init
        if (this->dimx != X)
        {
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    this->phi (X, j, k) = this->phi (X - 1, j, k);
                }
        }

        if (this->dimy != Y)
        {
            for (uint i = 0; i < this->dimx; ++i)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    this->phi (i, Y, k) = this->phi (i, Y - 1, k);
                }
        }

        if (this->dimz != Z)
        {
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                {
                    this->phi (i, j, Z) = this->phi (i, j, Z - 1);
                }
        }

        this->phi.change_range_of_intensity (this->b0, this->a0);
    }

    // 3rd CASE: init has incompatible dimensions to initialize phi, cube is used
    else
    {
        this->initialize_phi_with_cube();
        std::clog << "WARNING::apply: wrong dimensions of init, cube initial " <<
                  "contour is used instead." << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_heaviside (im3d::image3d<T>& Heps) const
{
    Heps.setdim (dimx, dimy, dimz);
    Heps.seth (hx, hy, hz);

    T center = (this->a0 + this->b0) / 2.;

    #pragma omp parallel for
    for (uint i = 0; i < dimx; ++i)
        for (uint j = 0; j < dimy; ++j)
            for (uint k = 0; k < dimz; ++k)
                Heps (i, j, k) =
                    0.5 * (1. + 2. / M_PI *
                           atan ( (this->phi (i, j, k) - center)
                                  / (this->epsilon/**std::abs(this->b0-this->a0)*/) ) );

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::edge_detect (im3d::image3d<T>& res, im3d::image3d<T> const& f) const
{
    im3d::image3d<T> aux (f);
    char choice = 'n';
    im3d::interface<T> shower;

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    if (this->verbosity == 1)
    {
        std::clog << "- Setting Edge Detector ..." << std::endl;
    }

    if (this->edge_detector == GRAD || this->edge_detector == GRAD_THRESHOLD || this->dimz != 1)
    {
        std::vector<im3d::image3d<T> > aux_vec;

        if (this->edge_detector_sigma != 0)
        {
            // Build filters of the correct dimensions
            conv::filtering<T>::build_gaussian_filter (aux,
                                                       this->edge_detector_sigma,
                                                       3,
                                                       this->gaussian_pixel_approach);
            cv (aux, f);

            conv::filtering<T>::build_gaussian_filter (aux,
                                                       this->sigma,
                                                       3,
                                                       this->gaussian_pixel_approach);
        }

        aux.grad (aux_vec);

        // case GRAD
        if (this->edge_detector == GRAD)
        {
            // aux = (Modulus of aux_vec)^2
            aux = aux_vec[0] * aux_vec[0];
            aux += aux_vec[1] * aux_vec[1];

            if (this->space_dim == 3)
            {
                aux += aux_vec[2] * aux_vec[2];
            }
        }

        // case GRAD_THRESHOLD
        else
        {
            // aux = Modulus of aux_vec
            im3d::vector_abs (res, aux_vec);

            shower.convertfromimage3d (res);
            std::cout << "showing modulus of gradient to help you to set " <<
                      "the correct threshold ..." << std::endl;
            shower.show_image();

            double t;

            while (choice == 'n')
            {
                std::cout << "Choose a threshold for modulus of gradient:" << std::endl;
                std::cout << "\tthreshold: ";
                cin >> t;
                res.im_to_black_and_white (aux, t);

                shower.convertfromimage3d (aux);
                std::cout << "showing output of chosen edge detector ..." << std::endl;
                shower.show_image();

                std::cout << "Do you accept this edge detector ";
                std::cout << "(if not you can recompute it)? [y/n]\t" << std::endl;
                cin >> choice;
            }
        }

    }// end GRAD/GRAD_THRESHOLD edge detector


    else if (this->edge_detector == LBP)
    {
        T constant;
        double t1, t2;

        while (choice == 'n')
        {
            std::cout << "Choose parameters of the local binary pattern" << std::endl;
            std::cout << "\tconstant T: ";
            cin >> constant;
            std::cout << "\tthreshold t1: ";
            cin >> t1;
            std::cout << "\tthreshold t2 (>t1): ";
            cin >> t2;
            // aux = Local Bynary Pattern Edge Detector Output
            f.local_binary_pattern_edge_detector (aux, constant, t1, t2);

            shower.convertfromimage3d (aux);
            std::cout << "showing output of chosen edge detector ..." << std::endl;

            shower.show_image();

            std::cout << "Do you accept this edge detector ";
            std::cout << "(if not you can recompute it)? [y/n]\t" << std::endl;
            cin >> choice;
        }

    }// end LBP edge detector

    // Build gforshrink
    aux *= this->beta;
    aux += static_cast<T> (1.);
    res = static_cast<T> (1.);
    res /= aux;
    res /= this->lamda;
    res *= this->nu;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::cut_phi ()
{
    for (uint k = 0; k < this->dimz; ++k)
        for (uint j = 0; j < this->dimy; ++j)
            for (uint i = 0; i < this->dimx; ++i)
            {
                if (this->phi (i, j, k) < this->a0)
                {
                    this->phi (i, j, k) = this->a0;
                }

                else if (this->phi (i, j, k) > this->b0)
                {
                    this->phi (i, j, k) = this->b0;
                }
            }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::shrink
( std::vector<im3d::image3d<T> >& res, std::vector<im3d::image3d<T> > const& f1, im3d::image3d<T> const& f2) const
{

    for (uint y = 0; y < this->space_dim; ++y)
        #pragma omp parallel for
        for (uint i = 0; i < dimx; ++i)
            for (uint j = 0; j < dimy; ++j)
                for (uint k = 0; k < dimz; ++k)
                {

                    //try{
                    res[y] (i, j, k) =
                        static_cast<T> ( static_cast<int> (f1[y] (i, j, k) > 0.) - static_cast<int> (f1[y] (i, j, k) < 0.) ) * // sign(f1[y])
                        max ( std::abs (f1[y] (i, j, k) ) - f2 (i, j, k), static_cast<T> (0.) ) ;
                    //test_fpe_exception();
                    //}
                    //catch (BadFPOper & x){std::clog << x.what() << std::endl;}
                    //catch (BadDivision & x){std::clog << x.what() << std::endl;}
                }
    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::internal_show ()
{
    if (this->space_dim == 2)
    {
        this->levelset.show_contour_with_background_image (this->image, this->alpha);
        //this->levelset.show_contour(this->alpha);
        this->levelset.show_image_and_contour (this->alpha);
    }
    else
    {
        this->levelset.show_contour_with_background_image (this->image, this->alpha);
        //this->levelset.show_contour(this->alpha);
        //this->levelset.show_image();
        this->levelset.show_image_and_contour (this->alpha);
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::extract_connected_component ()
{

    im3d::image3d<T> connected_component;

    // first case: manual extraction with user interaction
    if (!auto_extract_conn_comp)
    {
        char choice;

        std::cout << "Do you want to see a connected component of the level set? (y/n): ";
        cin >> choice;

        if (choice == 'y')
        {
            std::cout << "Choose pixel of desired connected component: " << std::endl;

            uint i, j, k;

            this->phi.im_to_black_and_white (connected_component, (this->alpha - this->a0) /
                                             (this->b0 - this->a0) );

            im3d::interface<T> helper (connected_component);

            helper.get_coordinates (i, j, k);

            this->phi.connected_component (connected_component, i, j, k);

            helper.convertfromimage3d (connected_component);
            helper.show_image_and_contour (0.5);

            std::cout << "Do you want to set this connected component " <<
                      "as current level set? (y/n): ";
            cin >> choice;

            if (choice == 'y')
            {
                this->phi = connected_component;
                this->phi.change_range_of_intensity (this->b0, this->a0);
                if (verbosity)
                {
                    std::clog << "-- Updating current levelset to a chosen connected component" << std::endl;
                }

                std::cout << "Do you want to initialize all other variables of the algorithm? (y/n): ";
                cin >> choice;
                if (choice == 'y')
                {
                    this->init_variables = true;
                }

            }
        }

    }

    // second case: automatic extraction using private pixel
    else
    {
        int i = floor ( this->cc_init_pixel_x / this->hx );

        int j = floor ( this->cc_init_pixel_y / this->hy );

        int k = floor ( this->cc_init_pixel_z / this->hz );

        if ( i < 0 || j < 0 || (this->dimz != 1 && k < 0) )
        {
            std::clog << "WARNING::apply::extract_connected_component: " <<
                      "coordinates outside allowed range";
            std::clog << std::endl;
            return;
        }

        uint I, J, K;
        I = static_cast<uint> (i);
        J = static_cast<uint> (j);
        K = static_cast<uint> (k);

        if ( I > this->dimx - 1 || J > this->dimy - 1 || (this->dimz != 1 && K > this->dimz - 1) )
        {
            std::clog << "WARNING::apply::extract_connected_component: " <<
                      "coordinates outside allowed range";
            std::clog  << std::endl;
            return;
        }

        this->phi.im_to_black_and_white (connected_component, (this->alpha - this->a0) / (this->b0 - this->a0) );

        // continue only if private pixel is a white pixel after b&w conversion
        if ( connected_component (I, J, K) == 0 )
        {
            return;
        }

        this->phi.connected_component (connected_component, I, J, K, 0.5, false, this->cc_binary_output);

        // first subcase: automatic setting of connected_component as current levelset
        if ( this->auto_extract_conn_comp_freq != 0 )
        {
            this->phi = connected_component;
            this->phi.change_range_of_intensity (this->b0, this->a0);
            if (this->cc_init_variables)
            {
                this->init_variables = true;
            }
            if (verbosity)
            {
                std::clog << "-- Updating current levelset to a chosen connected component" << std::endl;
            }
        }// end first subcase

        // second subcase: asking user if he wants to set it as current levelset
        else
        {
            im3d::interface<T> helper (connected_component);
            char choice;

            std::cout << "showing a connected component using private pixel to compute it" << std::endl;

            helper.convertfromimage3d (connected_component);
            helper.show_image_and_contour (0.5);

            std::cout << "Do you want to set this connected component as current level set? (y/n): ";
            cin >> choice;

            if (choice == 'y')
            {
                this->phi = connected_component;
                this->phi.change_range_of_intensity (this->b0, this->a0);
                if (verbosity)
                {
                    std::clog << "-- Updating current levelset to a chosen connected component" << std::endl;
                }

                std::cout << "Do you want to initialize all other variables of the algorithm? (y/n): ";
                cin >> choice;
                if (choice == 'y')
                {
                    this->init_variables = true;
                }

            }
        }// end second subcase

    }// end second case


    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::update_param_onthego ()
{
    char choice;

    if (this->onthego == true)
    {
        if (this->verbosity)
        {
            std::clog << "-- Updating parameters of algorithm from file " << this->getpotfile;
            std::clog << " (section splitbregman/onthego/)" << std::endl;
        }
        this->set_param_from_getpot ("onthego/");
    }
    else
    {
        std::cout << "Do you want to save current level set? (y/n): ";
        cin >> choice;

        if (choice == 'y')
        {
            this->save_current = true;
        }

        std::cout << "Do you want to terminate algortihm? (y/n): ";
        cin >> choice;

        if (choice == 'y')
        {
            this->end_now = true;
        }
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_param_from_getpot (std::string const& section)
{
    // preparing GetPot input
    std::string comment_start = "#";
    std::string comment_end   = "\n";
    std::string input_file  = this->getpotfile;

    GetPot ifile (input_file.c_str(), comment_start.c_str(), comment_end.c_str() );

    double ttemp;
    int itemp;
    std::string stemp;

    ifile.set_prefix ( ("splitbregman/" + section).c_str() );

    // General Parameters
    itemp = ifile ("verbosity", -1);
    if (itemp != -1)
    {
        this->set_verbosity ( static_cast<bool> (itemp) );
    }

    itemp = ifile ("maxiter", -1);
    if (itemp > 0)
    {
        this->set_maxiter (itemp);
    }

    itemp = ifile ("showfrequency", -1);
    if (itemp >= 0)
    {
        this->set_showfrequency (itemp);
    }

    itemp = ifile ("dumpfrequency", -1);
    if (itemp >= 0)
    {
        this->set_dumpfrequency (itemp);
    }

    stemp = ifile ("outputname", "");
    if (stemp.size() != 0 )
    {
        this->set_outputname (stemp);
    }

    itemp = ifile ("save_current", -1);
    if (itemp != -1)
    {
        this->save_current = static_cast<bool> (itemp);
    }

    itemp = ifile ("end_now", -1);
    if (itemp != -1)
    {
        this->end_now = static_cast<bool> (itemp);
    }

    // RSFE Parameters
    ttemp = ifile ("sigma", -1.);
    if (ttemp != -1.)
    {
        this->set_sigma (ttemp);
    }

    ttemp = ifile ("lamda", -1.);
    if (ttemp != -1.)
    {
        this->set_lamda (ttemp);
    }

    ttemp = ifile ("lamda1", -1.);
    if (ttemp != -1.)
    {
        this->set_lamda1 (ttemp);
    }

    ttemp = ifile ("lamda2", -1.);
    if (ttemp != -1.)
    {
        this->set_lamda2 (ttemp);
    }

    ttemp = ifile ("nu", -1.);
    if (ttemp != -1.)
    {
        this->set_nu (ttemp);
    }

    ttemp = ifile ("beta", -1.);
    if (ttemp != -1.)
    {
        this->set_beta (ttemp);
    }

    ttemp = ifile ("dt", -1.);
    if (ttemp != -1.)
    {
        this->set_dt (ttemp);
    }

    itemp = ifile ("ls_steps", -1);
    if (itemp != -1)
    {
        this->set_ls_steps ( static_cast<uint> (itemp) );
    }

    ttemp = ifile ("epsilon", -1.);
    if (ttemp != -1.)
    {
        this->set_epsilon (ttemp);
    }

    ttemp = ifile ("a0", 2.1e-6);
    if (ttemp != 2.1e-6)
    {
        this->set_a0 (ttemp);
    }

    ttemp = ifile ("b0", -2.1e-6);
    if (ttemp != -2.1e-6)
    {
        this->set_b0 (ttemp);
    }

    // Connected Component
    itemp = ifile ("auto_extract_conn_comp", -1);
    if (itemp != -1)
    {
        this->set_auto_extract_conn_comp (itemp);
    }

    itemp = ifile ("cc_binary_output", -1);
    if (itemp != -1)
    {
        this->set_cc_binary_output (itemp);
    }

    itemp = ifile ("cc_init_variables", -1);
    if (itemp != -1)
    {
        this->set_cc_init_variables (itemp);
    }

    ttemp = ifile ("cc_init_pixel_x", -1.);
    if (ttemp != -1.)
    {
        this->set_cc_init_pixel_x (ttemp);
    }

    ttemp = ifile ("cc_init_pixel_y", -1.);
    if (ttemp != -1.)
    {
        this->set_cc_init_pixel_y (ttemp);
    }

    ttemp = ifile ("cc_init_pixel_z", -1.);
    if (ttemp != -1.)
    {
        this->set_cc_init_pixel_z (ttemp);
    }

    itemp = ifile ("auto_extract_conn_comp_freq", -1);
    if (itemp >= 0)
    {
        this->set_auto_extract_conn_comp_freq (itemp);
    }

    // Expert
    itemp = ifile ("bc", -1);
    if (itemp != -1)
    {
        this->set_bc ( static_cast<bc_type> (itemp) );
    }

    ttemp = ifile ("ls_tol", -1.);
    if (ttemp != -1.)
    {
        this->set_ls_tol (ttemp);
    }

    itemp = ifile ("edge_detector", -1);
    if (itemp != -1)
    {
        this->set_edge_detector ( static_cast<ed_type> (itemp) );
    }

    ttemp = ifile ("edge_detector_sigma", -1.);
    if (ttemp != -1.)
    {
        this->set_edge_detector_sigma (ttemp);
    }

    itemp = ifile ("gaussian_pixel_approach", -1);
    if (itemp != -1)
    {
        this->set_gaussian_pixel_approach (itemp);
    }

    ttemp = ifile ("tol", -1.);
    if (ttemp != -1.)
    {
        this->set_tol (ttemp);
    }

    ttemp = ifile ("alpha", 1.e-9);
    if (ttemp != 1.e-9)
    {
        this->set_alpha (ttemp);
    }


    return;
}




// PUBLIC MEMBER IMPLEMENTATION

// MEMBERS TO SET PRIVATE PARAMETERS

// NOT CHANGABLE WITH GETPOTFILE
template <typename T>
void segm::rsfe_splitbregman<T>::set_getpotfile (std::string const& name)
{

    this->getpotfile = name;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_onthego (bool const onthego)
{
    this->onthego = onthego;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_logfilename (std::string const& name)
{
    this->logfilename = name;

    return;
}


template <typename T>
void segm::rsfe_splitbregman<T>::set_cv (conv::filtering<T> const cv)
{
    this->cv = cv;

    return;
}


template <typename T>
void segm::rsfe_splitbregman<T>::set_onestep_poisson
(lapl::unsteady_poisson_functor<T> const osl)
{
    this->onestep_poisson = osl;

    return;
}




// CHANGABLE WITH GETPOTFILE


template <typename T>
void segm::rsfe_splitbregman<T>::set_verbosity (bool const v)
{
    this->verbosity = v;
    if (verbosity)
    {
        std::clog << "---\tverbosity = " << this->verbosity << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_maxiter (uint const& maxiter)
{
    this->maxiter = maxiter;
    if (verbosity)
    {
        std::clog << "---\tmaxiter = " << this->maxiter << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_showfrequency (uint const& sf)
{
    this->showfrequency = sf;
    if (verbosity)
    {
        std::clog << "---\tshowfrequency = " << this->showfrequency << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_auto_extract_conn_comp_freq (uint const& cc_freq)
{
    this->auto_extract_conn_comp_freq = cc_freq;
    if (verbosity)
        std::clog << "---\tauto_extract_conn_comp_freq = " <<
                  this->auto_extract_conn_comp_freq << std::endl;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_dumpfrequency (uint const& dump)
{

    this->dumpfrequency = dump;
    if (verbosity)
    {
        std::clog << "---\tdumpfrequency = " << this->dumpfrequency << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_outputname (std::string const& name)
{

    this->outputname = name;
    if (verbosity)
    {
        std::clog << "---\toutputname = " << this->outputname << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_tol (T const& tol)
{

    if (tol >= 0)
    {
        this->tol = tol;
        if (verbosity)
        {
            std::clog << "---\ttol = " << this->tol << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_tol: tol must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_cc_init_pixel_x (double const& cc_x)
{
    if (cc_x >= 0)
    {
        this->cc_init_pixel_x = cc_x;
        if (verbosity)
        {
            std::clog << "---\tcc_init_pixel_x = " << this->cc_init_pixel_x << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_cc_init_pixel_x: coordinates must be positive" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_cc_init_pixel_y (double const& cc_y)
{
    if (cc_y >= 0)
    {
        this->cc_init_pixel_y = cc_y;
        if (verbosity)
        {
            std::clog << "---\tcc_init_pixel_y = " << this->cc_init_pixel_y << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_cc_init_pixel_y: coordinates must be positive" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_cc_init_pixel_z (double const& cc_z)
{
    if (cc_z >= 0)
    {
        this->cc_init_pixel_z = cc_z;
        if (verbosity)
        {
            std::clog << "---\tcc_init_pixel_z = " << this->cc_init_pixel_z << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_cc_init_pixel_z: coordinates must be positive" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_lamda (T const& lamda)
{
    if (lamda > 0)
    {
        this->lamda = lamda;
        if (verbosity)
        {
            std::clog << "---\tlamda = " << this->lamda << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_lamda: lamda must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_lamda1 (T const& lamda1)
{
    if (lamda1 > 0)
    {
        this->lamda1 = lamda1;
        if (verbosity)
        {
            std::clog << "---\tlamda1 = " << this->lamda1 << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_lamda1: lamda1 must be greater than zero" << std::endl;
    }

    return;
}




template <typename T>
void segm::rsfe_splitbregman<T>::set_lamda2 (T const& lamda2)
{
    if (lamda2 > 0)
    {
        this->lamda2 = lamda2;
        if (verbosity)
        {
            std::clog << "---\tlamda2 = " << this->lamda2 << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_lamda2: lamda2 must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_sigma (T const& sigma)
{
    if (sigma > 0)
    {
        this->sigma = sigma;
        if (verbosity)
        {
            std::clog << "---\tsigma = " << this->sigma << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_sigma: sigma must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_edge_detector (ed_type const ed)
{

    if (ed == 1 || ed == 2 || ed == 3)
    {
        this->edge_detector = ed;
    }
    else
    {
        std::clog << "WARNING::set_edge_detector: specify a proper ed_type condition" << std::endl;
    }

    return;

}



template <typename T>
void segm::rsfe_splitbregman<T>::set_edge_detector_sigma (T const& eds)
{
    if (eds >= 0)
    {
        this->edge_detector_sigma = eds;
        if (verbosity)
        {
            std::clog << "---\tedge_detector_sigma = " << this->edge_detector_sigma << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_edge_detector_sigma: " <<
                  "edge_detector_sigma must be greater than zero";
        std::clog << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_nu (T const& nu)
{
    if (nu >= 0)
    {
        this->nu = nu;
        if (verbosity)
        {
            std::clog << "---\tnu = " << this->nu << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_nu: nu must be non negative" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_beta (T const& beta)
{
    if (beta >= 0)
    {
        this->beta = beta;
        if (verbosity)
        {
            std::clog << "---\tbeta = " << this->beta << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_beta: beta must be non negative" << std::endl;
    }

    return;
}


template <typename T>
void segm::rsfe_splitbregman<T>::set_dt (T const& dt)
{

    if (dt > 0)
    {
        this->dt = dt;
        if (verbosity)
        {
            std::clog << "---\tdt = " << this->dt << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_dt: dt must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_a0 (T const& a0)
{
    this->a0 = a0;
    this->alpha = ( this->a0 + this->b0 ) / 2.;
    if (verbosity)
    {
        std::clog << "---\ta0 = " << this->a0 << std::endl;
        std::clog << "---\t  alpha = " << this->alpha << std::endl;
    }

    if (this->b0 < this->a0)
    {
        this->b0 = this->a0 + 1;
        this->alpha = ( this->a0 + this->b0 ) / 2.;

        std::clog << "WARNING::set_a0: b0 has to be greater than a0, hence its value is ";
        std::clog << "automatically modified" << std::endl;
        std::clog << "---\t  b0 = " << this->b0 << std::endl;
        std::clog << "---\t  alpha = " << this->alpha << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_b0 (T const& b0)
{
    this->b0 = b0;
    this->alpha = ( this->a0 + this->b0 ) / 2.;
    if (verbosity)
    {
        std::clog << "---\tb0 = " << this->b0 << std::endl;
        std::clog << "---\t  alpha = " << this->alpha << std::endl;
    }

    if (this->b0 < this->a0)
    {
        this->a0 = this->b0 - 1;
        this->alpha = ( this->a0 + this->b0 ) / 2.;

        std::clog << "WARNING::set_b0: a0 has to be lower than b0, hence its value is ";
        std::clog << "automatically modified" << std::endl;
        std::clog << "---\t  a0 = " << this->a0 << std::endl;
        std::clog << "---\t  alpha = " << this->alpha << std::endl;
    }
    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_epsilon (T const& epsilon)
{
    if (epsilon > 0)
    {
        this->epsilon = epsilon;
        if (verbosity)
        {
            std::clog << "---\tepsilon = " << this->epsilon << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_epsilon: epsilon must be greater than zero" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_bc (bc_type const bc)
{
    if (bc == NEU)
    {
        this->bc = bc;
        this->onestep_poisson.changef (lapl::NeuGaussSeidel);
        if (verbosity)
        {
            std::clog << "---\tbc = NEU" << std::endl;
        }
    }
    else if (bc == DIR)
    {
        this->bc = bc;
        this->onestep_poisson.changef (lapl::DirGaussSeidel);
        if (verbosity)
        {
            std::clog << "---\tbc = DIR" << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_bc: specify a proper bc_type condition" << std::endl;
    }

    return;

}



template <typename T>
void segm::rsfe_splitbregman<T>::set_ls_tol (T const& tol)
{
    if (tol >= 0)
    {
        this->ls_tol = tol;
        if (verbosity)
        {
            std::clog << "---\tls_tol = " << this->ls_tol << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_tol: tol must be positive" << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_ls_steps (uint const steps)
{
    if (steps > 0)
    {
        this->ls_steps = steps;
        if (verbosity)
        {
            std::clog << "---\tls_steps = " << ls_steps << std::endl;
        }
    }
    else
    {
        std::clog << "WARNING::set_ls_steps: ls_steps must be greater than 0" << std::endl;
    }

    return;

}



template <typename T>
void segm::rsfe_splitbregman<T>::set_gaussian_pixel_approach (bool const gpa)
{
    this->gaussian_pixel_approach = gpa;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_auto_extract_conn_comp (bool const auto_cc)
{
    this->auto_extract_conn_comp = auto_cc;
    if (verbosity)
    {
        std::clog << "---\tauto_extract_conn_comp = " << this->auto_extract_conn_comp << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_cc_binary_output (bool const cc_bo)
{
    this->cc_binary_output = cc_bo;
    if (verbosity)
    {
        std::clog << "---\tcc_binary_output = " << this->cc_binary_output << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_cc_init_variables (bool const cc_iv)
{
    this->cc_init_variables = cc_iv;
    if (verbosity)
    {
        std::clog << "---\tcc_init_variables = " << this->cc_init_variables << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_alpha (T const& alpha)
{
    if ( alpha > this->a0 && alpha < this->b0 )
    {
        this->alpha = alpha;
        if (verbosity)
        {
            std::clog << "---\talpha = " << this->alpha << std::endl;
        }
    }

    else
    {
        std::clog << "WARNING::set_alpha: alpha must be greater than " <<
                  "a0 and less than b0" << std::endl;
    }

    return;
}




// MEMBERS TO APPLY ALGORITHM

template <typename T>
void segm::rsfe_splitbregman<T>::initialize_contour_as_cube (im3d::image3d<T> const& image)
{
    im3d::interface<T> myim (image);

    uint x1, x2, y1, y2, z1, z2;

    int X (image.getdimx() ), Y (image.getdimy() ), Z (image.getdimz() );

    if (image.getdimz() == 1)
    {
        std::cout << "Choose the two vertices of initial rectangle you'd like to set." << std::endl;
    }
    else
    {
        std::cout << "Choose the two vertices of initial cube you'd like to set." << std::endl;
    }

    myim.get_coordinates (x1, y1, z1, x2, y2, z2);

    if (x1 > x2)
    {
        this->ix = floor (x2);
        this->fx = floor (x1);
    }
    else
    {
        this->ix = floor (x1);
        this->fx = floor (x2);
    }

    if (y1 > y2)
    {
        this->iy = floor (y2);
        this->fy = floor (y1);
    }
    else
    {
        this->iy = floor (y1);
        this->fy = floor (y2);
    }

    if (z1 > z2)
    {
        this->iz = floor (z2);
        this->fz = floor (z1);
    }
    else
    {
        this->iz = floor (z1);
        this->fz = floor (z2);
    }

    if ( this->ix < 0 || this->ix >= X )
    {
        this->ix = 0;
    }
    if ( this->iy < 0 || this->iy >= Y )
    {
        this->iy = 0;
    }
    if ( this->iz < 0 || this->iz >= Z )
    {
        this->iz = 0;
    }

    if ( this->fx < 0 || this->fx >= X )
    {
        this->fx = image.getdimx() - 1;
    }
    if ( this->fy < 0 || this->fy >= Y )
    {
        this->fy = image.getdimy() - 1;
    }
    if ( this->fz < 0 || this->fz >= Z )
    {
        this->fz = image.getdimz() - 1;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::apply (im3d::image3d<T> const& myim)
{
    this->dimx = myim.getdimx();
    this->dimy = myim.getdimy();
    this->dimz = myim.getdimz();

    this->initialize_phi_with_cube();

    this->apply (myim, this->phi);

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::apply (im3d::image3d<T> const& myim, im3d::image3d<T> const& init)
{
    // redirecting stderr if logfilename is set
    if (this->logfilename != "")
    {
        stderr = freopen (this->logfilename.c_str(), "a", stderr);
    }


    std::clog << "- RSFE/Split-Bregman Algorithm for 3d image segmentation" << std::endl;

    if (this->getpotfile != "")
    {
        if (verbosity)
        {
            std::clog << "-- Initializing parameters of algorithm from " <<
                      "file " << this->getpotfile << " (section splitbregman/init/)" << std::endl;
        }
        this->set_param_from_getpot ("init/");
    }
    else
    {
        this->onthego = false;
    }

    T tollerance = 1000, old_sigma = this->sigma;
    bool old_eds = this->edge_detector_sigma;
    uint t;

    this->dimx = myim.getdimx();
    this->dimy = myim.getdimy();
    this->dimz = myim.getdimz();

    this->hx = myim.gethx();
    this->hy = myim.gethy();
    this->hz = myim.gethz();

    if (this->dimz == 1)
    {
        this->space_dim = 2;
    }
    else
    {
        this->space_dim = 3;
    }

    // begin algorithm from im3d::image3d<T> init if its dimensions are coherent
    this->initialize_phi_with_init (init);

    // adjust phi spacing
    this->phi.seth (this->hx, this->hy, this->hz);

    std::clog << "- Sizeof(phi): " << sizeof (T) *dimx* dimy* dimz << " byte" << std::endl;
    std::clog << "- dimensions:\t" << this->dimx << " x " << this->dimy;
    std::clog << " x " << this->dimz << std::endl;
    std::clog << "- spacing:\t" << this->hx << " x " << this->hy << " x " << this->hz << std::endl;

    //construction of all need variables with correct dimensions and spacing
    im3d::image3d<T> divdb (this->dimx, this->dimy, this->dimz, this->hx, this->hy, this->hz);
    im3d::image3d<T> gforshrink (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> f1 (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> f2 (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> Heps (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> kones (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> kmyim (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> num (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> den (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> kr1 (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> kr2 (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> aux (0, 0, 0, this->hx, this->hy, this->hz);
    im3d::image3d<T> phiold (0, 0, 0, this->hx, this->hy, this->hz);

    std::vector<im3d::image3d<T> > b (this->space_dim, divdb), d;


    //linking intarface image to myim
    this->image.convertfromimage3d (myim);

    for (uint i = 0; i < this->space_dim; ++i)
    {
        b[i] = 0.;
    }

    // showing initial contour
    if ( this->showfrequency > 0 )
    {
        std::cout << "- showing initial contour..." << std::endl;

        this->levelset.convertfromimage3d (this->phi);
        this->levelset.show_contour_with_background_image (this->image, this->alpha);
    }

    // Build filters of the correct dimensions for the algorithm
    conv::filtering<T>::build_gaussian_filter (this->phi,
                                               this->sigma,
                                               3,//radius of filter
                                               this->gaussian_pixel_approach);

    // compute kones and kmyim at the beginning of the algorithm
    // (they don't change during execution)
    aux.setdim (this->dimx, this->dimy, this->dimz, 1.);
    this->cv (kones, aux);
    this->cv (kmyim, myim);
    aux.setdim (0, 0, 0);

    // compute gforshrink using edge_detect private member
    // (coefficient nu/lamda is included)
    this->edge_detect ( gforshrink , myim );


    // algorithm loop
    for (t = 1; tollerance > tol && t <= this->maxiter && end_now == false; ++t)
    {
        if (this->verbosity)
        {
            std::clog << "-- Iteration " << t << std::endl;
        }

        // saving previous iteration result to compute norm of errors
        phiold = this->phi;

        // uncomment all the commented code ending with '\\!!!!!' to compute other norms
        /*
        std::vector<im3d::image3d<T> > b_old(this->space_dim,b[0]);
        b_old[1]=b[1];
        b_old[2]=b[2];
        \\!!!!!
        */



        // -------------------------------------------
        // |   UPDATE RHS OF FIRST MIN PROB ON PHI   |
        // -------------------------------------------

        // if user has modified sigma onthego
        if (old_sigma != this->sigma)
        {
            old_sigma = this->sigma;
            // Build filters with new sigma
            conv::filtering<T>::build_gaussian_filter (this->phi,
                                                       this->sigma,
                                                       3,
                                                       this->gaussian_pixel_approach);
            // compute again kones and kmyim
            aux.setdim (this->dimx, this->dimy, this->dimz, 1.);
            this->cv (kones, aux);
            this->cv (kmyim, myim);
            aux.setdim (0, 0, 0);
        }

        // compute new heaviside using previous iteration result
        this->set_heaviside (Heps);

        // compute auxiliar variable num and den using convolution
        this->cv (num, Heps * myim);
        this->cv (den, Heps);
        // set the dimensions of Heps to zero to limit memory usage
        Heps.setdim (0, 0, 0);

        // compute f1 and f2
        f1 = num / den;
        num -= kmyim; // note: -= limits memory usage
        den -= kones;
        num /= den;
        // set the dimensions of den to zero to limit memory usage
        den.setdim (0, 0, 0);
        f2 = num; // f2 = (kmyim-num)/(kones-den)
        // set the dimensions of num to zero to limit memory usage
        num.setdim (0, 0, 0);

        // compute kr1 and kr2
        this->cv (kr2, f1 * this->lamda1 - f2 * this->lamda2); // factor = f1*lamda1 - f2*lamda2

        f1 *= f1;
        f1 *= this->lamda1; // note: *= limits memory usage
        f2 *= f2;
        f2 *= this->lamda2;
        f1 -= f2;
        // set the dimensions of f2 to zero to limit memory usage
        f2.setdim (0, 0, 0);
        this->cv (kr1, f1); // factor = f1*f1*lamda1 - f2*f2*lamda2
        // set the dimensions of f2 to zero to limit memory usage
        f1.setdim (0, 0, 0);

        // choose bc value depending on kind of boundary condition
        // set in bc private parameters
        double chosen_bc;
        if (this->bc == NEU)
        {
            chosen_bc = 0;
        }
        else if (this->bc == DIR)
        {
            chosen_bc = this->a0;
        }

        // compute one iteration of laplace unsteady equation with forcing function
        // depending on results of convolutions (kones, kmyim, kr1, kr2), myim, private
        // parameters (lamda1, lamda2, lamda) and shrink results at the previous iteration.

        // computing forcing function in kr2
        kr2 *= myim;
        kr2 *= 2;
        aux.setdim (this->dimx, this->dimy, this->dimz, this->lamda1 - this->lamda2 );
        aux *= myim;
        aux *= myim;
        aux *= kones;
        kr2 -= aux;
        kr2 -= kr1;
        kr2 /= this->lamda;
        kr2 -= divdb;
        // set the dimensions of kr1 to zero to limit memory usage
        kr1.setdim (0, 0, 0);



        // -----------------------------------------
        // |   SOLVE MINIMIZATION PROBLEM ON PHI   |
        // -----------------------------------------

        // solve linear system
        // forcing function = (2*kr2*myim-kones*myim*myim*(lamda1-lamda2)-kr1)/lamda-divdb
        T ls_err;
        this->ls_tol ? ls_err = 1000. : ls_err = -1.;
        uint s;
        for (s = 0; s < this->ls_steps || ls_err > this->ls_tol; ++s)
        {
            if ( (this->verbosity && s == this->ls_steps - 1) || this->ls_tol != 0)
            {
                aux = this->phi;
            }
            onestep_poisson (this->phi, kr2, this->dt, chosen_bc);
            if (this->ls_tol != 0 && s >= this->ls_steps - 1)
            {
                ls_err = (this->phi - aux).normL2() / aux.normL2();
            }
        }
        if (this->verbosity)
        {
            if (this->ls_tol == 0)
            {
                ls_err = (this->phi - aux).normL2() / aux.normL2();
            }
            std::clog << "\tlinear system error at " << s << "th (last) step: " << ls_err << std::endl;
        }

        // set the dimensions of kr2 and aux to zero to limit memory usage
        aux.setdim (0, 0, 0);
        kr2.setdim (0, 0, 0);

        // cut all values of phi outside the range [a0,b0]
        this->cut_phi();



        // ------------------------------------------
        // |   SET PHI AS ITS CONNECTED COMPONENT   |
        // ------------------------------------------

        // case total automatic connected component extraction
        if ( this->auto_extract_conn_comp && this->auto_extract_conn_comp_freq != 0 && t % this->auto_extract_conn_comp_freq == 0)
        {
            this->extract_connected_component();
        }


        // showing partial result and allowing user interaction every showfrequency
        if ( this->showfrequency != 0 && t % this->showfrequency == 0 && t != this->maxiter )
        {
            this->levelset.convertfromimage3d (this->phi);
            // showing
            this->internal_show();
            // if not case total automatic
            if ( ! (this->auto_extract_conn_comp && this->auto_extract_conn_comp_freq != 0) )
            {
                // allowing showing of a connected component of phi and setting it as current levelset
                this->extract_connected_component();
            }
            // update parameters of the algorithm on the go
            this->update_param_onthego();
        }

        // compute normL2 of error
        tollerance = (this->phi - phiold).normL2() / phiold.normL2();

        // compute norminf
        if (this->verbosity)
        {
            T inf_err = (this->phi - phiold).norminf() / phiold.norminf();

            std::clog << "\tnorm_inf: " << inf_err << "\tnorm_L2: " << tollerance << std::endl;
            if (this->logfilename != "")
            {
                std::cout << "\tnorm_inf: " << inf_err << "\tnorm_L2: " << tollerance << std::endl;
            }
        }

        // set the dimensions of phiold to zero to limit memory usage
        phiold.setdim (0, 0, 0);



        // --------------------------------------------------
        // |   L1 PART and SPLIT-BREGMAN VARIABLES UPDATE   |
        // --------------------------------------------------

        if (this->init_variables)
        {
            divdb = 0.;
            for (uint i = 0; i < this->space_dim; ++i)
            {
                b[i] = 0.;
            }
            this->init_variables = false;
        }
        else
        {
            // initialize splitbregman variable d:
            // d = grad(phi)
            this->phi.grad (d);

            // partial update splitbregman variable b:
            // b = b_old + d = b_old + grad(phi)
            for (uint i = 0; i < this->space_dim; ++i)
            {
                b[i] += d[i];
            }

            // Update gforshrink if user has just changed edge_detector_sigma
            if ( old_eds != this->edge_detector_sigma)
            {
                edge_detect ( gforshrink , myim );
            }
            old_eds = this->edge_detector_sigma;

            // update splitbregman variable d:
            // d_new = shrink(b, gforshrink) = shrink(b_old+grad(phi), gforshrink)
            shrink (d, b, gforshrink);

            /*
            if(this->verbosity)
            {
                im3d::vector_abs(aux,d);
                if (this->logfilename!="")
                    std::cout << std::endl << "\tnormL1 of |d|: " << aux.normL1();
                std::clog << std::endl << "\tnormL1 of |d|: " << aux.normL1();
            }\\!!!!!
            */

            for (uint i = 0; i < this->space_dim; ++i)
            {
                // update splitbregman variable b:
                // b_new = b-d_new = b_old+grad(phi)-shrink(b_old+grad(phi), gforshrink)
                b[i] -= d[i];
                // d = d_new - b_new
                d[i] -= b[i];
                //b_old[i] -= b[i];\\!!!!!
            }

            /*
            if(this->verbosity)
            {
                im3d::vector_abs(aux,b);

                if (this->logfilename!="")
                    std::cout << "\tnormL1 of |b|: " << aux.normL1() << std::endl;
                std::clog << "\tnormL1 of |b|: " << aux.normL1() << std::endl;

                T max(0.), max1(b_old[0].norminf()), max2(b_old[1].norminf()), max3(b_old[2].norminf());
                if (max1>max2 && max1>max3) max=max1;
                if (max2>max1 && max2>max3) max=max2;
                if (max3>max1 && max3>max2) max=max3;
                if (this->logfilename!="")
                    std::cout << "\tnorminf of |b-bold|: " << max << std::endl;
                std::clog << "\tnorminf of |b-bold|: " << max << std::endl;
            }
            aux.setdim(0,0,0);
            \\!!!!!
            */

            // compute div ( d_new - b_new )
            // this result will be used in the forcing function of unsteady laplace member
            // in next iteration
            im3d::div (divdb, d);

            // set the size of d to zero to limit memory usage
            d.resize (0);
        }


        // saving partial results depending on values of private parameters
        // dumpfrequency and save_current
        if ( (dumpfrequency != 0 && t % dumpfrequency == 0 && t != this->maxiter) ||
                save_current == true )
        {
            this->save_current = false;
            this->levelset.convertfromimage3d (this->phi);
            std::string  aux = this->outputname;
            std::ostringstream oss;
            oss << t;
            aux += "_";
            aux += oss.str();
            this->levelset.write (aux);
        }

    }//end algorithm loop


    this->levelset.convertfromimage3d (this->phi);

    // last show
    if (showfrequency > 0 && this->onthego == false)
    {
        std::cout << "--- showing result at the end of algorithm..." << std::endl;
        this->internal_show();
        this->extract_connected_component();
    }

    // saving final result
    if (dumpfrequency > 0)
    {
        std::string  aux = outputname;
        std::ostringstream oss;
        --t;
        oss << t;
        aux += "_";
        aux += oss.str();
        aux += "_final";
        this->levelset.write (aux);
    }

    // set coordinates of initial rectangle to default values again
    this->ix = -1;
    this->iy = -1;
    this->iz = -1;
    this->fx = -1;
    this->fy = -1;
    this->fz = -1;

    // setting stderr again to default
    if (this->logfilename != "")
    {
        fclose (stderr);
    }

    return;
}



#endif // SEGMENTATION_IMP_HXX_INCLUDED
