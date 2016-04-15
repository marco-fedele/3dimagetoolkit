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
  \file poisson_imp.hxx

 \brief file with the implementation of all the methods and functions of namespace \ref lapl
*/

#ifndef POISSON_IMP_HXX_INCLUDED
#define POISSON_IMP_HXX_INCLUDED



template <typename T>
void lapl::NeuGaussSeidel (im3d::image3d<T>& res, im3d::image3d<T> const& b,
                           T const& dt, T const& bc, im3d::image3d<T> const& input)
{

    // we impose external normal derivative equal to bc value on the border

    uint const X = res.getdimx(), Y = res.getdimy(), Z = res.getdimz();
    double const hx = res.gethx(), hy = res.gethy(), hz = res.gethz();

    if (Z == 1)
    {
        double htilde = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy);
        htilde = 1. / htilde;

        double htildei = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy);
        htildei = 1. / htildei;

        double htildej = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy);
        htildej = 1. / htildej;

        double htildeij = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy);
        htildeij = 1. / htildeij;

        // COMPUTE THE OBLIQUE NORMAL COMPONENTS
        T nxy (1.);
        nxy /= sqrt (2.);

        // the computation must be done in the correct order, in res there is the previous
        // step result Gauss-Seidel method is semi-implicit, it uses both current step
        // value and previous step values we read the previous step value in a pixel,
        // we compute the new value and write it in the same pixel the computation
        // involves pixels near to the one on which we write the result
        // (the central pixel) so the pixels with indexes smaller than the indexes of
        // the central pixel will be treated as current step values because at current
        // step we have already computed them the pixels with indexes greater than the
        // indexes of the central pixel will be treated as previous step values because
        // we have not computed them yet we start calculating the result from the pixel
        // (0,0,0)

        // corner i=0, j=0
        res (0, 0, 0) = ( ( res (1, 0, 0) - nxy * bc * hx ) / (hx * hx) +
                          ( res (0, 1, 0) - nxy * bc * hy ) / (hy * hy) +
                          b (0, 0, 0) + res (0, 0, 0) / dt ) * htildeij;

        // edge j=0 , i>0
        for (uint i = 1; i < X - 2; ++i)
            res (i, 0, 0) = ( (res (i + 1, 0, 0) + res (i - 1, 0, 0) ) / (hx * hx) +
                              (res (i, 1, 0) - bc * hy ) / (hy * hy) +
                              b (i, 0, 0) + res (i, 0, 0) / dt  ) * htildej;

        // corner i=X-1 j=0
        res (X - 1, 0, 0) = ( ( res (X - 2, 0, 0) + nxy * bc * hx ) / (hx * hx) +
                              ( res (X - 1, 1, 0) - nxy * bc * hy ) / (hy * hy) +
                              b (X - 1, 0, 0) + res (X - 1, 0, 0) / dt ) * htildeij;

        // INTERNAL NODES
        for (uint j = 1; j < Y - 2; ++j)
        {
            // edge i=0 j>0
            res (0, j, 0) =
                ( ( res (1, j, 0) - bc * hx ) / (hx * hx) + (res (0, j + 1, 0) + res (0, j - 1, 0) ) * (hy * hy)
                  + b (0, j, 0) + res (0, j, 0) / dt   ) * htildei;

            // internal nodes i>0, j>0
            for (uint i = 1; i < X - 1; ++i)
                res (i, j, 0) =
                    ( (res (i + 1, j, 0) + res (i - 1, j, 0) ) / (hx * hx) + (res (i, j + 1, 0) + res (i, j - 1, 0) ) /
                      (hy * hy) + b (i, j, 0) + res (i, j, 0) / dt  ) * htilde;

            // edge i=max, j>0
            res (X - 1, j, 0) =
                ( ( res (X - 2, j, 0) + bc * hx ) / (hx * hx) + ( res (X - 1, j + 1, 0) + res (X - 1, j - 1, 0) ) /
                  (hy * hy) + b (X - 1, j, 0) + res (X - 1, j, 0) / dt ) * htildei;
        }
        // corner i=0, j=Y-1
        res (0, Y - 1, 0) =
            ( ( res (1, Y - 1, 0) - nxy * bc * hx ) / (hx * hx) + ( res (0, Y - 2, 0) + nxy * bc * hy ) /
              (hy * hy) + b (0, Y - 1, 0) + res (0, Y - 1, 0) / dt ) * htildeij;

        //edge i>0, j=Y-1
        for (uint i = 1; i < X - 2; ++i)
            res (i, Y - 1, 0) =
                ( (res (i + 1, Y - 1, 0) + res (i - 1, Y - 1, 0) ) * (hx * hx) + ( res (i, Y - 2, 0) + bc * hy ) /
                  (hy * hy) + b (i, Y - 1, 0) + res (i, Y - 1, 0) / dt ) * htildej;

        // corner i=X-1, j=Y-1
        res (X - 1, Y - 1, 0) =
            ( ( res (X - 2, Y - 1, 0) + nxy * bc * hx ) / (hx * hx) + ( res (X - 1, Y - 2, 0) + nxy * bc * hy ) /
              (hy * hy) + b (X - 1, Y - 1, 0) + res (X - 1, Y - 1, 0) / dt ) * htildeij;

    }
    else
    {
        double htilde = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz);
        htilde = 1. / htilde;

        double htildei = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz);
        htildei = 1. / htildei;

        double htildej = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy) + 2. / (hz * hz);
        htildej = 1. / htildej;

        double htildek = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) + 1. / (hz * hz);
        htildek = 1. / htildek;

        double htildejk = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy) + 1. / (hz * hz);
        htildejk = 1. / htildejk;

        double htildeik = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy) + 1. / (hz * hz);
        htildeik = 1. / htildeik;

        double htildeij = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy) + 2. / (hz * hz);
        htildeij = 1. / htildeij;

        double htildeijk = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy) + 1. / (hz * hz);
        htildeijk = 1. / htildeijk;

        // the computation must be done in the correct order, in res there is the previous
        // step result Gauss-Seidel method is semi-implicit, it uses both curent step
        // value and previous step values we read the previous step value in a pixel,
        // we compute the new value and write it in the same pixel the computation
        // involves pixels near to the one on which we write the result
        // (the central pixel) so the pixels with indexes smaller than the indexes of the
        // central pixel will be treated as current step values because at current step we
        // have already computed them the pixels with indexes greater than the indexes of
        // the central pixel will be treated as previous step values because we have not
        // computed them yet

        // OBLIQUE NORMAL COMPUTATION
        // edge  [(1,1,0) (1,0,1) (0,1,1)]*(1/sqrt(2))
        // angle (1,1,1)*(1/sqrt(3))
        T n_edge (1.), n_angle (1.);
        n_edge /= sqrt (2.);
        n_angle /= sqrt (3.);

        // COMPUTE FACE k=0
        // corner i=0, j=0, k=0
        res (0, 0, 0) = ( (res (1, 0, 0) - n_angle * bc * hx) / (hx * hx) +
                          (res (0, 1, 0) - n_angle * bc * hy) / (hy * hy) +
                          (res (0, 0, 1) - n_angle * bc * hz) / (hz * hz) +
                          b (0, 0, 0) + res (0, 0, 0) / dt ) * htildeijk;
        // edge i>0, j=0, k=0
        for (uint i = 1; i < X - 1; ++i)
            res (i, 0, 0) =  ( (res (i + 1, 0, 0) + res (i - 1, 0, 0) ) / (hx * hx) +
                               (res (i, 1, 0) - n_edge * bc * hy) / (hy * hy) +
                               (res (i, 0, 1) - n_edge * bc * hz) / (hz * hz) +
                               b (i, 0, 0) + res (i, 0, 0) / dt ) * htildejk;
        // corner i=X-1, j=0, k=0
        res (X - 1, 0, 0) = ( (res (X - 2, 0, 0) + n_angle * bc * hx) / (hx * hx) +
                              (res (X - 1, 1, 0) - n_angle * bc * hy) / (hy * hy) +
                              (res (X - 1, 0, 1) - n_angle * bc * hz) / (hz * hz) +
                              b (X - 1, 0, 0) + res (X - 1, 0, 0) / dt ) * htildeijk;
        // COMPUTE INTERNAL NODES OF THE FACE k=0
        for (uint j = 1; j < Y - 1; ++j)
        {
            // edge i=0, j>0, k=0
            res (0, j, 0) = ( (res (1, j, 0) - n_edge * bc * hx) / (hx * hx) +
                              (res (0, j + 1, 0) + res (0, j - 1, 0) ) / (hy * hy) +
                              (res (0, j, 1) - n_edge * bc * hz) / (hz * hz) +
                              b (0, j, 0) + res (0, j, 0) / dt ) * htildeik;
            // face k=0
            // INTERNAL NODES i>0, j>0, k=0
            for (uint i = 1; i < X - 1; ++i)
                res (i, j, 0) = ( ( res (i + 1, j, 0) + res (i - 1, j, 0) ) / (hx * hx) +
                                  ( res (i, j + 1, 0) + res (i, j - 1, 0) ) / (hy * hy) +
                                  ( res (i, j, 1) - bc * hz) / (hz * hz) +
                                  b (i, j, 0) + res (i, j, 0) / dt ) * htildek;
            // edge i=X-1, j>0,  k=0
            res (X - 1, j, 0) = ( (res (X - 2, j, 0) + n_edge * bc * hx) / (hx * hx) +
                                  (res (X - 1, j + 1, 0) + res (X - 1, j - 1, 0) ) / (hy * hy) +
                                  (res (X - 1, j, 1) - n_edge * bc * hz) / (hz * hz) +
                                  b (X - 1, j, 0) + res (X - 1, j, 0) / dt ) * htildeik;
        }
        // corner i=0, j=Y-1, k=0
        res (0, Y - 1, 0) = ( (res (1, Y - 1, 0) - n_angle * bc * hx) / (hx * hx) +
                              (res (0, Y - 2, 0) + n_angle * bc * hy) / (hy * hy) +
                              (res (0, Y - 1, 1) - n_angle * bc * hz) / (hz * hz) +
                              b (0, Y - 1, 0) + res (0, Y - 1, 0) / dt ) * htildeijk;
        // edge i>0, j=Y-1, k=0
        for (uint i = 1; i < X - 1; ++i)
            res (i, Y - 1, 0) = ( ( res (i + 1, Y - 1, 0) + res (i - 1, Y - 1, 0) ) / (hx * hx) +
                                  ( res (i, Y - 2, 0) + n_edge * bc * hy) / (hy * hy) +
                                  ( res (i, Y - 1, 1) - n_edge * bc * hz) / (hz * hz) +
                                  b (i, Y - 1, 0) + res (i, Y - 1, 0) / dt ) * htildejk;
        // corner i=X-1, j=Y-1, k=0
        res (X - 1, Y - 1, 0) = ( (res (X - 2, Y - 1, 0) + n_angle * bc * hx) / (hx * hx) +
                                  (res (X - 1, Y - 2, 0) + n_angle * bc * hy) / (hy * hy) +
                                  (res (X - 1, Y - 1, 1) - n_angle * bc * hz) / (hz * hz) +
                                  b (X - 1, Y - 1, 0) + res (X - 1, Y - 1, 0) / dt ) * htildeijk;

        // COMPUTE THE INTERNAL NODES OF THE CUBE k>0
        for (uint k = 1; k < Z - 1; ++k)
        {
            // edge i=0, j=0, k>0
            res (0, 0, k) = ( (res (1, 0, k) - n_edge * bc * hx) / (hx * hx) +
                              (res (0, 1, k) - n_edge * bc * hy) / (hy * hy) +
                              (res (0, 0, k + 1) + res (0, 0, k - 1) ) / (hz * hz) +
                              b (0, 0, k) + res (0, 0, k) / dt ) * htildeij;
            // face i>0, j=0, k>0
            for (uint i = 1; i < X - 1; ++i)
                res (i, 0, k) = ( (res (i + 1, 0, k) + res (i - 1, 0, k) ) / (hx * hx) +
                                  (res (i, 1, k) - bc * hy) / (hy * hy) +
                                  (res (i, 0, k + 1) + res (i, 0, k - 1) ) / (hz * hz) +
                                  b (i, 0, k) + res (i, 0, k) / dt ) * htildej;
            // edge i=X-1, j=0, k>0
            res (X - 1, 0, k) = ( (res (X - 2, 0, k) + n_edge * bc * hx) / (hx * hx) +
                                  (res (X - 1, 1, k) - n_edge * bc * hy) / (hy * hy) +
                                  (res (X - 1, 0, k + 1) + res (X - 1, 0, k - 1) ) / (hz * hz) +
                                  b (X - 1, 0, k) + res (X - 1, 0, k) / dt ) * htildeij;
            for (uint j = 1; j < Y - 1; ++j)
            {
                // face i=0, j>0, k>0
                res (0, j, k) = ( (res (1, j, k) - bc * hx) / (hx * hx) +
                                  (res (0, j + 1, k) + res (0, j - 1, k) ) / (hy * hy) +
                                  (res (0, j, k + 1) + res (0, j, k - 1) ) / (hz * hz) +
                                  b (0, j, k) + res (0, j, k) / dt ) * htildei;
                // internal nodes i>0, j>0, k>0
                for (uint i = 1; i < X - 1; ++i)
                {
                    res ( i, j, k) = ( ( res (i + 1, j, k) + res (i - 1, j, k) ) / (hx * hx) +
                                       ( res (i, j + 1, k) + res (i, j - 1, k) ) / (hy * hy) +
                                       ( res (i, j, k + 1) + res (i, j, k - 1) ) / (hz * hz) +
                                       b (i, j, k) + res (i, j, k) / dt ) * htilde;
                }
                // face i=dimx-1, j>0, k>0
                res (X - 1, j, k) = ( (res (X - 2, j, k) + bc * hx) / (hx * hx) +
                                      (res (X - 1, j + 1, k) + res (X - 1, j - 1, k) ) / (hy * hy) +
                                      (res (X - 1, j, k + 1) + res (X - 1, j, k - 1) ) / (hz * hz) +
                                      b (X - 1, j, k) + res (X - 1, j, k) / dt ) * htildei;
            }
            // edge i=0, j=Y-1, k>0
            res (0, Y - 1, k) = ( (res (1, Y - 1, k) - n_edge * bc * hx) / (hx * hx) +
                                  (res (0, Y - 2, k) + n_edge * bc * hy) / (hy * hy) +
                                  (res (0, Y - 1, k + 1) + res (0, Y - 1, k - 1) ) / (hz * hz) +
                                  b (0, Y - 1, k) + res (0, Y - 1, k) / dt ) * htildeij;
            // face i>0, j=dimy-1, k>0
            for (uint i = 1; i < X - 1; ++i)
                res (i, Y - 1, k) = ( (res (i + 1, Y - 1, k) + res (i - 1, Y - 1, k) ) / (hx * hx) +
                                      (res (i, Y - 2, k) + bc * hy) / (hy * hy) +
                                      (res (i, Y - 1, k + 1) + res (i, Y - 1, k - 1) ) / (hz * hz) +
                                      b (i, Y - 1, k) + res (i, Y - 1, k) / dt ) * htildej;
            // edge i=X-1, j=Y-1, k>0
            res (X - 1, Y - 1, k) = ( (res (X - 2, Y - 1, k) + n_edge * bc * hx) / (hx * hx) +
                                      (res (X - 1, Y - 2, k) + n_edge * bc * hy) / (hy * hy) +
                                      (res (X - 1, Y - 1, k + 1) + res (X - 1, Y - 1, k - 1) ) / (hz * hz) +
                                      b (X - 1, Y - 1, k) + res (X - 1, Y - 1, k) / dt ) * htildeij;
        }

        // COMPUTE THE FACE k=Z-1

        // corner i=0, j=0, k=Z-1
        res (0, 0, Z - 1) = ( (res (1, 0, Z - 1) - n_angle * bc * hx) / (hx * hx) +
                              (res (0, 1, Z - 1) - n_angle * bc * hy) / (hy * hy) +
                              (res (0, 0, Z - 2) + n_angle * bc * hz) / (hz * hz) +
                              b (0, 0, Z - 1) + res (0, 0, Z - 1) / dt ) * htildeijk;

        // edge i>0, j=0, k=Z-1
        for (uint i = 1; i < X - 1; ++i)
            res (i, 0, Z - 1) = ( (res (i + 1, 0, Z - 1) + res (i - 1, 0, Z - 1) ) / (hx * hx) +
                                  (res (i, 1, Z - 1) - n_edge * bc * hy) / (hy * hy) +
                                  (res (i, 0, Z - 2) + n_edge * bc * hz) / (hz * hz) +
                                  b (i, 0, Z - 1) + res (i, 0, Z - 1) / dt ) * htildejk;

        // corner i=X-1, j=0, k=Z-1
        res (X - 1, 0, Z - 1) = ( (res (X - 2, 0, Z - 1) + n_angle * bc * hx) / (hx * hx) +
                                  (res (X - 1, 1, Z - 1) - n_angle * bc * hy) / (hy * hy) +
                                  (res (X - 1, 0, Z - 2) + n_angle * bc * hz) / (hz * hz) +
                                  b (X - 1, 0, Z - 1) + res (X - 1, 0, Z - 1) / dt ) * htildeijk;

        for (uint j = 1; j < Y - 1; ++j)
        {
            // edge i=0, j>0, k=Z-1
            res (0, j, Z - 1) = ( (res (1, j, Z - 1) - n_edge * bc * hx) / (hx * hx) +
                                  (res (0, j + 1, Z - 1) + res (0, j - 1, Z - 1) ) / (hy * hy) +
                                  (res (0, j, Z - 2) + n_edge * bc * hz) / (hz * hz) +
                                  b (0, j, Z - 1) + res (0, j, Z - 1) / dt ) * htildeik;


            // INTERNAL NODES OF THE FACE k=Z-1
            for (uint i = 1; i < X - 1; ++i)
                res (i, j, Z - 1) = ( (res (i + 1, j, Z - 1) + res (i - 1, j, Z - 1) ) / (hx * hx) +
                                      (res (i, j + 1, Z - 1) + res (i, j - 1, Z - 1) ) / (hy * hy) +
                                      (res (i, j, Z - 2) + bc * hz) / (hz * hz) +
                                      b (i, j, Z - 1) + res (i, j, Z - 1) / dt ) * htildek;
            // edge i=X-1, j>0, k=Z-1
            res (X - 1, j, Z - 1) = ( (res (X - 2, j, Z - 1) + n_edge * bc * hx) / (hx * hx) +
                                      (res (X - 1, j + 1, Z - 1) + res (X - 1, j - 1, Z - 1) ) / (hy * hy) +
                                      (res (X - 1, j, Z - 2) + n_edge * bc * hz) / (hz * hz) +
                                      b (X - 1, j, Z - 1) + res (X - 1, j, Z - 1) / dt) * htildeik;

        }

        // corner i=0, j=Y-1, k=Z-1
        res (0, Y - 1, Z - 1) = ( (res (1, Y - 1, Z - 1) - n_angle * bc * hx) / (hx * hx) +
                                  (res (0, Y - 2, Z - 1) + n_angle * bc * hy) / (hy * hy) +
                                  (res (0, Y - 1, Z - 2) + n_angle * bc * hz) / (hz * hz) +
                                  b (0, Y - 1, Z - 1) + res (0, Y - 1, Z - 1) / dt ) * htildeijk;

        for (uint i = 1; i < X - 1; ++i)
            // edge i>0, j=Y-1, k=Z-1
            res (i, Y - 1, Z - 1) = ( (res (i + 1, Y - 1, Z - 1) + res (i - 1, Y - 1, Z - 1) ) / (hx * hx) +
                                      (res (i, Y - 2, Z - 1) + n_edge * bc * hy) / (hy * hy) +
                                      (res (i, Y - 1, Z - 2) + n_edge * bc * hz) / (hz * hz) +
                                      b (i, Y - 1, Z - 1) + res (i, Y - 1, Z - 1) / dt ) * htildejk;
        // corner i=X-1, j=Y-1, k=Z-1
        res (X - 1, Y - 1, Z - 1) = ( (res (X - 2, Y - 1, Z - 1) + n_angle * bc * hx) / (hx * hx) +
                                      (res (X - 1, Y - 2, Z - 1) + n_angle * bc * hy) / (hy * hy) +
                                      (res (X - 1, Y - 1, Z - 2) + n_angle * bc * hz) / (hz * hz) +
                                      b (X - 1, Y - 1, Z - 1) + res (X - 1, Y - 1, Z - 1) / dt ) * htildeijk;

    }
    return;
}



template <typename T>
void lapl::DirGaussSeidel (im3d::image3d<T>& res, im3d::image3d<T> const& b,
                           T const& dt, T const& bc, im3d::image3d<T> const& input)
{
    uint const X = res.getdimx(), Y = res.getdimy(), Z = res.getdimz();
    double const hx = res.gethx(), hy = res.gethy(), hz = res.gethz();


    // BC ASSIGNEMENT

    // 2D case
    if (Z == 1)
    {
        // edges j=0, k=0 and j=Y-1, k=0
        for (uint i = 0; i < X; ++i)
        {
            res (i, 0, 0) = res (i, Y - 1, 0) = bc;
        }

        // edges i=0, k=0 and i=X-1, k=0
        for (uint j = 0; j < Y; ++j)
        {
            res (0, j, 0) = res (X - 1, j, 0) = bc;
        }

    }

    // 3D case
    else
    {
        // face i=0, i=X-1
        for (uint k = 0; k < Z; ++k)
            for (uint j = 0; j < Y; ++j)
            {
                res (0, j, k) = res (X - 1, j, k) = bc;
            }

        // face j=0, j=Y-1
        for (uint k = 0; k < Z; ++k)
            for (uint i = 0; i < X; ++i)
            {
                res ( i, 0, k) = res ( i, Y - 1, k) = bc;
            }

        // face k=0, k=Z-1
        for (uint j = 0; j < Y; ++j)
            for (uint i = 0; i < X; ++i)
            {
                res (i, j, 0) = res (i, j, Z - 1) = bc;
            }

    }


    // INTERNAL NODES COMPUTATION

    // 2D case
    if (Z == 1)
    {
        double const htilde = 1. / ( 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) );

        for (uint j = 1; j < Y - 1; ++j)
            for (uint i = 1; i < X - 1; ++i)
            {
                res (i, j, 0) =
                    ( (res (i + 1, j, 0) + res (i - 1, j, 0) ) / (hx * hx) +
                      (res (i, j + 1, 0) + res (i, j - 1, 0) ) / (hy * hy) +
                      b (i, j, 0)     + res (i, j, 0) / dt         ) * htilde;
            }
    }

    // 3D case
    else
    {
        double const htilde = 1. / ( 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz) );

        for (uint i = 1; i < X - 1; ++i)
            for (uint j = 1; j < Y - 1; ++j)
                for (uint k = 1; k < Z - 1; ++k)
                {
                    res (i, j, k) =
                        ( (res (i + 1, j, k) + res (i - 1, j, k) ) / (hx * hx) +
                          (res (i, j + 1, k) + res (i, j - 1, k) ) / (hy * hy) +
                          (res (i, j, k + 1) + res (i, j, k - 1) ) / (hz * hz) +
                          b (i, j, k)     + res (i, j, k) / dt         ) * htilde;
                }
    }

    return;
}



#endif // POISSON_IMP_HXX_INCLUDED


