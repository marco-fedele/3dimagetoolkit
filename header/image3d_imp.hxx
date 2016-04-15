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
  \file image3d_imp.hxx

   \brief file with the implementation of all the methods of class \ref im3d::image3d
*/

#ifndef IMAGE3D_IMP_HPP_INCLUDED
#define IMAGE3D_IMP_HPP_INCLUDED


#include "interface.hxx"



//CONSTRUCTORS

template <typename T>
im3d::image3d<T>::image3d () :
    rawimage (0), dimx (0), dimy (0), dimz (1),
    hx (0), hy (0), hz (1)
{   }


template <typename T>
im3d::image3d<T>::image3d (uint const& x, uint const& y, uint const& z,
                           double const& hx, double const& hy, double const& hz) :
    dimx (x), dimy (y), dimz (z), hx (hx), hy (hy), hz (hz)
{
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
}


template <typename T>
im3d::image3d<T>::image3d (image3d const& tocopy)
{
    this->setdim (tocopy.getdimx(), tocopy.getdimy(), tocopy.getdimz() );
    this->seth (tocopy.gethx(), tocopy.gethy(), tocopy.gethz() );

    #pragma omp parallel for
    for (uint i = 0; i < dimx; ++i)
        for (uint j = 0; j < dimy; ++j)
            for (uint k = 0; k < dimz; ++k)
            {
                (*this) (i, j, k) = tocopy (i, j, k);
            }
}



//MEMBERS TO GET AND SET PRIVATE PARAMETERS

template <typename T>
void im3d::image3d<T>::setdimx (uint const& x)
{
    this->dimx = x;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdimy (uint const& y)
{
    this->dimy = y;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdimz (uint const& z)
{
    this->dimz = z;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdim (uint const& x, uint const& y, uint const& z,
                               T const& value)
{
    this->dimx = x;
    this->dimy = y;
    this->dimz = z;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz, value);
    return;
}


template <typename T>
void im3d::image3d<T>::sethx (double const& hx)
{
    this->hx = hx;
    return;
}


template <typename T>
void im3d::image3d<T>::sethy (double const& hy)
{
    this->hy = hy;
    return;
}


template <typename T>
void im3d::image3d<T>::sethz (double const& hz)
{
    this->hz = hz;
    return;
}


template <typename T>
void im3d::image3d<T>::seth (double const& hx, double const& hy, double const& hz)
{
    this->hx = hx;
    this->hy = hy;
    this->hz = hz;
    return;
}

//USEFULL OVERLOADING OF VAROIUS OPERATORS

template <typename T>
T im3d::image3d<T>::operator() (uint const& i, uint const& j, uint const& k) const
{
    return rawimage[ this->dimy * this->dimz * i + this->dimz * j + k ];
}


template <typename T>
T& im3d::image3d<T>::operator() (uint const& i, uint const& j, uint const& k)
{
    return rawimage[ this->dimy * this->dimz * i + this->dimz * j + k ];
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator= (S const& toassign)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] = static_cast<T> (toassign);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator= (image3d<S> const& toassign)
{
    this->setdim (toassign.getdimx(), toassign.getdimy(), toassign.getdimz() );
    this->seth (toassign.gethx(), toassign.gethy(), toassign.gethz() );

    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for (uint k = 0; k < this->dimz; ++k)
            {
                (*this) (i, j, k) = static_cast<T> (toassign (i, j, k) );
            }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator += (S const& addend)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] += static_cast<T> (addend);
    }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator += (image3d<S> const& addend)
{
    if (this->dimx == addend.getdimx() &&
            this->dimy == addend.getdimy() &&
            this->dimz == addend.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) += static_cast<T> (addend (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator+= : dimensions must agree" << std::endl;
    }

    return *this;
}




template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator -= (S const& addend)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] -= static_cast<T> (addend);
    }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator -= (image3d<S> const& addend)
{
    if (this->dimx == addend.getdimx() &&
            this->dimy == addend.getdimy() &&
            this->dimz == addend.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) -= static_cast<T> (addend (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator-= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator *= (S const& factor)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] *= static_cast<T> (factor);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator *= (image3d<S> const& factor)
{
    if (this->dimx == factor.getdimx() &&
            this->dimy == factor.getdimy() &&
            this->dimz == factor.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) *= static_cast<T> (factor (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator*= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator /= (S const& factor)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] /= static_cast<T> (factor);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator /= (image3d<S> const& factor)
{
    if (this->dimx == factor.getdimx() &&
            this->dimy == factor.getdimy() &&
            this->dimz == factor.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) /= static_cast<T> (factor (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator/= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator + (image3d<S> const& addend1, R const& addend2)
{
    image3d<S> sum (addend1);
    return sum += addend2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator - (image3d<S> const& addend1, R const& addend2)
{
    image3d<S> sum (addend1);
    return sum -= addend2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator * (image3d<S> const& factor1, R const& factor2)
{
    image3d<S> prod (factor1);
    return prod *= factor2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator / (image3d<S> const& factor1, R const& factor2)
{
    image3d<S> prod (factor1);
    return prod /= factor2;
}



//MEMBERS WORKING ON IMAGE'S VALUES


template <typename S, typename R>
S const im3d::scalarprod ( image3d<S> const& factor1, image3d<R> const& factor2 )
{
    if (factor2.dimx == factor1.dimx &&
            factor2.dimy == factor1.dimy &&
            factor2.dimz == factor1.dimz )
    {
        S result (0);

        #pragma omp parallel for reduction (+:result)
        for (uint i = 0; i < factor1.dimx * factor1.dimy * factor1.dimz; ++i)
        {
            result += factor1.rawimage[i] * static_cast<S> (factor2.rawimage[i]);
        }

        return result;
    }
    else
    {
        std::cout << "WARNING::scalarprod: dimensions must agree, 0 is returned" << std::endl;
        return 0;
    }
}



template <typename S, typename R>
S const im3d::scalarprod_L2 ( image3d<S> const& factor1, image3d<R> const& factor2 )
{

    if (factor2.hx == factor1.hx &&
            factor2.hy == factor1.hy &&
            factor2.hz == factor1.hz )
    {
        return factor1.hx * factor1.hy * factor1.hz * scalarprod (factor1, factor2);
    }
    else
    {
        std::cout << "WARNING::scalarprod_L2: spacing must agree, 0 is returned" << std::endl;
        return 0;
    }
}



template <typename T>
void im3d::image3d<T>::crop (image3d<T>& res,
                             uint const& XSTART, uint const& YSTART, uint const& ZSTART,
                             uint const& XEND, uint const& YEND, uint const& ZEND) const
{

    uint xstart (XSTART), ystart (YSTART), zstart (ZSTART), xend (XEND), yend (YEND), zend (ZEND);

    // check if coordinates are inside the image
    if (xend > this->dimx)
    {
        xend = this->dimx - 1;
    }
    if (xstart > this->dimx)
    {
        xstart = this->dimx - 1;
    }

    if (ystart > this->dimy)
    {
        ystart = this->dimy - 1;
    }
    if (yend > this->dimy)
    {
        yend = this->dimy - 1;
    }


    // check if starts are smaller than ends
    if (xstart > xend)
    {
        uint aux = xstart;
        xstart = xend;
        xend = aux;
    }

    if (ystart > yend)
    {
        uint aux = ystart;
        ystart = yend;
        yend = aux;
    }

    // if 2d, set zstart=zend=0
    if (this->dimz == 1)
    {
        zstart = 0;
        zend = 0;
    }
    // check if z coordinates are inside the image and start is smaller than end
    else
    {
        if (zstart > this->dimz)
        {
            zstart = this->dimz - 1;
        }
        if (zend > this->dimz)
        {
            zend = this->dimz - 1;
        }
        if (zstart > zend)
        {
            uint aux = zstart;
            zstart = zend;
            zend = aux;
        }
    }

    int newdimx = xend - xstart + 1;
    int newdimy = yend - ystart + 1;
    int newdimz = zend - zstart + 1;

    res.setdim ( newdimx, newdimy, newdimz);
    res.seth ( this->hx, this->hy, this->hz);

    #pragma omp parallel for
    for (uint i = 0; i < res.getdimx(); ++i)
        for (uint j = 0; j < res.getdimy(); ++j)
            for (uint k = 0; k < res.getdimz(); ++k)
            {
                res (i, j, k) = (*this) (i + xstart, j + ystart, k + zstart) ;
            }

    return;
}



template <typename T>
void im3d::image3d<T>::crop (image3d<T>& res,
                             double const& XSTART, double const& YSTART, double const& ZSTART,
                             double const& XEND, double const& YEND, double const& ZEND) const
{
    double xstart (XSTART), ystart (YSTART), zstart (ZSTART), xend (XEND), yend (YEND), zend (ZEND);

    // if negative set it to zero
    if (xstart < 0)
    {
        xstart = 0;
    }
    if (ystart < 0)
    {
        ystart = 0;
    }
    if (zstart < 0)
    {
        zstart = 0;
    }
    if (xend < 0)
    {
        xend = 0;
    }
    if (yend < 0)
    {
        yend = 0;
    }
    if (zend < 0)
    {
        zend = 0;
    }

    //  uint xs = floor(dimx*xstart), ys = floor(dimy*ystart), zs = floor(dimz*zstart);
    //  uint xe = ceil(dimx*xend), ye = ceil(dimy*yend), ze = ceil(dimz*zend);

    uint xs = floor (xstart / hx), ys = floor (ystart / hy), zs = floor (zstart / hz);
    uint xe = ceil (xend / hx), ye = ceil (yend / hy), ze = ceil (zend / hz);

    // in crop will be checked the indexes
    this->crop (res, xs, ys, zs, xe, ye, ze);

    return;
}



template <typename T>
void im3d::image3d<T>::change_resolution (image3d<T>& res, uint ratio, bool increase) const
{

    // decreasing resolution
    if (!increase)
    {
        res.setdim (ceil ( static_cast<T> (this->dimx) / static_cast<T> (ratio) ),
                    ceil ( static_cast<T> (this->dimy) / static_cast<T> (ratio) ),
                    ceil ( static_cast<T> (this->dimz) / static_cast<T> (ratio) ) );

        res.seth ( this->hx * ratio, this->hy * ratio, this->hz * ratio );

        #pragma omp parallel for
        for (uint i = 0; i < res.getdimx(); ++i)
            for (uint j = 0; j < res.getdimy(); ++j)
                for (uint k = 0; k < res.getdimz(); ++k)
                {
                    res (i, j, k) = (*this) (i * ratio, j * ratio, k * ratio) ;
                }
    }

    // increasing resolution
    else
    {
        // make ratio a power of 2
        while ( (ratio % 2) != 0 )
        {
            ++ratio;
        }

        image3d<T> resold (*this);

        for (uint r = 2; r <= ratio; r *= 2)
        {
            res.setdim ( (this->dimx - 1) *r + 1, (this->dimy - 1) *r + 1, (this->dimz - 1) *r + 1 );
            res.seth (this->hx / static_cast<double> (r),
                      this->hy / static_cast<double> (r),
                      this->hz / static_cast<double> (r) );

            // copying resold elements every 2 pixels
            #pragma omp parallel for
            for (uint i = 0; i < res.getdimx(); i += 2)
                for (uint j = 0; j < res.getdimy(); j += 2)
                    for (uint k = 0; k < res.getdimz(); k += 2)
                    {
                        res (i, j, k) = resold (i / 2, j / 2, k / 2) ;
                    }


            // 3d case
            if (dimz != 1)
            {

                #pragma omp parallel for
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                        for (uint k = 0; k < res.getdimz() - 2; k += 2)
                        {

                            res (i, j, k + 1) = ( res (i, j, k) + res (i, j, k + 2) ) / 2.;

                            res (i + 1, j, k) = ( res (i, j, k) + res (i + 2, j, k) ) / 2.;

                            res (i, j + 1, k) = ( res (i, j, k) + res (i, j + 2, k) ) / 2.;

                            res (i + 1, j + 1, k + 1) = ( res (i, j, k) + res (i, j, k + 2) +
                                                          res (i, j + 2, k) + res (i, j + 2, k + 2) +
                                                          res (i + 2, j, k) + res (i + 2, j, k + 2) +
                                                          res (i + 2, j + 2, k) + res (i + 2, j + 2, k + 2) ) / 8.;

                            res (i + 1, j, k + 1) = ( res (i, j, k) + res (i, j, k + 2) +
                                                      res (i + 2, j, k) + res (i + 2, j, k + 2) ) / 4.;

                            res (i + 1, j + 1, k) = ( res (i, j, k) + res (i + 2, j, k) +
                                                      res (i, j + 2, k) + res (i + 2, j + 2, k) ) / 4.;

                            res (i, j + 1, k + 1) = ( res (i, j, k) + res (i, j + 2, k) +
                                                      res (i, j, k + 2) + res (i, j + 2, k + 2) ) / 4.;
                        }


                //  complete interpolation on the last three faces
                //  of the domain with an even index
                uint X (res.getdimx() - 1), Y (res.getdimy() - 1), Z (res.getdimz() - 1);

                // face i=X
                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    for (uint k = 0; k < res.getdimz() - 2; k += 2)
                    {
                        res (X, j, k + 1) = ( res (X, j, k) + res (X, j, k + 2) ) / 2.;
                        res (X, j + 1, k) = ( res (X, j, k) + res (X, j + 2, k) ) / 2.;
                        res (X, j + 1, k + 1) = ( res (X, j, k) + res (X, j + 2, k) +
                                                  res (X, j, k + 2) + res (X, j + 2, k + 2) ) / 4.;
                    }

                // face j=Y
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint k = 0; k < res.getdimz() - 2; k += 2)
                    {
                        res (i, Y, k + 1) = ( res (i, Y, k) + res (i, Y, k + 2) ) / 2.;
                        res (i + 1, Y, k) = ( res (i, Y, k) + res (i + 2, Y, k) ) / 2.;
                        res (i + 1, Y, k + 1) = ( res (i, Y, k) + res (i, Y, k + 2) +
                                                  res (i + 2, Y, k) + res (i + 2, Y, k + 2) ) / 4.;
                    }

                // face k=Z
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    {
                        res (i + 1, j, Z) = ( res (i, j, Z) + res (i + 2, j, Z) ) / 2.;
                        res (i, j + 1, Z) = ( res (i, j, Z) + res (i, j + 2, Z) ) / 2.;
                        res (i + 1, j + 1, Z) = ( res (i, j, Z) + res (i + 2, j, Z) +
                                                  res (i, j + 2, Z) + res (i + 2, j + 2, Z) ) / 4.;
                    }

                // three common edges
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                {
                    res (i + 1, Y, Z) = ( res (i, Y, Z) + res (i + 2, Y, Z) ) / 2.;
                }

                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                {
                    res (X, j + 1, Z) = ( res (X, j, Z) + res (X, j + 2, Z) ) / 2.;
                }

                for (uint k = 0; k < res.getdimz() - 2; k += 2)
                {
                    res (X, Y, k + 1) = ( res (X, Y, k) + res (X, Y, k + 2) ) / 2.;
                }

            }// end 3d case

            // 2d case
            else
            {
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    {

                        res (i + 1, j, 0) = ( res (i, j, 0) + res (i + 2, j, 0) ) / 2.;

                        res (i, j + 1, 0) = ( res (i, j, 0) + res (i, j + 2, 0) ) / 2.;

                        res (i + 1, j + 1, 0) = ( res (i, j, 0) + res (i + 2, j, 0) +
                                                  res (i, j + 2, 0) + res (i + 2, j + 2, 0) ) / 4.;
                    }

                //  complete interpolation on the last two edges
                //  of the domain with an even index
                uint X (res.getdimx() - 1), Y (res.getdimy() - 1);

                //  edge i=X
                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                {
                    res (X, j + 1, 0) = ( res (X, j, 0) + res (X, j + 2, 0) ) / 2.;
                }

                //  edge j=Y
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                {
                    res (i + 1, Y, 0) = ( res (i, Y, 0) + res (i + 2, Y, 0) ) / 2.;
                }

            }// end 2d case

            resold = res;

        }// end for on r = 2...ratio


    }// end increasing resolution
    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::grad (std::vector<image3d<S> >& res) const
{

    // CASE OF 3D IMAGES
    if (dimz > 1)
    {

        if (res.size() != 3)
        {
            res.resize (3);
        }
        for (uint i = 0; i < 3; ++i)
        {
            res[i].setdim (this->dimx, this->dimy, this->dimz);
            res[i].seth (this->hx, this->hy, this->hz);
        }

        #pragma omp parallel
        {
            // starting parallel section

            #pragma omp for
            for (uint i = 1; i < this->dimx - 1; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    for (uint k = 0; k < this->dimz; ++k)
                        res[0] (i, j, k) =
                            (static_cast<S> ( (*this) (i + 1, j, k) ) -
                             static_cast<S> ( (*this) (i - 1, j, k) ) ) / (2.*hx) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 1; j < this->dimy - 1; ++j)
                    for (uint k = 0; k < this->dimz; ++k)
                        res[1] (i, j, k) =
                            (static_cast<S> ( (*this) (i, j + 1, k) ) -
                             static_cast<S> ( (*this) (i, j - 1, k) ) ) / (2.*hy) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    for (uint k = 1; k < this->dimz - 1; ++k)
                        res[2] (i, j, k) =
                            ( static_cast<S> ( (*this) (i, j, k + 1) ) -
                              static_cast<S> ( (*this) (i, j, k - 1) ) ) / (2.*hz) ;

        }   // ending parallel section

        // boundary values

        for (uint j = 0; j < this->dimy; ++j)
            for (uint k = 0; k < this->dimz; ++k)
            {
                res[0] (0, j, k) = ( 4 * static_cast<S> ( (*this) (1, j, k) ) -
                                     3 * static_cast<S> ( (*this) (0, j, k) ) -
                                     static_cast<S> ( (*this) (2, j, k) ) ) / (2.*hx) ;
                res[0] (this->dimx - 1, j, k) = (3 * static_cast<S> ( (*this) (dimx - 1, j, k) ) -
                                                 4 * static_cast<S> ( (*this) (dimx - 2, j, k) ) +
                                                 static_cast<S> ( (*this) (dimx - 3, j, k) ) ) / (2.*hx);
            }

        for (uint i = 0; i < this->dimx; ++i)
            for (uint k = 0; k < this->dimz; ++k)
            {
                res[1] (i, 0, k) = ( 4 * static_cast<S> ( (*this) (i, 1, k) ) -
                                     3 * static_cast<S> ( (*this) (i, 0, k) ) -
                                     static_cast<S> ( (*this) (i, 2, k) ) ) / (2.*hy) ;
                res[1] (i, this->dimy - 1, k) = (3 * static_cast<S> ( (*this) (i, dimy - 1, k) ) -
                                                 4 * static_cast<S> ( (*this) (i, dimy - 2, k) ) +
                                                 static_cast<S> ( (*this) (i, dimy - 3, k) ) ) / (2.*hy);
            }

        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
            {
                res[2] (i, j, 0) = ( 4 * static_cast<S> ( (*this) (i, j, 1) ) -
                                     3 * static_cast<S> ( (*this) (i, j, 0) ) -
                                     static_cast<S> ( (*this) (i, j, 2) ) ) / (2.*hz) ;
                res[2] (i, j, this->dimz - 1) = (3 * static_cast<S> ( (*this) (i, j, dimz - 1) ) -
                                                 4 * static_cast<S> ( (*this) (i, j, dimz - 2) ) +
                                                 static_cast<S> ( (*this) (i, j, dimz - 3) ) ) / (2.*hz);
            }
    } // end if(dimz>1)

    // CASE OF 2D IMAGES
    else
    {

        if (res.size() != 2)
        {
            res.resize (2);
        }
        for (uint i = 0; i < 2; ++i)
        {
            res[i].setdim (this->dimx, this->dimy, this->dimz);
            res[i].seth (this->hx, this->hy, this->hz);
        }

        #pragma omp parallel
        {
            // starting parallel section

            #pragma omp for
            for (uint i = 1; i < this->dimx - 1; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    res[0] (i, j, 0) =
                        ( static_cast<S> ( (*this) (i + 1, j, 0) ) -
                          static_cast<S> ( (*this) (i - 1, j, 0) ) ) / (2.*hx) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 1; j < this->dimy - 1; ++j)
                    res[1] (i, j, 0) =
                        ( static_cast<S> ( (*this) (i, j + 1, 0) ) -
                          static_cast<S> ( (*this) (i, j - 1, 0) ) ) / (2.*hy) ;

        }   // ending parallel section


        // boundary values

        for (uint j = 0; j < this->dimy; ++j)
        {
            res[0] (0, j, 0) = ( 4 * static_cast<S> ( (*this) (1, j, 0) ) -
                                 3 * static_cast<S> ( (*this) (0, j, 0) ) -
                                 static_cast<S> ( (*this) (2, j, 0) ) ) / (2.*hx) ;
            res[0] (this->dimx - 1, j, 0) = ( 3 * static_cast<S> ( (*this) (dimx - 1, j, 0) ) -
                                              4 * static_cast<S> ( (*this) (dimx - 2, j, 0) ) +
                                              static_cast<S> ( (*this) (dimx - 3, j, 0 ) ) ) / (2.*hx);
        }

        for (uint i = 0; i < this->dimx; ++i)
        {
            res[1] (i, 0, 0) = ( 4 * static_cast<S> ( (*this) (i, 1, 0) ) -
                                 3 * static_cast<S> ( (*this) (i, 0, 0) ) -
                                 static_cast<S> ( (*this) (i, 2, 0) ) ) / (2.*hy) ;
            res[1] (i, this->dimy - 1, 0) = ( 3 * static_cast<S> ( (*this) (i, dimy - 1, 0) ) -
                                              4 * static_cast<S> ( (*this) (i, dimy - 2, 0) ) +
                                              static_cast<S> ( (*this) (i, dimy - 3, 0) ) ) / (2.*hy);
        }

    }// end else ( CASE OF 2D IMAGES)

    return;
}



template <typename S, typename R>
void im3d::vector_abs (image3d<S>& res, std::vector<image3d<R> > const& fun)
{
    // check dimensions of input vector
    if (fun.size() != 2 && fun.size() != 3)
    {
        std::cout << "In image3d::vector_abs: input vector could be only 2d or 3d, ";
        std::cout << fun.size() << " is a not allowed dimension." << std::endl;
        return;
    }

    res.setdim (fun[0].dimx, fun[0].dimy, fun[0].dimz );
    res.seth (fun[0].hx, fun[0].hy, fun[0].hz );

    // 3d case
    if (fun.size() == 3 &&
            fun[0].dimx == fun[1].dimx && fun[1].dimx == fun[2].dimx  &&
            fun[0].dimy == fun[1].dimy && fun[1].dimy == fun[2].dimy  &&
            fun[0].dimz == fun[1].dimz && fun[1].dimz == fun[2].dimz )
    {
        #pragma omp parallel for
        for (uint i = 0; i < res.dimx * res.dimy * res.dimz; ++i)
            res.rawimage[i] =
                static_cast<S> (sqrt ( fun[0].rawimage[i] * fun[0].rawimage[i] +
                                       fun[1].rawimage[i] * fun[1].rawimage[i] +
                                       fun[2].rawimage[i] * fun[2].rawimage[i] ) );
        //note: sqrt works only with float, double or long double
    }

    // 2d case
    else if (fun.size() == 2 &&
             fun[0].dimx == fun[1].dimx &&
             fun[0].dimy == fun[1].dimy &&
             fun[0].dimz == fun[1].dimz )
    {
        #pragma omp parallel for
        for (uint i = 0; i < res.dimx * res.dimy * res.dimz; ++i)
            res.rawimage[i] =
                static_cast<S> (sqrt ( fun[0].rawimage[i] * fun[0].rawimage[i] +
                                       fun[1].rawimage[i] * fun[1].rawimage[i] ) );
        //note: sqrt works only with float, double or long double
    }

    else
    {
        std::cout << "In image3d::vector_abs: dimensions of image3d " <<
                  "inside input vector must agree" << std::endl;
        return;
    }

    return;
}



template <typename S, typename R>
void im3d::div (image3d<S>& res, std::vector<image3d<R> > const& fun)
{

    // check dimensions of input vector
    if (fun.size() != 2 && fun.size() != 3)
    {
        std::cout << "In image3d::div: input vector could be only 2d or 3d, ";
        std::cout << fun.size() << " is a not allowed dimension." << std::endl;
        return;
    }

    res.setdim (fun[0].getdimx(), fun[0].getdimy(), fun[0].getdimz() );
    res.seth (fun[0].gethx(), fun[0].gethy(), fun[0].gethz() );

    double hx = res.gethx(), hy = res.gethy(), hz = res.gethz();
    uint X = res.getdimx(), Y = res.getdimy(), Z = res.getdimz();

    // 2d case
    if (fun.size() == 2 &&
            fun[0].getdimx() == fun[1].getdimx() &&
            fun[0].getdimy() == fun[1].getdimy() &&
            fun[0].getdimz() == fun[1].getdimz() )
    {

        #pragma omp parallel for
        for (uint i = 1; i < X - 1; ++i)
            for (uint j = 1; j < Y - 1; ++j)
                res (i, j, 0) =
                    ( (static_cast<S> (fun[0] (i + 1, j, 0) ) -
                       static_cast<S> (fun[0] (i - 1, j, 0) ) ) / (2.*hx) +
                      (static_cast<S> (fun[1] (i, j + 1, 0) ) -
                       static_cast<S> (fun[1] (i, j - 1, 0) ) ) / (2.*hy) );

        //edges i-variabili
        for (uint i = 1; i < X - 1; ++i)
        {
            //j=0, k=0
            res (i, 0, 0) =
                ( (static_cast<S> (fun[0] (i + 1, 0, 0) ) -
                   static_cast<S> (fun[0] (i - 1, 0, 0) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (i, 1, 0) ) -
                   3 * static_cast<S> (fun[1] (i, 0, 0) ) -
                   static_cast<S> (fun[1] (i, 2, 0) ) ) / (2.*hy) );

            //j=Y-1, k=0
            res (i, Y - 1, 0) =
                ( (static_cast<S> (fun[0] (i + 1, Y - 1, 0) ) -
                   static_cast<S> (fun[0] (i - 1, Y - 1, 0) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (i, Y - 1, 0) ) -
                   4 * static_cast<S> (fun[1] (i, Y - 2, 0) ) +
                   static_cast<S> (fun[1] (i, Y - 3, 0) ) ) / (2.*hy) );
        }

        //edges j variabili
        for (uint j = 1; j < Y - 1; ++j)
        {
            //i=0, k=0
            res (0, j, 0) =
                ( (4 * static_cast<S> (fun[0] (1, j, 0) ) -
                   3 * static_cast<S> (fun[0] (0, j, 0) ) -
                   static_cast<S> (fun[0] (2, j, 0) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (0, j + 1, 0) ) -
                   static_cast<S> (fun[1] (0, j - 1, 0) ) ) / (2.*hy) );

            //i=X-1, k=0
            res (X - 1, j, 0) =
                ( (3 * static_cast<S> (fun[0] (X - 1, j, 0) ) -
                   4 * static_cast<S> (fun[0] (X - 2, j, 0) ) +
                   static_cast<S> (fun[0] (X - 3, j, 0) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (X - 1, j + 1, 0) ) -
                   static_cast<S> (fun[1] (X - 1, j - 1, 0) ) ) / (2.*hy) );
        }

        //angles

        //i=0, j=0, k=0
        res (0, 0, 0) =
            ( (4 * static_cast<S> (fun[0] (1, 0, 0) ) -
               3 * static_cast<S> (fun[0] (0, 0, 0) ) -
               static_cast<S> (fun[0] (2, 0, 0) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (0, 1, 0) ) -
               3 * static_cast<S> (fun[1] (0, 0, 0) ) -
               static_cast<S> (fun[1] (0, 2, 0) ) ) / (2.*hy) );

        //i=X-1, j=0, k=0
        res (X - 1, 0, 0) =
            ( (3 * static_cast<S> (fun[0] (X - 1, 0, 0) ) -
               4 * static_cast<S> (fun[0] (X - 2, 0, 0) ) +
               static_cast<S> (fun[0] (X - 3, 0, 0) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (X - 1, 1, 0) ) -
               3 * static_cast<S> (fun[1] (X - 1, 0, 0) ) -
               static_cast<S> (fun[1] (X - 1, 2, 0) ) ) / (2.*hy) );

        //i=0, j=Y-1, k=0
        res (0, Y - 1, 0) =
            ( (4 * static_cast<S> (fun[0] (1, Y - 1, 0) ) -
               3 * static_cast<S> (fun[0] (0, Y - 1, 0) ) -
               static_cast<S> (fun[0] (2, Y - 1, 0) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (0, Y - 1, 0) ) -
               4 * static_cast<S> (fun[1] (0, Y - 2, 0) ) +
               static_cast<S> (fun[1] (0, Y - 3, 0) ) ) / (2.*hy) );

        //i=X-1, j=Y-1, k=0
        res (X - 1, Y - 1, 0) =
            ( (3 * static_cast<S> (fun[0] (X - 1, Y - 1, 0) ) -
               4 * static_cast<S> (fun[0] (X - 2, Y - 1, 0) ) +
               static_cast<S> (fun[0] (X - 3, Y - 1, 0) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (X - 1, Y - 1, 0) ) -
               4 * static_cast<S> (fun[1] (X - 1, Y - 2, 0) ) +
               static_cast<S> (fun[1] (X - 1, Y - 3, 0) ) ) / (2.*hy) );
    }

    // 3d case
    else if (fun.size() == 3 &&
             fun[0].getdimx() == fun[1].getdimx() && fun[1].getdimx() == fun[2].getdimx()  &&
             fun[0].getdimy() == fun[1].getdimy() && fun[1].getdimy() == fun[2].getdimy()  &&
             fun[0].getdimz() == fun[1].getdimz() && fun[1].getdimz() == fun[2].getdimz() )
    {

        #pragma omp parallel for
        for (uint i = 1; i < X - 1; ++i)
            for (uint j = 1; j < Y - 1; ++j)
                for (uint k = 1; k < Z - 1; ++k)
                    res (i, j, k) =
                        ( (static_cast<S> (fun[0] (i + 1, j, k) ) -
                           static_cast<S> (fun[0] (i - 1, j, k) ) ) / (2.*hx) +
                          (static_cast<S> (fun[1] (i, j + 1, k) ) -
                           static_cast<S> (fun[1] (i, j - 1, k) ) ) / (2.*hy) +
                          (static_cast<S> (fun[2] (i, j, k + 1) ) -
                           static_cast<S> (fun[2] (i, j, k - 1) ) ) / (2.*hz) );

        //face i=0 and i=X-1
        for (uint j = 1; j < Y - 1; ++j)
            for (uint k = 1; k < Z - 1; ++k)
            {
                res (0, j, k) =
                    ( (4 * static_cast<S> (fun[0] (1, j, k) ) -
                       3 * static_cast<S> (fun[0] (0, j, k) ) -
                       static_cast<S> (fun[0] (2, j, k) ) ) / (2.*hx) +
                      (static_cast<S> (fun[1] (0, j + 1, k) ) -
                       static_cast<S> (fun[1] (0, j - 1, k) ) ) / (2.*hy) +
                      (static_cast<S> (fun[2] (0, j, k + 1) ) -
                       static_cast<S> (fun[2] (0, j, k - 1) ) ) / (2.*hz) );

                res (X - 1, j, k) =
                    ( (3 * static_cast<S> (fun[0] (X - 1, j, k) ) -
                       4 * static_cast<S> (fun[0] (X - 2, j, k) ) +
                       static_cast<S> (fun[0] (X - 3, j, k) ) ) / 2.* (hx) +
                      (static_cast<S> (fun[1] (X - 1, j + 1, k) ) -
                       static_cast<S> (fun[1] (X - 1, j - 1, k) ) ) / (2.*hy) +
                      (static_cast<S> (fun[2] (X - 1, j, k + 1) ) -
                       static_cast<S> (fun[2] (X - 1, j, k - 1) ) ) / (2.*hz) );
            }

        //face j=0 and j=Y-1
        for (uint i = 1; i < X - 1; ++i)
            for (uint k = 1; k < Z - 1; ++k)
            {
                res (i, 0, k) =
                    ( (static_cast<S> (fun[0] (i + 1, 0, k) ) -
                       static_cast<S> (fun[0] (i - 1, 0, k) ) ) / (2.*hx) +
                      (4 * static_cast<S> (fun[1] (i, 1, k) ) -
                       3 * static_cast<S> (fun[1] (i, 0, k) ) -
                       static_cast<S> (fun[1] (i, 2, k) ) ) / (2.*hy) +
                      (static_cast<S> (fun[2] (i, 0, k + 1) ) -
                       static_cast<S> (fun[2] (i, 0, k - 1) ) ) / (2.*hz) );

                res (i, Y - 1, k) =
                    ( (static_cast<S> (fun[0] (i + 1, Y - 1, k) ) -
                       static_cast<S> (fun[0] (i - 1, Y - 1, k) ) ) / (2.*hx) +
                      (3 * static_cast<S> (fun[1] (i, Y - 1, k) ) -
                       4 * static_cast<S> (fun[1] (i, Y - 2, k) ) +
                       static_cast<S> (fun[1] (i, Y - 3, k) ) ) / (2.*hy) +
                      (static_cast<S> (fun[2] (i, Y - 1, k + 1) ) -
                       static_cast<S> (fun[2] (i, Y - 1, k - 1) ) ) / (2.*hz) );
            }

        //face k=0 and k=Z-1
        for (uint i = 1; i < X - 1; ++i)
            for (uint j = 1; j < Y - 1; ++j)
            {
                res (i, j, 0) =
                    ( (static_cast<S> (fun[0] (i + 1, j, 0) ) -
                       static_cast<S> (fun[0] (i - 1, j, 0) ) ) / (2.*hx) +
                      (static_cast<S> (fun[1] (i, j + 1, 0) ) -
                       static_cast<S> (fun[1] (i, j - 1, 0) ) ) / (2.*hy) +
                      (4 * static_cast<S> (fun[2] (i, j, 1) ) -
                       3 * static_cast<S> (fun[2] (i, j, 0) ) -
                       static_cast<S> (fun[2] (i, j, 2) ) ) / (2.*hz) );

                res (i, j, Z - 1) =
                    ( (static_cast<S> (fun[0] (i + 1, j, Z - 1) ) -
                       static_cast<S> (fun[0] (i - 1, j, Z - 1) ) ) / (2.*hx) +
                      (static_cast<S> (fun[1] (i, j + 1, Z - 1) ) -
                       static_cast<S> (fun[1] (i, j - 1, Z - 1) ) ) / (2.*hy) +
                      (3 * static_cast<S> (fun[2] (i, j, Z - 1) ) -
                       4 * static_cast<S> (fun[2] (i, j, Z - 2) ) +
                       static_cast<S> (fun[2] (i, j, Z - 3) ) ) / (2.*hz) );
            }

        //edges i-variabili
        for (uint i = 1; i < X - 1; ++i)
        {

            //j=0, k=0
            res (i, 0, 0) =
                ( (static_cast<S> (fun[0] (i + 1, 0, 0) ) -
                   static_cast<S> (fun[0] (i - 1, 0, 0) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (i, 1, 0) ) -
                   3 * static_cast<S> (fun[1] (i, 0, 0) ) -
                   static_cast<S> (fun[1] (i, 2, 0) ) ) / (2.*hy) +
                  (4 * static_cast<S> (fun[2] (i, 0, 1) ) -
                   3 * static_cast<S> (fun[2] (i, 0, 0) ) -
                   static_cast<S> (fun[2] (i, 0, 2) ) ) / (2.*hz) );

            //j=Y-1, k=0
            res (i, Y - 1, 0) =
                ( (static_cast<S> (fun[0] (i + 1, Y - 1, 0) ) -
                   static_cast<S> (fun[0] (i - 1, Y - 1, 0) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (i, Y - 1, 0) ) -
                   4 * static_cast<S> (fun[1] (i, Y - 2, 0) ) +
                   static_cast<S> (fun[1] (i, Y - 3, 0) ) ) / (2.*hy) +
                  (4 * static_cast<S> (fun[2] (i, Y - 1, 1) ) -
                   3 * static_cast<S> (fun[2] (i, Y - 1, 0) ) -
                   static_cast<S> (fun[2] (i, Y - 1, 2) ) ) / (2.*hz) );

            //j=Y-1, k=Z-1
            res (i, Y - 1, Z - 1) =
                ( (static_cast<S> (fun[0] (i + 1, Y - 1, Z - 1) ) -
                   static_cast<S> (fun[0] (i - 1, Y - 1, Z - 1) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (i, Y - 1, Z - 1) ) -
                   4 * static_cast<S> (fun[1] (i, Y - 2, Z - 1) ) +
                   static_cast<S> (fun[1] (i, Y - 3, Z - 1) ) ) / (2.*hy) +
                  (3 * static_cast<S> (fun[2] (i, Y - 1, Z - 1) ) -
                   4 * static_cast<S> (fun[2] (i, Y - 1, Z - 2) ) +
                   static_cast<S> (fun[2] (i, Y - 1, Z - 3) ) ) / (2.*hz) );

            //j=0, k=Z-1
            res (i, 0, Z - 1) =
                ( (static_cast<S> (fun[0] (i + 1, 0, Z - 1) ) -
                   static_cast<S> (fun[0] (i - 1, 0, Z - 1) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (i, 1, Z - 1) ) -
                   3 * static_cast<S> (fun[1] (i, 0, Z - 1) ) -
                   static_cast<S> (fun[1] (i, 2, Z - 1) ) ) / (2.*hy) +
                  (3 * static_cast<S> (fun[2] (i, 0, Z - 1) ) -
                   4 * static_cast<S> (fun[2] (i, 0, Z - 2) ) +
                   static_cast<S> (fun[2] (i, 0, Z - 2) ) ) / (2.*hz) );
        }

        //edges j variabili
        for (uint j = 1; j < Y - 1; ++j)
        {

            //i=0, k=0
            res (0, j, 0) =
                ( (4 * static_cast<S> (fun[0] (1, j, 0) ) -
                   3 * static_cast<S> (fun[0] (0, j, 0) ) -
                   static_cast<S> (fun[0] (2, j, 0) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (0, j + 1, 0) ) -
                   static_cast<S> (fun[1] (0, j - 1, 0) ) ) / (2.*hy) +
                  (4 * static_cast<S> (fun[2] (0, j, 1) ) -
                   3 * static_cast<S> (fun[2] (0, j, 0) ) -
                   static_cast<S> (fun[2] (0, j, 2) ) ) / (2.*hz) );

            //i=X-1, k=0
            res (X - 1, j, 0) =
                ( (3 * static_cast<S> (fun[0] (X - 1, j, 0) ) -
                   4 * static_cast<S> (fun[0] (X - 2, j, 0) ) +
                   static_cast<S> (fun[0] (X - 3, j, 0) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (X - 1, j + 1, 0) ) -
                   static_cast<S> (fun[1] (X - 1, j - 1, 0) ) ) / (2.*hy) +
                  (4 * static_cast<S> (fun[2] (X - 1, j, 1) ) -
                   3 * static_cast<S> (fun[2] (X - 1, j, 0) ) -
                   static_cast<S> (fun[2] (X - 1, j, 2) ) ) / (2.*hz) );

            //i=X-1, k=Z-1
            res (X - 1, j, Z - 1) =
                ( (3 * static_cast<S> (fun[0] (X - 1, j, Z - 1) ) -
                   4 * static_cast<S> (fun[0] (X - 2, j, Z - 1) ) +
                   static_cast<S> (fun[0] (X - 3, j, Z - 1) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (X - 1, j + 1, Z - 1) ) -
                   static_cast<S> (fun[1] (X - 1, j - 1, Z - 1) ) ) / (2.*hy) +
                  (3 * static_cast<S> (fun[2] (X - 1, j, Z - 1) ) -
                   4 * static_cast<S> (fun[2] (X - 1, j, Z - 2) ) +
                   static_cast<S> (fun[2] (X - 1, j, Z - 3) ) ) / (2.*hz) );

            //i=0, k=Z-1
            res (0, j, Z - 1) =
                ( (4 * static_cast<S> (fun[0] (1, j, Z - 1) ) -
                   3 * static_cast<S> (fun[0] (0, j, Z - 1) ) -
                   static_cast<S> (fun[0] (2, j, Z - 1) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (0, j + 1, Z - 1) ) -
                   static_cast<S> (fun[1] (0, j - 1, Z - 1) ) ) / (2.*hy) +
                  (3 * static_cast<S> (fun[2] (0, j, Z - 1) ) -
                   4 * static_cast<S> (fun[2] (0, j, Z - 2) ) +
                   static_cast<S> (fun[2] (0, j, Z - 3) ) ) / (2.*hz) );
        }

        //edges k variabili
        for (uint k = 1; k < Z - 1; ++k)
        {

            //i=0, j=0
            res (0, 0, k) =
                ( (4 * static_cast<S> (fun[0] (1, 0, k) ) -
                   3 * static_cast<S> (fun[0] (0, 0, k) ) -
                   static_cast<S> (fun[0] (2, 0, k) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (0, 1, k) ) -
                   3 * static_cast<S> (fun[1] (0, 0, k) ) -
                   static_cast<S> (fun[1] (0, 2, k) ) ) / (2.*hy) +
                  (static_cast<S> (fun[2] (0, 0, k + 1) ) -
                   static_cast<S> (fun[2] (0, 0, k - 1) ) ) / (2.*hz) );

            //i=X-1, j=0
            res (X - 1, 0, k) =
                ( (3 * static_cast<S> (fun[0] (X - 1, 0, k) ) -
                   4 * static_cast<S> (fun[0] (X - 2, 0, k) ) +
                   static_cast<S> (fun[0] (X - 3, 0, k) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (X - 1, 1, k) ) -
                   3 * static_cast<S> (fun[1] (X - 1, 0, k) ) -
                   static_cast<S> (fun[1] (X - 1, 2, k) ) ) / (2.*hy) +
                  (static_cast<S> (fun[2] (X - 1, 0, k + 1) ) -
                   static_cast<S> (fun[2] (X - 1, 0, k - 1) ) ) / (2.*hz) );

            //i=X-1, j=Y-1
            res (X - 1, Y - 1, k) =
                ( (3 * static_cast<S> (fun[0] (X - 1, Y - 1, k) ) -
                   4 * static_cast<S> (fun[0] (X - 2, Y - 1, k) ) +
                   static_cast<S> (fun[0] (X - 3, Y - 1, k) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (X - 1, Y - 1, k) ) -
                   4 * static_cast<S> (fun[1] (X - 1, Y - 2, k) ) +
                   static_cast<S> (fun[1] (X - 1, Y - 3, k) ) ) / (2.*hy) +
                  (static_cast<S> (fun[2] (X - 1, Y - 1, k + 1) ) -
                   static_cast<S> (fun[2] (X - 1, Y - 1, k - 1) ) ) / (2.*hz) );

            //i=0, j=Y-1
            res (0, Y - 1, k) =
                ( (4 * static_cast<S> (fun[0] (1, Y - 1, k) ) -
                   3 * static_cast<S> (fun[0] (0, Y - 1, k) ) -
                   static_cast<S> (fun[0] (2, Y - 1, k) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (0, Y - 1, k) ) -
                   4 * static_cast<S> (fun[1] (0, Y - 2, k) ) +
                   static_cast<S> (fun[1] (0, Y - 3, k) ) ) / (2.*hy) +
                  (static_cast<S> (fun[2] (0, Y - 1, k + 1) ) -
                   static_cast<S> (fun[2] (0, Y - 1, k - 1) ) ) / (2.*hz) );

        }

        //angles

        //i=0, j=0, k=0
        res (0, 0, 0) =
            ( (4 * static_cast<S> (fun[0] (1, 0, 0) ) - 3 * static_cast<S> (fun[0] (0, 0, 0) ) -
               static_cast<S> (fun[0] (2, 0, 0) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (0, 1, 0) ) - 3 * static_cast<S> (fun[1] (0, 0, 0) ) -
               static_cast<S> (fun[1] (0, 2, 0) ) ) / (2.*hy) +
              (4 * static_cast<S> (fun[2] (0, 0, 1) ) - 3 * static_cast<S> (fun[2] (0, 0, 0) ) -
               static_cast<S> (fun[2] (0, 0, 2) ) ) / (2.*hz) );

        //i=X-1, j=0, k=0
        res (X - 1, 0, 0) =
            ( (3 * static_cast<S> (fun[0] (X - 1, 0, 0) ) - 4 * static_cast<S> (fun[0] (X - 2, 0, 0) ) +
               static_cast<S> (fun[0] (X - 3, 0, 0) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (X - 1, 1, 0) ) - 3 * static_cast<S> (fun[1] (X - 1, 0, 0) ) -
               static_cast<S> (fun[1] (X - 1, 2, 0) ) ) / (2.*hy) +
              (4 * static_cast<S> (fun[2] (X - 1, 0, 1) ) - 3 * static_cast<S> (fun[2] (X - 1, 0, 0) ) -
               static_cast<S> (fun[2] (X - 1, 0, 2) ) ) / (2.*hz) );

        //i=0, j=Y-1, k=0
        res (0, Y - 1, 0) =
            ( (4 * static_cast<S> (fun[0] (1, Y - 1, 0) ) - 3 * static_cast<S> (fun[0] (0, Y - 1, 0) ) -
               static_cast<S> (fun[0] (2, Y - 1, 0) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (0, Y - 1, 0) ) - 4 * static_cast<S> (fun[1] (0, Y - 2, 0) ) +
               static_cast<S> (fun[1] (0, Y - 3, 0) ) ) / (2.*hy) +
              (4 * static_cast<S> (fun[2] (0, Y - 1, 1) ) - 3 * static_cast<S> (fun[2] (0, Y - 1, 0) ) -
               static_cast<S> (fun[2] (0, Y - 1, 2) ) ) / (2.*hz) );

        //i=0, j=0, k=Z-1
        res (0, 0, Z - 1) =
            ( (4 * static_cast<S> (fun[0] (1, 0, Z - 1) ) - 3 * static_cast<S> (fun[0] (0, 0, Z - 1) ) -
               static_cast<S> (fun[0] (2, 0, Z - 1) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (0, 1, Z - 1) ) - 3 * static_cast<S> (fun[1] (0, 0, Z - 1) ) -
               static_cast<S> (fun[1] (0, 2, Z - 1) ) ) / (2.*hy) +
              (3 * static_cast<S> (fun[2] (0, 0, Z - 1) ) - 4 * static_cast<S> (fun[2] (0, 0, Z - 2) ) +
               static_cast<S> (fun[2] (0, 0, Z - 3) ) ) / (2.*hz) );

        //i=0, j=Y-1, k=Z-1
        res (0, Y - 1, Z - 1) =
            ( (4 * static_cast<S> (fun[0] (1, Y - 1, Z - 1) ) - 3 * static_cast<S> (fun[0] (0, Y - 1, Z - 1) ) -
               static_cast<S> (fun[0] (2, Y - 1, Z - 1) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (0, Y - 1, Z - 1) ) - 4 * static_cast<S> (fun[1] (0, Y - 2, Z - 1) ) +
               static_cast<S> (fun[1] (0, Y - 3, Z - 1) ) ) / (2.*hy) +
              (3 * static_cast<S> (fun[2] (0, Y - 1, Z - 1) ) - 4 * static_cast<S> (fun[2] (0, Y - 1, Z - 2) ) +
               static_cast<S> (fun[2] (0, Y - 1, Z - 3) ) ) / (2.*hz) );

        //i=X-1, j=Y-1, k=0
        res (X - 1, Y - 1, 0) =
            ( (3 * static_cast<S> (fun[0] (X - 1, Y - 1, 0) ) - 4 * static_cast<S> (fun[0] (X - 2, Y - 1, 0) ) +
               static_cast<S> (fun[0] (X - 3, Y - 1, 0) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (X - 1, Y - 1, 0) ) - 4 * static_cast<S> (fun[1] (X - 1, Y - 2, 0) ) +
               static_cast<S> (fun[1] (X - 1, Y - 3, 0) ) ) / (2.*hy) +
              (4 * static_cast<S> (fun[2] (X - 1, Y - 1, 1) ) - 3 * static_cast<S> (fun[2] (X - 1, Y - 1, 0) ) -
               static_cast<S> (fun[2] (X - 1, Y - 1, 2) ) ) / (2.*hz) );

        //i=X-1, j=0, k=Z-1
        res (X - 1, 0, Z - 1) =
            ( (3 * static_cast<S> (fun[0] (X - 1, 0, Z - 1) ) - 4 * static_cast<S> (fun[0] (X - 2, 0, Z - 1) ) +
               static_cast<S> (fun[0] (X - 3, 0, Z - 1) ) ) / (2.*hx) +
              (4 * static_cast<S> (fun[1] (X - 1, 1, Z - 1) ) - 3 * static_cast<S> (fun[1] (X - 1, 0, Z - 1) ) -
               static_cast<S> (fun[1] (X - 1, 2, Z - 1) ) ) / (2.*hy) +
              (3 * static_cast<S> (fun[2] (X - 1, 0, Z - 1) ) - 4 * static_cast<S> (fun[2] (X - 1, 0, Z - 2) ) +
               static_cast<S> (fun[2] (X - 1, 0, Z - 3) ) ) / (2.*hz) );

        //i=X-1, j=Y-1, k=Z-1
        res (X - 1, Y - 1, Z - 1) =
            ( (3 * static_cast<S> (fun[0] (X - 1, Y - 1, Z - 1) ) -
               4 * static_cast<S> (fun[0] (X - 2, Y - 1, Z - 1) ) +
               static_cast<S> (fun[0] (X - 3, Y - 1, Z - 1) ) ) / (2.*hx) +
              (3 * static_cast<S> (fun[1] (X - 1, Y - 1, Z - 1) ) -
               4 * static_cast<S> (fun[1] (X - 1, Y - 2, Z - 1) ) +
               static_cast<S> (fun[1] (X - 1, Y - 3, Z - 1) ) ) / (2.*hy) +
              (3 * static_cast<S> (fun[2] (X - 1, Y - 1, Z - 1) ) -
               4 * static_cast<S> (fun[2] (X - 1, Y - 1, Z - 2) ) +
               static_cast<S> (fun[2] (X - 1, Y - 1, Z - 3) ) ) / (2.*hz) );

    }

    else
    {
        std::cout << "In image3d::div: dimensions of image3d inside " <<
                  "input vector must agree" << std::endl;
        return;
    }

    return;
}



template <typename T>
T const im3d::image3d<T>::norm1 () const
{
    T norm = 0;

    #pragma omp parallel for reduction (+:norm)
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        norm += std::abs (this->rawimage[i]);
    }

    return norm;
}



template <typename T>
T const im3d::image3d<T>::normL1 () const
{
    return this->norm1() * this->hx * this->hy * this->hz;
}



template <typename T>
T const im3d::image3d<T>::norm2 () const
{
    return sqrt (scalarprod (*this, *this) );
}



template <typename T>
T const im3d::image3d<T>::normL2 () const
{
    return sqrt (scalarprod_L2 (*this, *this) );
}



template <typename T>
T const im3d::image3d<T>::norminf () const
{
    T norminf = std::abs ( (*this) (0, 0, 0) );

    #pragma omp parallel
    {

        T private_norminf = norminf;

        #pragma omp for nowait
        for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
            if ( private_norminf < this->rawimage[i] )
            {
                private_norminf = std::abs (this->rawimage[i]);
            }

        #pragma omp critical
        {
            norminf = (private_norminf > norminf) ? private_norminf : norminf ;
        }

    }

    return norminf;
}


template <typename T>
T const im3d::image3d<T>::max () const
{
    T shared_max ( (*this) (0, 0, 0) );

    #pragma omp parallel
    {

        T max (shared_max);

        #pragma omp for nowait
        for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
            if ( max < this->rawimage[i] )
            {
                max = this->rawimage[i];
            }

        #pragma omp critical
        {
            shared_max = (max > shared_max) ? max : shared_max ;
        }

    }

    return shared_max;
}



template <typename T>
T const im3d::image3d<T>::min () const
{
    T shared_min ( (*this) (0, 0, 0) );

    #pragma omp parallel
    {

        T min (shared_min);

        #pragma omp for nowait
        for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
            if ( min > this->rawimage[i] )
            {
                min = this->rawimage[i];
            }

        #pragma omp critical
        {
            shared_min = (min < shared_min) ? min : shared_min ;
        }

    }

    return shared_min;
}



template <typename T>
void im3d::image3d<T>::histogram_equalization (image3d<T>& res,
                                               uint const& quantization) const
{

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    // finding max and min to create bands
    T max (this->max() ), min (this->min() );

    // initializing a vector to store frequency of each band
    std::vector<int> frequency (quantization, 0);

    // computing width of each band
    double band_width (static_cast<double> (floor (max - min) + 1) / static_cast<double> (quantization) );

    // computing frequency of each band
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for ( uint k = 0; k < this->dimz; ++k)
            {
                uint index =
                    static_cast<uint> ( floor (static_cast<double> ( (*this) (i, j, k) - min)
                                               / band_width) );
                #pragma omp atomic
                ++ (frequency[index]);
            }

    // transforming it in the vector of cumulative frequency
    for (uint x = 1; x < quantization; ++x)
    {
        frequency[x] += frequency[x - 1];
    }

    // computing cumulative probability
    int elem_num (this->dimx * this->dimy * this->dimz);
    std::vector<double> cumulative (quantization, 0.);

    #pragma omp parallel for
    for (uint x = 0; x < quantization; ++x)
    {
        cumulative[x] = static_cast<double> (frequency[x]) / static_cast<double> (elem_num);
    }

    // modifying original image replacing original intensity with the new one computed
    // using cumulative probability vector
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for ( uint k = 0; k < this->dimz; ++k)
            {
                uint index =
                    static_cast<uint> ( floor (static_cast<double> ( (*this) (i, j, k) - min)
                                               / band_width) );

                res (i, j, k) =
                    static_cast<T> (cumulative[index] * static_cast<double> (quantization) );
            }
    return;
}



template <typename T>
void im3d::image3d<T>::histogram_equalization (image3d<T>& res) const
{

    T max (this->max() ), min (this->min() );

    this->histogram_equalization (res, floor ( static_cast<double> (max - min) ) + 1 );

    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::im_to_black_and_white (image3d<S>& res,
                                              double threshold, bool negative) const
{

    if (threshold > 1 || threshold < 0)
    {
        std::cout << "WARNING::im_to_black_and_white: threshold must be between 0 and 1, ";
        std::cout << "otherwise default value 0.5 is applied" << std::endl;
        threshold = 0.5;
    }

    T max (this->max() ), min (this->min() ), pivot (min + (max - min) *threshold);

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for ( uint k = 0; k < this->dimz; ++k)
                if ( (*this) (i, j, k) < pivot )
                {
                    res (i, j, k) = 0 + static_cast<int> (negative);
                }
                else
                {
                    res (i, j, k) = 1 - static_cast<int> (negative);
                }

    return;
}



template <typename T>
void im3d::image3d<T>::change_range_of_intensity (T const& max, T const& min)
{
    (*this) -= this->min();
    (*this) /= this->max() / (max - min);
    (*this) += min;

    return;
}



template <typename T>
void im3d::image3d<T>::select_range_of_intensity (image3d<T>& res, T const& lowerbound,
                                                  T const& upperbound, int type, T lowervalue, T uppervalue) const
{
    if ( lowerbound > upperbound )
    {
        std::cout << "WARNING::select_range_of_intensity: upperbound must be greater";
        std::cout << " than lowerbound" << std::endl;
        return;
    }

    T low_external_value = lowerbound;
    T up_external_value = upperbound;

    if (type > 0)
    {
        low_external_value = upperbound;
    }
    else if (type == -1)
    {
        up_external_value = lowerbound;
    }
    else if (type == -2)
    {
        low_external_value = lowervalue;
        up_external_value = uppervalue;
    }

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for ( uint k = 0; k < this->dimz; ++k)
            {
                if ( (*this) (i, j, k) < lowerbound )
                {
                    res (i, j, k) = low_external_value;
                }
                else if ( (*this) (i, j, k) > upperbound )
                {
                    res (i, j, k) = up_external_value;
                }
                else
                {
                    res (i, j, k) = (*this) (i, j, k);
                }
            }

    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::connected_component (image3d<S>& res,
                                            uint const& i, uint const& j, uint const& k,
                                            double threshold, bool full_connected, bool binary_output) const
{

    image3d<S> bw;
    uint I (i), J (j), K (k);

    // initialize res and resold as a black image and norm1!=0
    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    if (I == 0)
    {
        ++I;
    }
    if (J == 0)
    {
        ++J;
    }
    if (I >= this->dimx - 1)
    {
        I = this->dimx - 2;
    }
    if (J >= this->dimy - 1)
    {
        J = this->dimy - 2;
    }
    if (this->dimz == 1)
    {
        K = 0;
    }
    else
    {
        if (K == 0)
        {
            ++K;
        }
        if (K >= this->dimz - 1)
        {
            K = this->dimz - 2;
        }
    }

    this->im_to_black_and_white (bw, threshold);

    if ( bw (I, J, K) == 0 )
    {
        this->im_to_black_and_white (bw, threshold, true);
    }

    res (I, J, K) = 1;

    this->connected_component (res, bw, full_connected);

    if (!binary_output)
    {
        this->cc_not_binary_output (res, threshold);
    }

    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::connected_component (image3d<S>& res, image3d<S> const& init,
                                            double threshold, bool full_connected, bool binary_output) const
{

    if (init.getdimx() != this->dimx ||
            init.getdimy() != this->dimy ||
            init.getdimz() != this->dimz )
    {
        std::cout << "WARNING: in image3d::connected_component: dimensions of init must be ";
        std::cout << "the same of the image" << std::endl;
        return;
    }

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    image3d<S> bw;

    this->im_to_black_and_white (bw, threshold);

    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for (uint k = 0; k < this->dimz; ++k)
                if ( init (i, j, k) && bw (i, j, k) )
                {
                    res (i, j, k) = static_cast<S> (1);
                }

    this->connected_component (res, bw, full_connected);

    if (!binary_output)
    {
        this->cc_not_binary_output (res, threshold);
    }

    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::cc_not_binary_output (image3d<S>& res, double threshold) const
{
    image3d<S> bw_negative;
    this->im_to_black_and_white (bw_negative, threshold, true);

    res += bw_negative; // image equal to 1 in cc and in negative output of the black and white threshold

    res *= *this - this->min();
    res += this->min();

    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::connected_component (image3d<S>& res, image3d<S>& bw,
                                            bool const& full_connected) const
{
    image3d<S> resold (this->dimx, this->dimy, this->dimz);

    uint norm1 = 1;

    // 2d case
    if (this->dimz == 1)
    {
        while ( norm1 != 0 )
        {
            for (uint i = 1; i < this->dimx - 1; ++i)
                for (uint j = 1; j < this->dimy - 1; ++j)
                    if ( res (i, j, 0) == 1)
                    {
                        // east
                        res (i + 1, j, 0) = 1 * bw (i + 1, j, 0);
                        // north
                        res (i, j + 1, 0) = 1 * bw (i, j + 1, 0);
                        if (full_connected)
                        {
                            // north-east
                            res (i + 1, j + 1, 0) = 1 * bw (i + 1, j + 1, 0);
                        }
                    }

            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = dimy - 2; j > 0; --j)
                    if ( res (i, j, 0) == 1)
                    {
                        // west
                        res (i - 1, j, 0) = 1 * bw (i - 1, j, 0);
                        // south
                        res (i, j - 1, 0) = 1 * bw (i, j - 1, 0);
                        if (full_connected)
                        {
                            // south-west
                            res (i - 1, j - 1, 0) = 1 * bw (i - 1, j - 1, 0);
                        }
                    }

            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = 1; j < dimy - 1; ++j)
                    if ( res (i, j, 0) == 1)
                    {
                        // west
                        res (i - 1, j, 0) = 1 * bw (i - 1, j, 0);
                        // north
                        res (i, j + 1, 0) = 1 * bw (i, j + 1, 0);
                        if (full_connected)
                        {
                            // north-west
                            res (i - 1, j + 1, 0) = 1 * bw (i - 1, j + 1, 0);
                        }
                    }

            for (uint i = 1; i < dimx - 1; ++i)
                for (uint j = dimy - 2; j > 0; --j)
                    if ( res (i, j, 0) == 1)
                    {
                        // east
                        res (i + 1, j, 0) = 1 * bw (i + 1, j, 0);
                        // south
                        res (i, j - 1, 0) = 1 * bw (i, j - 1, 0);
                        if (full_connected)
                        {
                            // south-east
                            res (i + 1, j - 1, 0) = 1 * bw (i + 1, j - 1, 0);
                        }
                    }

            norm1 = (res - resold).norm1();

            resold = res;

        }//end while

    }// end 2d case

    // 3d case
    else
    {
        while ( norm1 != 0 )
        {
            // up-north-east
            for (uint i = 1; i < dimx - 1; ++i)
                for (uint j = 1; j < dimy - 1; ++j)
                    for (uint k = 1; k < dimz - 1; ++k)
                        if ( res (i, j, k) == 1 )
                        {
                            // east
                            res (i + 1, j, k) = 1 * bw (i + 1, j, k);
                            // north
                            res (i, j + 1, k) = 1 * bw (i, j + 1, k);
                            // up
                            res (i, j, k + 1) = 1 * bw (i, j, k + 1);
                            if (full_connected)
                            {
                                // north-east
                                res (i + 1, j + 1, k) = 1 * bw (i + 1, j + 1, k);
                                // up-east
                                res (i + 1, j, k + 1) = 1 * bw (i + 1, j, k + 1);
                                // up-north
                                res (i, j + 1, k + 1) = 1 * bw (i, j + 1, k + 1);
                                // up-north-east
                                res (i + 1, j + 1, k + 1) = 1 * bw (i + 1, j + 1, k + 1);
                            }
                        }

            // up-north-west
            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = 1; j < dimy - 1; ++j)
                    for (uint k = 1; k < dimz - 1; ++k)
                        if ( res (i, j, k) == 1 )
                        {
                            // west
                            res (i - 1, j, k) = 1 * bw (i - 1, j, k);
                            // north
                            res (i, j + 1, k) = 1 * bw (i, j + 1, k);
                            // up
                            res (i, j, k + 1) = 1 * bw (i, j, k + 1);
                            if (full_connected)
                            {
                                // north-west
                                res (i - 1, j + 1, k) = 1 * bw (i - 1, j + 1, k);
                                // up-west
                                res (i - 1, j, k + 1) = 1 * bw (i - 1, j, k + 1);
                                // up-north
                                res (i, j + 1, k + 1) = 1 * bw (i, j + 1, k + 1);
                                // up-north-west
                                res (i - 1, j + 1, k + 1) = 1 * bw (i - 1, j + 1, k + 1);
                            }
                        }

            // up-south-west
            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = dimy - 2; j > 0; --j)
                    for (uint k = 1; k < dimz - 1; ++k)
                        if ( res (i, j, k) == 1 )
                        {
                            // west
                            res (i - 1, j, k) = 1 * bw (i - 1, j, k);
                            // south
                            res (i, j - 1, k) = 1 * bw (i, j - 1, k);
                            // up
                            res (i, j, k + 1) = 1 * bw (i, j, k + 1);
                            if (full_connected)
                            {
                                // south-west
                                res (i - 1, j - 1, k) = 1 * bw (i - 1, j - 1, k);
                                // up-west
                                res (i - 1, j, k + 1) = 1 * bw (i - 1, j, k + 1);
                                // up-south
                                res (i, j - 1, k + 1) = 1 * bw (i, j - 1, k + 1);
                                // up-south-west
                                res (i - 1, j - 1, k + 1) = 1 * bw (i - 1, j - 1, k + 1);
                            }
                        }

            // up-south-east
            for (uint i = 1; i < dimx - 1; ++i)
                for (uint j = dimy - 2; j > 0; --j)
                    for (uint k = 1; k < dimz - 1; ++k)
                        if ( res (i, j, k) == 1 )
                        {
                            // east
                            res (i + 1, j, k) = 1 * bw (i + 1, j, k);
                            // south
                            res (i, j - 1, k) = 1 * bw (i, j - 1, k);
                            // up
                            res (i, j, k + 1) = 1 * bw (i, j, k + 1);
                            if (full_connected)
                            {
                                // south-east
                                res (i + 1, j - 1, k) = 1 * bw (i + 1, j - 1, k);
                                // up-east
                                res (i + 1, j, k + 1) = 1 * bw (i + 1, j, k + 1);
                                // up-south
                                res (i, j - 1, k + 1) = 1 * bw (i, j - 1, k + 1);
                                // up-south-east
                                res (i + 1, j - 1, k + 1) = 1 * bw (i + 1, j - 1, k + 1);
                            }
                        }

            // down-north-east
            for (uint i = 1; i < dimx - 1; ++i)
                for (uint j = 1; j < dimy - 1; ++j)
                    for (uint k = dimz - 2; k > 0; --k)
                        if ( res (i, j, k) == 1 )
                        {
                            // east
                            res (i + 1, j, k) = 1 * bw (i + 1, j, k);
                            // north
                            res (i, j + 1, k) = 1 * bw (i, j + 1, k);
                            // down
                            res (i, j, k - 1) = 1 * bw (i, j, k - 1);
                            if (full_connected)
                            {
                                // north-east
                                res (i + 1, j + 1, k) = 1 * bw (i + 1, j + 1, k);
                                // down-east
                                res (i + 1, j, k - 1) = 1 * bw (i + 1, j, k - 1);
                                // down-north
                                res (i, j + 1, k - 1) = 1 * bw (i, j + 1, k - 1);
                                // down-north-east
                                res (i + 1, j + 1, k - 1) = 1 * bw (i + 1, j + 1, k - 1);
                            }
                        }

            // down-north-west
            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = 1; j < dimy - 1; ++j)
                    for (uint k = dimz - 2; k > 0; --k)
                        if ( res (i, j, k) == 1 )
                        {
                            // west
                            res (i - 1, j, k) = 1 * bw (i - 1, j, k);
                            // north
                            res (i, j + 1, k) = 1 * bw (i, j + 1, k);
                            // down
                            res (i, j, k - 1) = 1 * bw (i, j, k - 1);
                            if (full_connected)
                            {
                                // north-west
                                res (i - 1, j + 1, k) = 1 * bw (i - 1, j + 1, k);
                                // down-west
                                res (i - 1, j, k - 1) = 1 * bw (i - 1, j, k - 1);
                                // down-north
                                res (i, j + 1, k - 1) = 1 * bw (i, j + 1, k - 1);
                                // down-north-west
                                res (i - 1, j + 1, k - 1) = 1 * bw (i - 1, j + 1, k - 1);
                            }
                        }

            // down-south-west
            for (uint i = dimx - 2; i > 0; --i)
                for (uint j = dimy - 2; j > 0; --j)
                    for (uint k = dimz - 2; k > 0; --k)
                        if ( res (i, j, k) == 1 )
                        {
                            // west
                            res (i - 1, j, k) = 1 * bw (i - 1, j, k);
                            // south
                            res (i, j - 1, k) = 1 * bw (i, j - 1, k);
                            // down
                            res (i, j, k - 1) = 1 * bw (i, j, k - 1);
                            if (full_connected)
                            {
                                // south-west
                                res (i - 1, j - 1, k) = 1 * bw (i - 1, j - 1, k);
                                // down-west
                                res (i - 1, j, k - 1) = 1 * bw (i - 1, j, k - 1);
                                // down-south
                                res (i, j - 1, k - 1) = 1 * bw (i, j - 1, k - 1);
                                // down-south-west
                                res (i - 1, j - 1, k - 1) = 1 * bw (i - 1, j - 1, k - 1);
                            }
                        }

            // down-south-east
            for (uint i = 1; i < dimx - 1; ++i)
                for (uint j = dimy - 2; j > 0; --j)
                    for (uint k = dimz - 2; k > 0; --k)
                        if ( res (i, j, k) == 1 )
                        {
                            // east
                            res (i + 1, j, k) = 1 * bw (i + 1, j, k);
                            // south
                            res (i, j - 1, k) = 1 * bw (i, j - 1, k);
                            // down
                            res (i, j, k - 1) = 1 * bw (i, j, k - 1);
                            if (full_connected)
                            {
                                // south-east
                                res (i + 1, j - 1, k) = 1 * bw (i + 1, j - 1, k);
                                // down-east
                                res (i + 1, j, k - 1) = 1 * bw (i + 1, j, k - 1);
                                // down-south
                                res (i, j - 1, k - 1) = 1 * bw (i, j - 1, k - 1);
                                // down-south-east
                                res (i + 1, j - 1, k - 1) = 1 * bw (i + 1, j - 1, k - 1);
                            }
                        }

            norm1 = (res - resold).norm1();

            resold = res;

        }//end while

    }// end 3d case


    return;
}



template <typename T>
void im3d::image3d<T>::median_filter (image3d<T>& res, int const& radius) const
{
    // set proper dimension and spacing of the output
    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    int dimx = static_cast<int> (this->dimx);
    int dimy = static_cast<int> (this->dimy);
    int dimz = static_cast<int> (this->dimz);

    res = 0;
    // checking the radius
    if ( radius > dimx || radius > dimy || (radius > dimz && dimz != 1)
            || radius <= 0)
    {
        std::cout << "WARNING: in image3d::median_filter: radius is bigger than the smaller ";
        std::cout << "dimension or is negative" << std::endl;
        return;
    }

    // computing dimension of the mask
    int dimmask = (2 * radius + 1) * (2 * radius + 1);

    if (this->dimz == 1)
    {
        // initializing mask
        std::vector<T> mask;
        mask.reserve (dimmask);

        // internal nodes
        #pragma omp parallel for private (mask)
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)

                    {
                        mask.push_back ( (*this) (I + i, J + j, 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // 4 edges: north, east, south, west
        // north
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) (I + i, dimy - 1 - std::abs (dimy - (J + j) ), 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = radius; J < dimy - radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) (dimx - 1 - std::abs (dimx - (I + i) ), J + j, 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // south
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = 0; J < radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) (I + i, std::abs (J + j), 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // west
        for (int I = 0; I < radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) (std::abs (I + i), J + j, 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // 4 corners: north-west, north-east, south-east, south-west
        // north-west
        for (int I = 0; I < radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) (std::abs (I + i), dimy - 1 - std::abs (dimy - (J + j) ), 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // north-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                        mask.push_back ( (*this) (dimx - 1 - std::abs (dimx - (I + i) ),
                                                  dimy - 1 - std::abs (dimy - (J + j) ),
                                                  0) );

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // south-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = 0 ; J < radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) , std::abs (J + j) , 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

        // south-west
        for (int I = 0; I < radius; ++I)
            for (int J = 0 ; J < radius; ++J)
            {
                for (int i = -radius; i < radius + 1; ++i)
                    for (int j = -radius; j < radius + 1; ++j)
                    {
                        mask.push_back ( (*this) ( std::abs (I + i) , std::abs (J + j) , 0) );
                    }

                sort ( mask.begin() , mask.end() );
                res (I, J, 0) = mask[dimmask / 2];
                mask.resize (0);
            }

    }// end 2d version

    // 3d case
    else
    {
        dimmask *= (2 * radius + 1);

        // initialize mask
        std::vector<T> mask;
        mask.reserve (dimmask);

        // internal nodes
        #pragma omp parallel for private (mask)
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (I + i, J + j, K + k) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // 6 faces
        // down
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (I + i, J + j, std::abs (K + k) ) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // south
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (I + i, std::abs (J + j), K + k) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // west
        for (int I = 0; I < radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (std::abs (I + i), J + j, K + k) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // up
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) ( I + i, J + j, dimz - 1 - std::abs (dimz - (K + k) ) ) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // north
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (I + i, dimy - 1 - std::abs (dimy - (J + j) ), K + k) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                            {
                                mask.push_back ( (*this) (dimx - 1 - std::abs (dimx - (I + i) ), J + j, K + k) );
                            }

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // 12 adges

        // 4 down
        // south
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( I + i ,
                                                           std::abs (J + j) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // west
        for (int I = 0; I < radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           J + j ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // north
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( I + i ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           J + j ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // 4 up
        // south
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( I + i ,
                                                           std::abs (J + j) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // west
        for (int I = 0; I < radius; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           J + j ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // north
        for (int I = radius; I < dimx - radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( I + i ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = radius; J < dimy - radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           J + j ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // 4 center
        // south-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           std::abs (J + j) ,
                                                           K + k ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // south-west
        for (int I = 0; I < radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           K + k ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // north-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           K + k ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
        // north-west
        for (int I = 0; I < radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = radius; K < dimz - radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           K + k ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // 8 corners
        // down-south-west corner flip
        for (int I = 0; I < radius; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           std::abs (J + j) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // down-south-east side flip
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           std::abs (J + j) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // down-north-west side flip
        for (int I = 0; I < radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i),
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // down-north-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = 0; K < radius; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ),
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           std::abs (K + k) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }



        // up-south-west
        for (int I = 0; I < radius; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           std::abs (J + j) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // up-south-east
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = 0; J < radius; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           std::abs (J + j) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );

                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // up-nort-west
        for (int I = 0; I < radius; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( std::abs (I + i) ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );
                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }

        // up-north-east side flip
        for (int I = dimx - radius; I < dimx; ++I)
            for (int J = dimy - radius; J < dimy; ++J)
                for (int K = dimz - radius; K < dimz; ++K)
                {
                    for (int i = -radius; i < radius + 1; ++i)
                        for (int j = -radius; j < radius + 1; ++j)
                            for (int k = -radius; k < radius + 1; ++k)
                                mask.push_back ( (*this) ( dimx - 1 - std::abs (dimx - (I + i) ) ,
                                                           dimy - 1 - std::abs (dimy - (J + j) ) ,
                                                           dimz - 1 - std::abs (dimz - (K + k) ) ) );


                    sort ( mask.begin() , mask.end() );
                    res (I, J, K) = mask[dimmask / 2];
                    mask.resize (0);
                }
    }// end 3d version

    return;
}



template <typename T>
void im3d::image3d<T>::local_binary_pattern_edge_detector
(image3d<T>& res, T const& constant, double const& t1, double const& t2) const
{
    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    if ( this->dimx == 0 || this->dimy == 0 || this->dimz == 0 )
    {
        std::cout << "WARNING::local_binary_pattern_edge_detector: actual dimensions " <<
                  "generate an empty image" << std::endl;

        return;
    }

    T c (constant);

    if (c < 0)
    {
        c = std::abs (c);
    }

    res = 0;

    // initializing mask
    std::vector<bool> mask;

    int const r = 1, X (this->dimx - r), Y (this->dimy - r); //Z(this->dimz-r);

    // 3d case
    if ( this->dimz > 1 )
    {
        std::cout << "WARNING::local_binary_pattern_edge_detector: 3d version is not yet ";
        std::cout << "implemented" << std::endl;
        return;
    }

    // 2d case
    else
    {
        image3d<T> f (this->dimx, this->dimy, this->dimz), avg, std;
        mask.resize (8);

        avg = std = f;

        uint aux;

        // compute special avg
        #pragma omp parallel for
        for (int i = r; i < X; ++i)
            for (int j = r; j < Y; ++j)
            {
                avg (i, j, 0) += (*this) (i + r, j, 0);
                avg (i, j, 0) += (*this) (i + r, j - r, 0);
                avg (i, j, 0) += (*this) (i, j - r, 0);
                avg (i, j, 0) += (*this) (i - r, j - r, 0);
                avg (i, j, 0) += (*this) (i - r, j, 0);
                avg (i, j, 0) += (*this) (i - r, j + r, 0);
                avg (i, j, 0) += (*this) (i, j + r, 0);
                avg (i, j, 0) += (*this) (i + r, j + r, 0);
                avg (i, j, 0) /= 8;
            }

        // internal nodes
        for (int i = r; i < X; ++i)
            for (int j = r; j < Y; ++j)
            {
                mask[0] = (*this) (i + r, j, 0) > (*this) (i, j, 0) + c;
                mask[1] = (*this) (i + r, j - r, 0) > (*this) (i, j, 0) + c;
                mask[2] = (*this) (i, j - r, 0) > (*this) (i, j, 0) + c;
                mask[3] = (*this) (i - r, j - r, 0) > (*this) (i, j, 0) + c;
                mask[4] = (*this) (i - r, j, 0) > (*this) (i, j, 0) + c;
                mask[5] = (*this) (i - r, j + r, 0) > (*this) (i, j, 0) + c;
                mask[6] = (*this) (i, j + r, 0) > (*this) (i, j, 0) + c;
                mask[7] = (*this) (i + r, j + r, 0) > (*this) (i, j, 0) + c;

                aux = 0;
                if (mask[1] != mask[0])
                {
                    ++aux;
                }
                if (mask[2] != mask[1])
                {
                    ++aux;
                }
                if (mask[3] != mask[2])
                {
                    ++aux;
                }
                if (mask[4] != mask[3])
                {
                    ++aux;
                }
                if (mask[5] != mask[4])
                {
                    ++aux;
                }
                if (mask[6] != mask[5])
                {
                    ++aux;
                }
                if (mask[7] != mask[6])
                {
                    ++aux;
                }
                if (mask[7] != mask[0])
                {
                    ++aux;
                }


                // 1st STEP: LBP with a threshold (compute res)
                res (i, j, 0) = static_cast<T> (static_cast<int> (mask[0]) * 1 +
                                                static_cast<int> (mask[1]) * 2 +
                                                static_cast<int> (mask[2]) * 4 +
                                                static_cast<int> (mask[3]) * 8 +
                                                static_cast<int> (mask[4]) * 16 +
                                                static_cast<int> (mask[5]) * 32 +
                                                static_cast<int> (mask[6]) * 64 +
                                                static_cast<int> (mask[7]) * 128 );

                // 2nd STEP: compute auxiliar binary image to noise suppression (compute f)
                if ( aux == 2 && res (i, j, 0) != 1 && res (i, j, 0) != 2 && res (i, j, 0) != 4 &&
                        res (i, j, 0) != 8 && res (i, j, 0) != 16 && res (i, j, 0) != 32 &&
                        res (i, j, 0) != 64 && res (i, j, 0) != 128 && res (i, j, 0) != 127 &&
                        res (i, j, 0) != 191 && res (i, j, 0) != 223 && res (i, j, 0) != 239 &&
                        res (i, j, 0) != 247 && res (i, j, 0) != 251 && res (i, j, 0) != 253 &&
                        res (i, j, 0) != 254 )
                {
                    f (i, j, 0) = static_cast<T> (1);
                }

                // 3nd STEP: compute gradient as standard deviation * f
                std (i, j, 0) +=
                    ( (*this) (i + r, j, 0) - avg (i, j, 0) ) * ( (*this) (i + r, j, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i + r, j - r, 0) - avg (i, j, 0) ) * ( (*this) (i + r, j - r, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i, j - r, 0) - avg (i, j, 0) ) * ( (*this) (i, j - r, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i - r, j - r, 0) - avg (i, j, 0) ) * ( (*this) (i - r, j - r, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i - r, j, 0) - avg (i, j, 0) ) * ( (*this) (i - r, j, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i - r, j + r, 0) - avg (i, j, 0) ) * ( (*this) (i - r, j + r, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i, j + r, 0) - avg (i, j, 0) ) * ( (*this) (i, j + r, 0) - avg (i, j, 0) );
                std (i, j, 0) +=
                    ( (*this) (i + r, j + r, 0) - avg (i, j, 0) ) * ( (*this) (i + r, j + r, 0) - avg (i, j, 0) );
                std (i, j, 0) /= 8;
                std (i, j, 0) = sqrt ( std (i, j, 0) );

                std (i, j, 0) *= f (i, j, 0);
            }

        // 4nd STEP: edge-detection with two threshold
        std.im_to_black_and_white (avg, t2);

        res = 0;

        std.connected_component (res, avg, t1, true);

        std::cout << std::endl << "\tPixels survived to LBP class selection: ";
        std::cout << 100 * f.norm1() / (this->dimx * this->dimy) << " %" << std::endl;
        std::cout << "\tPixels survived to t2 threshold: ";
        std::cout << 100 * avg.norm1() / (this->dimx * this->dimy) << " %" << std::endl;
        std::cout << "\tPixels survived at the end: ";
        std::cout << 100 * res.norm1() / (this->dimx * this->dimy) << " %" << std::endl << std::endl;

    }// end 2d case

    return;
}

#endif // IMAGE3D_IMP_HPP_INCLUDED


