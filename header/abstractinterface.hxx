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
 \file abstractinterface.hxx

 \brief file with class \ref im3d::abstractinterface
 */

#ifndef ABSTRACTINTERFACE_HXX_INCLUDED
#define ABSTRACTINTERFACE_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"


/*!
 * \brief A namespace in which all objects linked to an image 3d are defined.
 */
namespace im3d
{
/*!
 * \brief A pure virtual class conceived as a linker between an \ref image3d
 * and a template \ref imageptr usefull to show, load and save an image.
 *
 * This pure virtual class is conceived as a sort of guide for the progammer.
 * He has to implement all of the abstract members of the class to make the
 * library work. The goal is to have an interface class (the son class)
 * specialized in communication between an \ref image3d and a template \ref
 * imageptr. This template \ref imageptr has to be in a special format in which is
 * easier to show a 3d image on the screen or write it in a file. For all
 * this reasons the son class \ref interface is implemented using
 * vtk libraries. To make the library work without vtk libraries user has to
 * implement itself another son class using an alternative way to implement
 * all abstract members of this class.
 */
template <typename T, typename S>
class abstractinterface
{
protected:
    /*!
     * \brief A template variable containing the type chosen to represent the image.
     *
     * We suggest to implement the child class in this way: \n
     * class \ref interface : public \ref abstractinterface < * "imagetype" > \n
     * where "imagetype" is the type choosed to represent the image.
     * In this way T becomes a pointer to an "imagetype".
     */
    T imageptr;

public:
    /*!
     \brief default constructor
     */
    abstractinterface () {};

    /*!
     \brief constructor that load an image from a file.

     \param imagename is name of the file
     */
    abstractinterface (std::string const& imagename) {};

    /*!
     \brief constructor that loads an image from an \ref im3d::image3d

     \param myim is the image3d to link with this abstractinterface
     */
    abstractinterface (image3d<S> const& myim) {};

    /*!
     \brief default destructor
     */
    virtual ~abstractinterface () {};

    /*!
     \brief Member to get the \ref imageptr allowing reading only.
     */
    T const& getimageptr()
    {
        return this->imageptr;
    };

    /*!
     * \brief Member to convert an \ref im3d::image3d into a T
     * \ref im3d::abstractinterface::imageptr.
     *
     * \param myim is the input image3d
     */
    virtual void convertfromimage3d (image3d<S> const& myim) = 0;

    /*!
     * \brief Member to convert a T \ref im3d::abstractinterface::imageptr into
     * an \ref im3d::image3d.
     *
     * \param myim is the output image3d
     */
    virtual void convert2image3d (image3d<S>& myim) const = 0;

    /*!
     * \brief Member to show only a level of a T \ref imageptr.
     *
     * \param level is the level of the image to show
     */
    virtual void show_contour (double const& level = 0) = 0;

    /*!
     * \brief Member to show a T \ref imageptr using three orthogonal planes.
     */
    virtual void show_image () = 0;

    /*!
     * \brief Member to show a T \ref imageptr using three orthogonal plane and a
     * selected level.
     *
     * \param level is the level of the image to show
     */
    virtual void show_image_and_contour (double const& level = 0) = 0;

    /*!
     * \brief Member to show a level of a T \ref imageptr with another image as a
     * background.
     *
     * \param backgroundimage is another abstractinterface storing the background
     * image
     * \param level is the level of the image to show
     */
    virtual void show_contour_with_background_image
    (abstractinterface& backgroundimage, double const& level = 0) {};

    /*!
     * \brief Member to get coordinates of a pixel after showing image.
     *
     * \param i is the pixel index of the first dimension coordinate
     * \param j is the pixel index of the second dimension coordinate
     * \param k is the pixel index of the third dimension coordinate
     */
    virtual void get_coordinates (uint& i, uint& j, uint& k) = 0;

    /*!
     * \brief Member to get coordinates of a pair of pixels after showing image.
     *
     * \param i1 is the first pixel index of the first dimension coordinate
     * \param j1 is the first pixel index of the second dimension coordinate
     * \param k1 is the first pixel index of the third dimension coordinate
     * \param i2 is the second pixel index of the first dimension coordinate
     * \param j2 is the second pixel index of the second dimension coordinate
     * \param k2 is the second pixel index of the third dimension coordinate
     */
    virtual void get_coordinates (uint& i1, uint& j1, uint& k1,
                                  uint& i2, uint& j2, uint& k2) = 0;

    /*!
     * \brief Member to write an image in a file.
     *
     * \param imagename contains the name of the file
     *  \param extension contains the extension of the file including dot (e.g. ".mhd")
     */
    virtual void write (std::string const& imagename, std::string const& extension) const = 0;

}; // end class abstractinterface

}// end namespace im3d

#endif // ABSTRACTINTERFACE_HXX_INCLUDED
