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
 \file   segmentation.hxx

   \brief header file of class \ref segm::segmentation
 */

#ifndef SEGMENTATION_HXX_INCLUDED
#define SEGMENTATION_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"
#include "interface.hxx"

using namespace std;

/*!
 \brief A namespace in which are stored all classes or function related to a segmentation
 algorithm.
 */
namespace segm
{

/*!
 \brief A template class conceived as a container for the fundamental tools and
 properties related to the segmentation of an image using a level set method.

 This class implements some tools needed to manipulate an image during a segmentation
 algorithm. In particular it contains members that allow to visualize the image as
 well as members that are able to perform the segmentation algorithm and all the
 related computations and finally also the level set function itself.
 This class does not implement a particular algorithm on purpose, leaving the user
 the liberty to exploit it in complete freedom.
 Any particular algorithm should be implemented in a child
 class, e.g. \ref rsfe_splitbregman.

 \tparam T should solely be used to set the precision of the pixel values: float,
 double...
 */
template <typename T>
class segmentation
{

protected:

    /*!
     \brief a number representing the level of the levelset in which contour is
     expected to be.
     */
    T alpha;

    /*!
     \brief member to handle the visualization process. It represents the image to be
     segmented.
     */
    im3d::interface<T> image;

    /*!
     \brief Member to handle the visualization process of the level set function.
     It is linked to \ref phi
     */
    im3d::interface<T> levelset;

    /*!
     \brief Member representing the level set function. This member should evolve
     during the course of the algorithm.
     */
    im3d::image3d<T> phi;

public:
    //CONSTRUCTOR

    //! Default constructor
    segmentation() {};

    //! Member to set private attribute \ref alpha
    virtual inline void set_alpha (T const& alpha)
    {
        this->alpha = alpha;
    };

    /*!
     * \brief Pure virtual member to apply algorithm starting from an initial image
     *
     * \param myim is the image to segment
     * \param init is the initial value with whom levelset \ref phi is initialized
     */
    virtual void apply (im3d::image3d<T> const& myim, im3d::image3d<T> const& init) = 0;

    /*!
     * \brief Pure virtual member to apply algorithm
     *
     * \param myim is the image to segment
     */
    virtual void apply (im3d::image3d<T> const& myim) = 0;

    //! Member to get output of the algorithm (\ref phi)
    inline im3d::image3d<T> getoutput () const
    {
        return this->phi;
    }

    //! Member to show level \ref alpha of the levelset \ref phi
    virtual inline void show_contour ();

    //! Member to show levelset \ref phi with planes
    virtual inline void show_levelset ();

    //! Member to show levelset \ref phi with planes and with its level \ref alpha
    virtual inline void show_levelset_and_contour();

    /*!
     * \brief Member to show level \ref alpha of the levelset \ref phi with the image
     * to segment as background
     */
    virtual inline void show_levelset_with_background_image();

}; //close class segmentation

} // close namespace segm


// IMPLEMENTATIONS

template <typename T>
void segm::segmentation<T>::show_contour()
{
    levelset.convertfromimage3d (this->phi);
    levelset.show_contour (this->alpha);
    return;
}

template <typename T>
void segm::segmentation<T>::show_levelset()
{
    levelset.convertfromimage3d (this->phi);
    levelset.show_image();
    return;
}

template <typename T>
void segm::segmentation<T>::show_levelset_and_contour()
{
    levelset.convertfromimage3d (this->phi);
    levelset.show_image_and_contour (this->alpha);
    return;
}

template <typename T>
void segm::segmentation<T>::show_levelset_with_background_image()
{
    levelset.convertfromimage3d (this->phi);
    levelset.show_contour_with_background_image (this->image, this->alpha);
    return;
}


#endif // SEGMENTATION_HXX_INCLUDED

