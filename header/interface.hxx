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
 \file  interface.hxx

 \brief header file of class \ref im3d::interface
 */

#ifndef INTERFACE_HXX_INCLUDED
#define INTERFACE_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"
#include "abstractinterface.hxx"

#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkImageViewer.h>
#include <vtkMetaImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkJPEGReader.h>
#include <vtkBMPReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkMarchingCubes.h>
#include <vtkContourFilter.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkCamera.h>
#include <vtkImageCast.h>



/*!
 * \brief A namespace in which all objects linked to
 * an image 3d are defined.
 */
namespace im3d
{

// STATIC VARIABLES
static vtkSmartPointer<vtkRenderer> RENDERER = NULL;
static vtkSmartPointer<vtkRenderWindowInteractor> RENWININT = NULL;

/*!
 * \brief A class inheriting from \ref abstractinterface conceived as a linker
 * between an \ref image3d and an \ref imageptr of type vtkSmartPointer<vtkImageData>.
 *
 * This class follows the guideline of \ref abstractinterface providing the user with
 * an interface class specialized in communication between an \ref image3d and
 * elements from the vtk library.
 * The vtk library is an open source project for 3D computer graphics, see www.vtk.org
 * for more details.
 * The use of this library makes the visualization process easy and aesthetically
 * satisfying, with planes the allow to explore the image as well as the possibility
 * to freely rotate it in the 3-dimensional space.
 * To make the program work without vtk libraries the user has to implement another
 * child class using an alternative way to implement abstract members.
 */
template <typename S>
class interface : public abstractinterface < vtkSmartPointer<vtkImageData>, S >
{

protected:

    //! member to define the background color of the visualization.
    double colour[3];
    //! member to define the opacity of the image during visualization.
    double opacity;

public:

    // CONSTRUCTORS

    /*!
     \brief default constructor
     */
    interface ();

    /*!
     \brief constructor from an image name.

     Will construct the class automatically calling the load function with the desired
     filename.

     \param imagename a string containing the path to the image to be processed,
     including the file's extension.
     Supported extensions are .mhd, .dcm, .jpg and .bmp.
     */
    interface (std::string const& imagename);

    /*!
     \brief copy constructor
     */
    interface (image3d<S> const& myim);

    /*!
     \brief default destructor
     */
    virtual ~interface () {};

    // CONVERSION

    /*!
     \brief function to convert an image3d into a T image
     */
    void convertfromimage3d (image3d<S> const& myim);

    /*!
     \brief function to convert T image into an image3d
     */
    void convert2image3d (image3d<S>& myim) const;

    // SHOWING

    /*!
     \brief member to set background color.

     the color is specified in RGB values.
     \param red the Red value.
     \param green the Green value.
     \param blue the Blue value.
     */
    void setcolour (double const& red = 1., double const& green = 1., double const& blue = 0.)
    {
        colour[0] = red;
        colour[1] = green;
        colour[2] = blue;
        return;
    };

    /*!
     \brief member to set opacity

     \param opacity should be a value between 0 and 1.
     */
    inline void setopacity (double const& opacity = 0.5)
    {
        this->opacity = opacity;
        return;
    };

    /*!
     \brief member to show only a level of a T \ref imageptr.

     This member will produce an image where only the pixels with a value equal to
     level are highlighted. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations.

     \param level the value of the level to be shown.
     */
    void show_contour (double const& level = 0);

    /*!
     \brief Member to show an image using three orthogonal planes.

     This member will show the image in a box made of three movable planes.
     These planes may be directed using the mouse, by pressing both buttons at the
     same time and moving it. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations and the addition of planes.
     */
    void show_image ();

    /*!
     \brief Member to show an image using three orthogonal planes and a selected level.

     This function is a combination of \ref show_contour
     and \ref show_image. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations and the addition of planes.
     \param level the value of the level to be shown.
     */
    void show_image_and_contour (double const& level = 0);

    /*!
     \brief Member to show a contour of an image with another one as background.

     This function is just like \ref show_image_and_contour except that the image from
     which the level is extracted and the background image need not coincide.
     An extensive use of the vtk library is made, concerning in particular rendering
     and/or shading operations and the addition of planes.

     \param level the value of the level to be shown.
     \param backgroundimage a reference to an interface class linked to the background
     image.
     */
    void show_contour_with_background_image (interface& backgroundimage,
                                             double const& level = 0);

    // COORDINATES

    /*!
     \brief Member to get coordinates of a pixel.

     This function helps user to find coordinates of a pixel showing an image using
     \ref show_image and asking user for coordinates.
     It returns real coordinates considering spacing between pixels. The coordinates
     returned are the same that user visualizes during execution of the member.

     \param x is the pixel index of the first dimension coordinate.
     \param y is the pixel index of the second dimension coordinate.
     \param z is the pixel index of the third dimension coordinate.
     */
    void get_coordinates (double& x, double& y, double& z);

    /*!
     \brief Member to get coordinates of a pixel.

     This function helps user to find coordinates of a pixel showing an image using
     \ref show_image and asking user for coordinates.
     It returns coordinates in pixel (integer coordinates).

     \param i is the pixel index of the first dimension coordinate.
     \param j is the pixel index of the second dimension coordinate.
     \param k is the pixel index of the third dimension coordinate.
     */
    void get_coordinates (uint& i, uint& j, uint& k);

    /*!
     \brief Member to get coordinates of a pair of pixels.

     This function helps user to find coordinates of a pair of pixels showing an image
     using \ref show_image and asking user for coordinates.
     It returns real coordinates considering spacing between pixels. The coordinates
     returned are the same that user visualizes during execution of the member.

     \param x1 is the first pixel index of the first dimension coordinate.
     \param y1 is the first pixel index of the second dimension coordinate.
     \param z1 is the first pixel index of the third dimension coordinate.
     \param x2 is the second pixel index of the first dimension coordinate.
     \param y2 is the second pixel index of the second dimension coordinate.
     \param z2 is the second pixel index of the third dimension coordinate.
     */
    void get_coordinates (double& x1, double& y1, double& z1,
                          double& x2, double& y2, double& z2);

    /*!
     \brief Member to get coordinates of a pair of pixels.

     This function helps user to find coordinates of a pair of pixels showing an image
     using \ref show_image and asking user for coordinates.
     It returns coordinates in pixel (integer coordinates).

     \param i1 is the first pixel index of the first dimension coordinate.
     \param j1 is the first pixel index of the second dimension coordinate.
     \param k1 is the first pixel index of the third dimension coordinate.
     \param i2 is the second pixel index of the first dimension coordinate.
     \param j2 is the second pixel index of the second dimension coordinate.
     \param k2 is the second pixel index of the third dimension coordinate.
     */
    void get_coordinates (uint& i1, uint& j1, uint& k1,
                          uint& i2, uint& j2, uint& k2);

    // WRITING

    /*!
     \brief function to write a T image in <imagename>.<extension>

     This function will write a file containing all the information relative to the
     image linked to the class.
     \param imagename the name of the file.
     \param extension the extension of the file.
     */
    void write (std::string const& imagename = "image3d",
                std::string const& extension = ".mhd") const;

    // SUPPORT FUNCTIONS FOR VTK USAGE

    /*!
     \brief Member to add a contour to the renderer
     */
    void addcontour2renderer (vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                              double const& level = 0);

    /*!
     \brief Member to add a contour to the renderer
     */
    void addcontour2renderer (double const& level,
                              vtkSmartPointer<vtkRenderer>& renderer = RENDERER);

    /*!
     \brief Member to add planes to the image view
     */
    void addplanes2renderwindowinteractor
    (vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
     vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT) const;

    /*!
     \brief Member to create a new rendering window
     */
    static void newrendering (vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                              vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT);

    /*!
     \brief Member to delete a rendering window
     */
    static void deleterendering ( vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                                  vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT);

    /*!
     \brief Member to reset a rendering window
     */
    static void reset (vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                       vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT);

    /*!
     \brief Member to initiate the rendering process

     This member will actually begin showing the image.
     */
    static void startshowing (vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                              vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT,
                              int const& windowwidth = 800, int const& windowheight = 480);

    /*!
     \brief Member to initiate the rendering process

     This member will actually begin showing the image.
     */
    static void startshowing (int const& windowwidth, int const& windowheight,
                              vtkSmartPointer<vtkRenderer>& renderer = RENDERER,
                              vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt = RENWININT);


};//end of class interface


}//end namespace im3d

#include "interface_imp.hxx"


#endif // INTERFACE_HXX_INCLUDED
