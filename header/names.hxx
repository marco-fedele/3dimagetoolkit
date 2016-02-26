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
 \file names.hxx
 
 \brief file with some include and typedef used by the whole library
 */

#ifndef NAMES_HPP_INCLUDED
#define NAMES_HPP_INCLUDED

// include files used by all the others in the library
#include<iostream>
#include<vector>
#include<cmath>
#include<omp.h>




/*!
 * \brief definition of uint as unsigned int
 * Made only to guarantee the compatibility with those compiler without uint defined
 */
typedef unsigned int uint;

#endif // NAMES_HPP_INCLUDED

/**
 \mainpage

 \author  Marco Fedele, fedele.marco@gmail.com
 \author  Luca Barbarotta, luca.barbarotta@gmail.com
 \author  Francesco Cremonesi, francesco.cremonesi0@gmail.com
 \author  Elena Faggiano, elena.faggiano@gmail.com
 
 \copyright Copyright 2013 Marco Fedele, Luca Barbarotta, Francesco Cremonesi, Elena Faggiano. All rights reserved. Read details of the license in the file license.txt provided with the library.

 <B>Download Link</B>: under construction!


 3DImageToolkit is a C++ library focused on medical image processing and segmentation.
 The main goal of the project is to implement a segmentation algorithm, useful to separate
 objects inside a 3d image, working good with medical images in wich usually the strong
 inhomogeneity and the proximity of the various organs makes it difficult to obtain good
 results.

 The library consists of a main class called \ref im3d::image3d (to store an image and to
 do some basic operations between two images) and others classes to do specifical
 operations on an image. We can group the features of the library into three families:
 visualization, pre/post-processing and segmentation.

 \section Visualization
 We create a class called \ref im3d::interface to read, visualize or write in a file a 3d
 image.\n
 First of all this class can read a file and transform it into an \ref im3d::image3d.
 It's allowed to work with files of the following extensions:
 - .mhd/.mha, a Meta Image header linked witk a .raw file in which datas are stored.
 - DICOM files (.dcm), a very popular extension used for CAT scan, MRI or other medical
 images. In this case you can select a single .dcm file (to read only a 2d slice of
 your datas) or the full directory with all slices transforming it into a three
 dimensional image.
 - .jpg and .bmp are supported too, because the whole library works well also with
 standard bidimensional images.

 Secondly it can write an \ref im3d::image3d into a .mhd-.raw file.\n
 Finally it can visualize a 3d image with movable planes or showing a single intensity
 level allowing user customizations on colour, opacity, background, etc...\n
 This class is built using some useful classes of VTK open-source library, please have a
 look at this link for further information: http://www.vtk.org .

 \section Processing
 In this library there are many ways to manipulate an image before segmentation algorithm.
 The main goal is to obtain an image with the required property to make the algorithm
 works better or faster.
 It could be interesting also modify the output image of the algorithm, for instance
 extracting only an object from it.\n
 We have implemented, as members of class \ref im3d::image3d, methods to crop an image,
 change its resolution, equalize its histogram, apply a median filter,
 extract a connected component and so on.\n
 We have also created a class called \ref conv::filtering to filter images with a
 pre-specified filter using convolution.
 This allows user to make, for instance, gaussian blur of an image or to apply
 a laplacian filter. To make convolution faster when dealing with high dimensions filter,
 we use frequency domain approach using FFTW free library to obtain a fast fourier
 transform of images.
 Please have a look at this link for further information: http://www.fftw.org .\n
 All this features can be easily applied to an image using executable generated from file
 \ref image_toolkit.cxx . Here are some pictures that show some of the effects quickly
 described above.

\image html preprocessing.png "A preprocessing example: original (left), range of intensity values (center), median filter (right)." width=10px

 \section Segmentation
 We implemented a general abstract class called \ref segm::segmentation from which
 all other classes linked to a specific algorithm derive.\n
 In this first release we only implement a region-based algorithm with a level set
 approach consisting on the merging of a Region Scalable Fitting Energy Method and the
 Split Bregman algorithm in order to obtain a robust and efficient technique suited to
 segment images with inhomogeneous intensity (ADD REFERENCES TO THE ARTICLE).
 We also customized the algorithm allowing user to extract at any time
 a connected component of partial result. In this way it is possible to avoid the
 connection between nearby objects inside an image and to obtain a better result of
 the part of interest of your image.
 We organize this algorithm in \ref segm::rsfe_splitbregman class allowing user to
 easily set all parameters of the algorithm before applying it, for further information
 have a look to class documentation.

\image html tac_result_zoom.png "A result of segmentation at different iterations" width=20cm

 */
