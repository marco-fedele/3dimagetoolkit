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
 \file image3d.hxx

 \brief header file of class \ref im3d::image3d
 */

#ifndef IMAGE3D_HPP_INCLUDED
#define IMAGE3D_HPP_INCLUDED

#include "names.hxx"

#include <algorithm>

namespace im3d
{

/*!
 \brief A template class conceived as a container of data and typical operations
 of a 3D image.

 This class is a container for data types that can be reassembled to a 3D image
 (mainly in size), making it possible to manipulate them, for instance in an
 algorithm.
 In particular it overloads arithmetic operators and implements various methods
 to access the data in both reading and writing mode.
 The main usage of this class is to represent a three dimensional image or an
 associate function defined on its uniform structured grid of pixels.
 It could be also used to represent a more general 3D signals, for instance a 3d audio
 signal or a video using the third dimension as time.
 It works good also with 2D images or signals if \ref image3d::dimz is set to 1.

 \tparam T is the type of data contained in the class, it can be any kind of number
 variable (int, float, double...) to represent black and white images.

 \warning This class is not conceived to represent RGB images. To make it work with
 RGB images user has to implement an appropriate son class reimplementing most of
 members.
 */
template <typename T = double>
class image3d
{

protected:
    //! a template vector containing the data of type T of the image.
    std::vector<T> rawimage;
    //! number of pixels in X direction.
    uint dimx;
    //! number of pixels in Y direction.
    uint dimy;
    //! number of pixels in Z direction.
    uint dimz;
    //! spacing between every pixel in X direction.
    double hx;
    //! spacing between every pixel in Y direction.
    double hy;
    //! spacing between every pixel in Z direction.
    double hz;

    /*!
     \brief private member used by both version of public member
     \ref connected_component.
     */
    template <typename S>
    void connected_component (image3d<S>& res, image3d<S>& bw,
                              bool const& full_connected) const;

    /*!
     \brief private member used by both version of public member
     \ref connected_component.
     */
    template <typename S>
    void cc_not_binary_output (image3d<S>& res, double threshold) const;

public:

    //CONSTRUCTORS

    //! Default Constructor to build an empty \ref image3d.
    image3d ();

    /*!
     \brief Variable Constructor to build a black (all zeros) \ref image3d of
     specified size and spacing.

     \param x set X dimension \ref dimx.
     \param y set Y dimension \ref dimy.
     \param z set Z direction \ref dimz.
     \param hx set X spacing \ref hx.
     \param hy set Y spacing \ref hy.
     \param hz set Z spacing \ref hz.
     */
    image3d (uint const& x, uint const& y, uint const& z,
             double const& hx = 1., double const& hy = 1., double const& hz = 1.);

    /*!
     \brief Copy Constructor.

     \param tocopy is the \ref image3d type object to be copied.
     */
    image3d (image3d const& tocopy);

    //!Default Denstructor
    virtual ~image3d() {};


    // MEMBERS TO GET AND SET PRIVATE PARAMETERS

    // Spacing

    //! Member to get X dimension \ref dimx.
    inline uint getdimx () const
    {
        return this->dimx;
    };

    //! Member to get Y dimension \ref dimy.
    inline uint getdimy () const
    {
        return this->dimy;
    };

    //! Member to get Z dimension \ref dimz.
    inline uint getdimz () const
    {
        return this->dimz;
    };

    /*!
     \brief Member to set X dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param x is the desired \ref dimx.
     */
    inline void setdimx (uint const& x);

    /*!
     \brief Member to set Y dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param y is the desired \ref dimy.
     */
    inline void setdimy (uint const& y);

    /*!
     \brief Member to set Z dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param z is the desired \ref dimz.
     */
    inline void setdimz (uint const& z);

    /*!
     \brief Member to set dimensions.

     This member also allows to choose the constant value which reinitializes the
     \ref image3d.

     \warning Changing dimensions reinitializes this \ref image3d as a constant
     image. To change dimensions without modify values of the image use \ref crop
     member.

     \param x is the desired \ref dimx.
     \param y is the desired \ref dimy.
     \param z is the desired \ref dimz.
     \param value is the constant desired value of the reinitialized image.
     */
    inline void setdim (uint const& x, uint const& y, uint const& z,
                        T const& value = 0);

    // Spacing

    //! Member to get X direction spacing \ref hx.
    inline double gethx () const
    {
        return this->hx;
    };

    //! Member to get Y direction spacing \ref hy.
    inline double gethy () const
    {
        return this->hy;
    };

    //! Member to get Z direction spacing \ref hz.
    inline double gethz () const
    {
        return this->hz;
    };

    /*!
     \brief Member to set X direction spacing.

     \param hx is the desired \ref hx.
     */
    inline void sethx (double const& hx);

    /*!
     \brief Member to set Y direction spacing.

     \param hy is the desired \ref hy.
     */
    inline void sethy (double const& hy);

    /*!
     \brief Member to set Z direction spacing.

     \param hz is the desired \ref hz.
     */
    inline void sethz (double const& hz);

    /*!
     \brief Member to set spacing.

     \param hx is the desired \ref hx.
     \param hy is the desired \ref hy.
     \param hz is the desired \ref hz.
     */
    inline void seth (double const& hx, double const& hy, double const& hz);


    //USEFULL OVERLOADING OF VAROIUS OPERATORS

    /*!
     \brief Operator() [const version] to quickly access to an image value with
     chosen coordinates.

     \param i is X coordinate of the desired value.
     \param j is Y coordinate of the desired value.
     \param k is Z coordinate of the desired value.
     \returns the T type scalar value
     */
    virtual inline T operator() (uint const& i, uint const& j, uint const& k) const;

    /*!
     \brief Operator() to quickly assign a value at desired coordinates of the image.

     \param i is X coordinate of the desired value.
     \param j is Y coordinate of the desired value.
     \param k is Z coordinate of the desired value.
     \returns the T type scalar value
     */
    virtual inline T& operator() (uint const& i, uint const& j, uint const& k);

    /*!
     \brief Overloading of = operator between a scalar and an \ref image3d.

     It is used to assign a scalar value to all coordinates of an \ref image3d with
     previously set dimensions.
     */
    template <typename S>
    image3d<T>& operator = (S const& toassign);

    /*!
     \brief Overloading of = operator between two \ref image3d.

     This operator copies all values of an \ref image3d after setting same dimensions
     and same spacing.
     */
    template <typename S>
    image3d<T>& operator = (image3d<S> const& toassign);

    /*!
     \brief Overloading of += operator between a scalar and an \ref image3d.

     This operator increases every element of current \ref image3d by the same scalar
     value.
     */
    template <typename S>
    image3d<T>& operator += (S const& addend);

    /*!
     \brief Overloading of += operator element by element of an \ref image3d.

     This operator sums every element with the same cooordinates of two \ref image3d.
     \warning It works only with two \ref image3d of the same dimensions.
     */
    template <typename S>
    image3d<T>& operator += (image3d<S> const& addend);

    /*!
     \brief Overloading of -= operator between a scalar and an \ref image3d.

     This operator decreases each element of current \ref image3d by the same scalar
     value.
     */
    template <typename S>
    image3d<T>& operator -= (S const& addend);

    /*!
     \brief Overloading of -= operator element by element of an \ref image3d.

     This operator subtracts each element with the same cooordinates of two
     \ref image3d.

     \warning It works only with two \ref image3d of the same dimensions.
     */
    template <typename S>
    image3d<T>& operator -= (image3d<S> const& addend);

    /*!
     \brief Overloading of *= operator between a scalar and an \ref image3d.

     This operator multiplies each element of current \ref image3d by the same scalar
     value.
     */
    template <typename S>
    image3d<T>& operator *= (S const& factor);

    /*!
     \brief Overloading of *= operator element by element of an \ref image3d.

     This operator multiplies each element with the same cooordinates of two
     \ref image3d.

     \warning It works only with two \ref image3d of the same dimensions.
     */
    template <typename S>
    image3d<T>& operator *= (image3d<S> const& factor);

    /*!
     \brief Overloading of /= operator between a scalar and an \ref image3d.

     This operator divides each element of current \ref image3d by the same scalar
     value.
     */
    template <typename S>
    image3d<T>& operator /= (S const& factor);

    /*!
     \brief Overloading of /= operator element by element of an \ref image3d.

     This operator divides each element with the same cooordinates of two
     \ref image3d.

     \warning It works only with two \ref image3d of the same dimensions.
     */
    template <typename S>
    image3d<T>& operator /= (image3d<S> const& factor);

    /*!
     \brief Overloading of + operator implementing sum between two images or an image
     and a scalar.

     It works also between two \ref image3d with different template parameters and it
     gives in output an \ref image3d with the same template parameter of the first
     addend
     */
    template <typename S, typename R>
    friend image3d<S> const operator + (image3d<S> const& addend1, R const& addend2);

    /*!
     \brief Overloading of - operator implementing subtraction between two images or
     an image and a scalar.

     It works also between two \ref image3d with different template parameters and it
     gives in output an \ref image3d with the same template parameter of the first
     addend
     */
    template <typename S, typename R>
    friend image3d<S> const operator - (image3d<S> const& addend1, R const& addend2);

    /*!
     \brief Overloading of * operator implementing multiplication between two images
     or an image and a scalar.

     It works also between two \ref image3d with different template parameter and it
     gives in output an \ref image3d with the same template parameter of the first
     factor
     */
    template <typename S, typename R>
    friend image3d<S> const operator * (image3d<S> const& factor1, R const& factor2);

    /*!
     \brief Overloading of / operator implementing division between two images or an
     image and a scalar.

     It works also between two \ref image3d with different template parameter and it
     gives in output an \ref image3d with the same template parameter of the first
     factor.
     */
    template <typename S, typename R>
    friend image3d<S> const operator / (image3d<S> const& factor1, R const& factor2);


    //MEMBERS WORKING ON VALUES OF THE IMAGE

    /*!
     \brief Friend function that computes scalar product between image3d

     This member computes scalar product considering an \ref image3d like a matrix,
     so it doesn't consider spacing between pixels/voxels.

     \param factor1 is the first factor of the product
     \param factor2 is the second factor of the product
     \return result
     */
    template <typename S, typename R>
    friend S const scalarprod (image3d<S> const& factor1, image3d<R> const& factor2);

    /*!
     \brief Friend function that computes L2 scalar product between image3d

     This member computes scalar product L2 considering \ref image3d like a constant
     function in each pixel/voxel.

     \warning The two factors have to be not only of the same dimensions, but also of
     the same spacing to make this function works properly.

     \param factor1 is the first factor of the product
     \param factor2 is the second factor of the product
     \return result
     */
    template <typename S, typename R>
    friend S const scalarprod_L2 (image3d<S> const& factor1, image3d<R> const& factor2);

    /*!
     \brief Member that computes gradient of an \ref image3d

     \param res is a 2D/3D vector depending on dimensions of the original image.
     */
    template <typename S>
    void grad (std::vector<image3d<S> >& res) const;

    /*!
     \brief Friend function that computes modulus of a 2d or 3d vector of \ref image3d
     */
    template <typename S, typename R>
    friend void vector_abs (image3d<S>& res, std::vector<image3d<R> > const& fun);

    /*!
     \brief Friend function that computes divergence of an \ref image3d
     */
    template <typename S, typename R>
    friend void div (image3d<S>& res, std::vector<image3d<R> > const& fun);

    /*!
     \brief Member that computes norm 1 of an \ref image3d

     This member computes norm 1 considering an \ref image3d like a matrix, so it
     doesn't consider spacing between pixels/voxels.

     \return the norm value
     */
    T const norm1 () const;

    /*!
     \brief Member that computes norm L1 of an \ref image3d

     This member computes norm L1 considering \ref image3d like a constant function
     in each pixel/voxel.

     \return the norm value
     */
    T const normL1 () const;

    /*!
     \brief Member that computes norm 2 of an \ref image3d

     This member computes norm 2 considering an \ref image3d like a matrix, so it
     doesn't consider spacing between pixels/voxels.

     \return the norm value
     */
    T const norm2 () const;

    /*!
     \brief Member that computes norm L2 of an \ref image3d

     This member computes norm L2 considering \ref image3d like a constant function
     in each pixel/voxel.

     \return the norm value
     */
    T const normL2 () const;

    /*!
     \brief Member that computes norm infinite of an \ref image3d

     \return the norm value
     */
    T const norminf () const;

    /*!
     \brief Member that computes maximum value of an \ref image3d
     \return the max value
     */
    T const max () const;

    /*!
     \brief Member that computes minimum value of an \ref image3d
     \return the min value
     */
    T const min () const;


    // PROCESSING

    /*!
     * \name Processing
     * \brief Public member to modify an image useful for pre/post processing
     */
    ///@{
    /*!
     \brief Member to crop an \ref image3d in a smaller one.

     Member to obtain in output a smaller image made of the same values of the
     original one stored in a selected range of pixels.

     \param res is the output where it's written the crop image
     \param xstart is the X direction index in which we start the crop
     \param xend is the X direction index in which we finish the crop
     \param ystart is the Y direction index in which we start the crop
     \param yend is the Y direction index in which we finish the crop
     \param zstart is the Z direction index in which we start the crop
     \param zend is the Z direction index in which we finish the crop
     */
    void crop (image3d<T>& res,
               uint const& xstart, uint const& ystart, uint const& zstart,
               uint const& xend, uint const& yend, uint const& zend) const;

    /*!
     \brief Member to crop an \ref image3d in a smaller one.

     Member to obtain in output a smaller image made of a selected fraction of
     the original one.

     \param res is the output where it's written the crop image
     \param xstart is the fraction of the X dimension of the image in which we start
     the crop.
     \param xend is the fraction of the X dimension of the image in which we finish
     the crop
     \param ystart is the fraction of the Y dimension in which we start the crop
     \param yend is the fraction of the Y dimension of the image in which we finish
     the crop
     \param zstart is the fraction of the Z dimension in which we start the crop
     \param zend is the fraction of the Z dimension of the image in which we finish
     the crop
     */
    void crop (image3d<T>& res,
               double const& xstart, double const& ystart, double const& zstart,
               double const& xend, double const& yend, double const& zend) const;

    /*!
     \brief Member to change resolution of an \ref image3d.

     This member gives in output an \ref image3d with a smaller/bigger resolution.
     If boolean variable increase
     is false, it generates a smaller image selecting only some pixel of
     the original one.
     Otherwise, if increase is true it generates a bigger image making a linear
     interpolation
     of the value of the pixel of the original image.
     \warning If you want to increase resolution, ratio could be only a power of 2.
     If it isn't, this member automatically change ratio value to nearest bigger
     power of 2.

     \param res is the output image
     \param ratio represents the desired ratio between original image and output image
     \param increase is true to generate a bigger image in output, otherwise it's false
     */
    void change_resolution (image3d<T>& res, uint ratio = 2, bool increase = false) const;

    /*!
     \brief Member that gives in output an image with an equalized histogram of
     intensity values.

     It allows user to choose quantization.

     \param res is the output image.
     \param quantization is the desired quantization of the output image
     (the number of values that intensity could assume).
     */
    void histogram_equalization (image3d<T>& res, uint const& quantization) const;

    /*!
     \brief Member that gives in output an image with an equalized histogram of
     intensity values.

     It chooses as quantization the difference between maximum and minimum value of
     the image intensity.

     \param res is the output image.
     */
    void histogram_equalization (image3d<T>& res) const;

    /*!
     \brief Member that gives in output a black & white version of the image.

     This member gives in output a b&w image using a percentage threshold to decide
     which pixels should be white and which black. First of all it calculates minimum
     and maximum value of intensity, then it assumes as black all pixels between
     minimum value and min+(max-min)*threshold and as white others.

     \param res is the output image.
     \param threshold indicates the percentage threshold used to divide pixels into
     black & white and must be between 0 and 1.
     \param negative if is true gives in output the negative image of the b&w image
     (default false)
     */
    template<typename S>
    void im_to_black_and_white (image3d<S>& res,
                                double threshold = 0.5, bool negative = false) const;

    /*!
     \brief Member that modifies image changing range of intensity between a min value
     and a max value selected from user.

     This member pratically rescales every pixel value in order to set minimal value
     and maximal value of the image to desired values chosen by user.
     If user chooses only maximal value, minimal value is automatically set to zero.

     \param max is desired maximal value of the new image.
     \param min is desired minimal value of the new image.
     */
    void change_range_of_intensity (T const& max, T const& min = 0);

    /*!
     \brief Member that gives in output an image with only pixels of input image that
     has a selected range of intensity.

     This member selects only pixels with intensity between lowerbound and upperbound.
     Value of other pixels depends on value of variable type.

     \param res is the output image.
     \param lowerbound is the desired minimun intensity value of output image.
     \param upperbound is the desired maximum intensity value of output image.
     \param type sets the way to substitute values out of range. If it is negative
     all values out of desired range are set to lowerbound, if it is positive they are
     set to upperbound, if it is zero (default) value smaller than lowerbound are set
     to lowerbound and value bigger than upperbound are set to upperbound
     */
    void select_range_of_intensity (image3d<T>& res, T const& lowerbound,
                                    T const& upperbound, int type = 0, T lowervalue = 0, T uppervalue = 0) const;

    /*!
     \brief Member that gives in output a black & white image with only a connected
     component.

     First of all it transforms the image in a b&w image using chosen threshold
     (the default value 0.5 is suggested),
     then it extracts the connected component of the image containing the input
     coordinates. To choose desired coordinates it could be useful
     to use member \ref interface::get_coordinates.

     \param res is the output image.
     \param i is the pixel index of the first dimension coordinate.
     \param j is the pixel index of the second dimension coordinate.
     \param k is the pixel index of the third dimension coordinate.
     \param threshold indicates the percentage threshold used to divide pixels
     into black & white and must be between 0 and 1.
     \param full_connected is a boolean value (default false) that, if it's set true,
     considers connected also pixels sharing only corners. Otherwise pixels are
     considered connected only if they share an edge in 2d or a face in 3d.
     \param binary_output if false (default true) the method return an image equal to the
     original image everywhere exept in the not connected white part of the black and white image.
     In this part the image is set to its minimum.
     */
    template <typename S>
    void connected_component (image3d<S>& res,
                              uint const& i, uint const& j, uint const& k, double threshold = 0.5,
                              bool full_connected = false, bool binary_output = true) const;

    /*!
     \brief Member that gives in output a black & white image with all connected
     components linked to an initial image.

     First of all it transforms the image in a b&w image using chosen threshold
     (the default value 0.5 is suggested), then it extracts all pixels connected with
     white pixels of input image init (every pixels not equal to 0 will be consired
     white from this algorithm).

     \param res is the output image.
     \param init is the input boolean image and has to be of the same dimensions.
     \param threshold indicates the percentage threshold used to divide pixels
     into black
     & white and must be between 0 and 1.
     \param full_connected is a boolean value (default false) that, if it's set true,
     considers connected also pixels sharing only corners. Otherwise pixels are
     considered connected only if they share an edge in 2d or a face in 3d.
     \param binary_output if false (default true) the method return an image equal to the
     original image everywhere exept in the not connected white part of the black and white image.
     In this part the image is set to its minimum.
     */
    template <typename S>
    void connected_component (image3d<S>& res,
                              image3d<S> const& init, double threshold = 0.5,
                              bool full_connected = false, bool binary_output = true) const;

    /*!
     \brief Member that gives in output the image filtered with a median filter

     \param res is the output image.
     \param radius is the dimension in pixel of the filter in every direction
    */
    void median_filter (image3d<T>& res, int const& radius = 1) const;

    /*!
     \brief Member that gives in output the image filtered with the local binary
     pattern algorithm for edge detection [ADD BIB]

     \warning This member works only with 2D images

     \param res is the output image
     \param constant is the constant to compute Local Binary Pattern in a
     circumference of pixels
     \param t1 is the lower threshold executes at the end of the algorithm
     \param t2 is the threshold executes after gradient computation with local binary
     pattern
     */
    void local_binary_pattern_edge_detector (image3d<T>& res,
                                             T const& constant = 0,
                                             double const& t1 = 0.25,
                                             double const& t2 = 0.5) const;

    ///@}

};// end of class image3d



// FRIEND MEMBER DECLARATION

template <typename S, typename R>
image3d<S> const operator + (image3d<S> const& addend1, R const& addend2);

template <typename S, typename R>
image3d<S> const operator - (image3d<S> const& addend1, R const& addend2);

template <typename S, typename R>
image3d<S> const operator * (image3d<S> const& factor1, R const& factor2);

template <typename S, typename R>
image3d<S> const operator / (image3d<S> const& factor1, R const& factor2);

template <typename S, typename R>
S const scalarprod ( image3d<S> const& factor1, image3d<R> const& factor2 );

template <typename S, typename R>
S const scalarprod_L2 ( image3d<S> const& factor1, image3d<R> const& factor2 );

template <typename S, typename R>
void div (image3d<S>& res, std::vector<image3d<R> > const& fun);

template <typename S, typename R>
void vector_abs (image3d<S>& res, std::vector<image3d<R> > const& fun);


}//end namespace im3d


#include "image3d_imp.hxx"


#endif // IMAGE3D_HPP_INCLUDED


