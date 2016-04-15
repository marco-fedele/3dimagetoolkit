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
 \file convolution.hxx

 \brief file with functor \ref conv::filtering and the linked functions
 */

#ifndef CONVOLUTION_HXX_INCLUDED
#define CONVOLUTION_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"

#include <fftw3.h>

/*!
 * \brief A namespace grouping all objects or functions linked to filtering an
 * \ref im3d::image3d using convolution.
 */
namespace conv
{

/*!
 \brief A template functor conceived to filter an \ref im3d::image3d with a previously
 chosen filter using convolution.

 This class implements various ways to do various kinds of convolution between an
 \ref im3d::image3d that the user can pass directly to the functor and an assigned
 filter saved in static members of the class in both time and frequency version.\n
 First of all, the user has to decide what kind of filter he'd like to build choosing
 from those available (using the related static function). All implemented filters are
 real and even.\n
 Then he can construct a functor choosing the kind of convolution desired. He can
 choose between aperiodic discrete convolution (both in frequency or time domain) and
 periodic discrete convolution (in frequency domain).
 Convolution in frequency domain uses fast
 fourier transform algorithm, thus is more efficient if desired filter is of relatively
 big dimension (see more details on function documentation
 \ref acyclic_fftw_convolution).
 The user can also implement by himself his own version of convolution
 (exploiting static filters) passing it to the class constructor.\n
 Some example of usage of class \ref filtering are shown in file
 \ref image_toolkit.cxx.

 \warning It's important to remember that before using functor on a new
 \ref im3d::image3d with different dimensions than the first one, the user has to
 rebuild filters (or at least
 frequency filters, if he doesn't want to change time filter).
 Check \ref build_freq_filters static functions for more details.

 \note This functor could filter also bidimensional images, in this case the user has
 to pass to the functor an \ref im3d::image3d with \ref im3d::image3d::dimz equal to 1.

 \tparam T is the same template parameter of the \ref im3d::image3d to filter and it
 characterizes also the filters built inside the class. It could be any kind of
 floating number.
 Integer number (for instance int or unsigned int) are not allowed, hence if you
 want to filter an \ref im3d::image3d with T=int you have to transform it
 (before passing it to the functor) in an \ref im3d::image3d with T equal to a
 floating number type using im3d::image3d::operator = ().
 */
template<typename T>
class filtering
{

private:

    //! Private member storing a function pointer to the chosen convolution function.
    void (*f) (im3d::image3d<T>& result, im3d::image3d<T> const& factor);

public:

    // CONSTRUCTOR AND FUNCTOR LINKER MEMBER

    /*!
     \brief Constructor taking a function pointer to assign to the private member
     function pointer the chosen convolution function.
     */
    filtering (void (*f) (im3d::image3d<T>& result, im3d::image3d<T> const& factor) )
        : f (f) {};

    /*!
     \brief Operator() implemented to make class \ref filtering a functor.

     If you call it, class \ref filtering will automatically call the function
     associated to the object when it was built or after the last call of member
     \ref changef.

     \param result is the first argument of the private function where result of
     chosen convolution is written.
     \param factor is the second argument of the private function where the factor
     is stored.
     */
    void operator () (im3d::image3d<T>& result, im3d::image3d<T> const& factor) const
    {
        this->f (result, factor);
    };

    /*!
     \brief This member allows user to change private function associated to the
     functor.

     After changing the private function associated to the functor using this member,
     operator()( ) will call the new function every time you use it.

     \param newf is the new function you are setting and it has to take same
     parameters of the private function (two \ref im3d::image3d, the second const).
     */
    void changef ( void (*newf) (im3d::image3d<T>& result, im3d::image3d<T> const& factor) )
    {
        this->f = newf;
    };

    //FILTER DEFINITION

    /*!
     \brief Static member to store the time version of the filter.

     The kind of filter stored in this member depends on what related static
     function has been called to build it. However, every kind of \ref time_filter
     has got some common properties:\n
     -# It's real and even and therefore its number of elements is odd.
     This ensures that also its fourier transform is real and even.\n
     -# The middle position in each direction represents the zero of the cartesian
     axis. Hence its maximum value is usually (not necessarily) in the middle position.

     Since this member is public, the user can directly initialize it as he likes.
     He has only to respect the two proprerties just shown to ensure a proper
     functioning of the functor.
     */
    static im3d::image3d<T> time_filter;

    /*!
     \brief Static member to store the frequency version of the filter used to perform
     cyclic convolution.

     In this member the fourier transform of the time filter is stored.
     Since fft considers the value in the first position on each direction as the zero
     of the cartesian axes, before performing fft, the negative part of time filter
     is shifted at the end of a local vector.\n
     Dimensions of the filter are equal to the dimensions of the convolution factor or
     to an integer multiple of it if \ref time_filter is bigger than factor.
     (Check \ref cyclic_fftw_convolution documentation for more details.)
     Therefore, between the positive and negative parts of the time filter are added
     as many zeros as needed to reach these dimensions.
     All these actions are performed by \ref build_freq_filters function.
     \warning This member is public since expert users could choose to initialize
     directly a \ref cycl_freq_filter (implementing an own function) without
     previously inizialize a \ref time_filter and without calling
     \ref build_freq_filters function. This approach
     is suggested only if you know what is the correct way of building a filter in a
     frequency domain to perform a cyclic convolution. Otherwise it's better to use
     this member "as a private member", calling \ref build_freq_filters function to
     build it after \ref time_filter is already built.
     */
    static im3d::image3d<T> cycl_freq_filter;

    /*!
     \brief Static member to store the frequency version of the filter used to perform
     acyclic convolution.

     In this member the fourier transform of the time filter is stored. Since fft
     considers the value in the first position on each direction as the zero of the
     cartesian axes, before performing fft, the negative part of the time filter is
     shifted at the end of a local vector. Furthermore, between the positive and
     negative parts of the time filter are added as many zeros as needed to reach a
     minimal dimension that ensures a correct performance of the aperiodic
     convolution, thus avoiding aliasing. All these actions are performed by
     \ref build_freq_filters function.
     \warning This member is public since expert users could choose
     to initialize directly a \ref acycl_freq_filter (implementing an own function)
     without previously initializing a \ref time_filter and without calling
     \ref build_freq_filters function. This approach
     is suggested only if you know what is the correct way to build a filter in a
     frequency domain to perform an acyclic convolution. Otherwise it's better to use
     this member "as a private member", calling \ref build_freq_filters function to
     build it after \ref time_filter is already built.
     */
    static im3d::image3d<T> acycl_freq_filter;

    /*!
     \brief Static member to build static filters as a Gaussian kernels.

     This member builds the \ref time_filter as a Gaussian kernel and automatically
     calls the \ref build_freq_filters member to build the corrisponding
     \ref cycl_freq_filter and \ref acycl_freq_filter (note that the Fourier transform
     of a Gaussian kernel is still a Gaussian kernel with a different sigma).
     Have a look at \ref image_toolkit.cxx file for an example of usage.

     \param f is the factor you'd like to filter and is passed to this function only
     to get its size and spacing.
     \param sigma is the standard deviation of the time domain Gaussian kernel
     (in pixel).
     \param radius is a parameter to indicate the distance, as multiples of sigma,
     in which the Gaussian kernel could be considered null.
     Default value is equal to 3. If you want to improve accuracy make it greater,
     if you want to decrease computational costs make it smaller (not too small!).
     \param pixel_approach is a boolean variable that changes way of building the
     gaussian kernel. If it is set to false (default) this member considers spacing
     between pixels to build kernel, hence it builds a \ref time_filter with bigger
     dimensions if spacing of factor f are smaller
     (because you need more pixels to obtain a kernel of a certain radius if spacing
     are smaller) and viceversa. On the contrary, if it
     is set to true it doesn't care about spacing of the factor and it builds a
     \ref time_filter considering spacing equal to 1 in each direction.
     Hence, this alternative way allows user to make a gaussian blur based on number
     of pixels.
     \warning this member builds a discrete version of a continous Gaussian kernel,
     for a good approximation of the continous one is suggested to use sigma greater
     than 3*h (where h is the smallest spacing between \ref im3d::image3d::hx,
     \ref im3d::image3d::hy and \ref im3d::image3d::hz of the factor to convolve).
     By the way, to ensure that convolving with this continuous kernel is a
     "good mean operation", norm1 of \ref filtering::time_filter is manually set to
     one.
     */
    static void build_gaussian_filter (im3d::image3d<T> const& f,
                                       T const& sigma,
                                       double const& radius = 3,
                                       bool pixel_approach = false);

    /*!
     \brief Static member to build \ref time_filter as an high-pass filter.

     This member builds the \ref time_filter as an high-pass filter and automatically
     calls the \ref build_freq_filters member to build the
     corrisponding \ref cycl_freq_filter
     and \ref acycl_freq_filter. Practically it builds \ref time_filter as a sort of
     laplacian mask, but with the central value greater than the sum of others to make
     average value unitary.

     \param f is the factor you'd like to filter and is passed to this function only
     to get its size and spacing.
     \param central_node_coeff is a coefficient that the more is greater than one the
     more the central node has an higher value compared to the others (default 1).
     \param pixel_approach is a boolean variable that changes way of building the
     high-pass filter. If it is set to false (default) this member considers spacing
     between pixels to build kernel giving higher coefficients in directions with
     higher spacing (supposing spacing less than 1, otherwise the opposite).
     On the contrary, if it is set to true it doesn't care about spacing of the factor
     and it builds a \ref time_filter considering spacing equal to 1 in each direction.
     */
    static void build_high_pass_filter (im3d::image3d<T> const& f,
                                        double const& central_node_coeff = 1,
                                        bool pixel_approach = false);

    /*!
     \brief Static member to build \ref time_filter as a laplacian filter.

     This member builds the \ref time_filter as a laplacian mask and automatically
     calls the \ref build_freq_filters member to build the corrisponding
     \ref cycl_freq_filter and \ref acycl_freq_filter.

     \param f is the factor you'd like to filter and is passed to this function only
     to get its size and spacing.
     \param pixel_approach is a boolean variable that changes way of building the
     high-pass filter. If it is set to false (default) this member considers spacing
     between pixels to build kernel giving higher coefficients in directions with
     higher spacing. On the contrary, if it is set to true it doesn't care about
     spacing of the factor and it builds a \ref time_filter considering spacing
     equal to 1 in each direction.
     */
    static void build_laplacian_filter (im3d::image3d<T> const& f,
                                        bool pixel_approach = false);

    /*!
     \brief Static member to build both cyclic and acyclic frequency filters, starting
     from time filter.

     This member builds \ref acycl_freq_filter and \ref cycl_freq_filter starting from
     \ref time_filter and giving them the correct dimensions and structure to perform
     a convolution in the frequency domain using in a correct way the convolution
     theorem. Check \ref acycl_freq_filter and \ref cycl_freq_filter documentation for
     more details.
     \param f is the factor you'd like to filter and is passed to this function only
     to get its size.
     \warning If you perform a convolution in the frequency domain with a factor f and
     then you want to perform this kind of convolution again but with a new factor of
     different dimensions, you have to call this function another time before
     performing it again.
     If you don't, both \ref cycl_freq_filter and \ref acycl_freq_filter
     cuold be of incorrect dimensions and hence \ref cyclic_fftw_convolution or
     \ref acyclic_fftw_convolution could not works.
     On the contrary, if the new factor has got the same dimensions as the previous
     one it's not necessary to call again this function
     (In this case it's suggested not to call it, because it's useless and
     computationally expensive). Have a look at the
     exemple shown on \ref acyclic_fftw_convolution documentation for more details.
     */
    static void build_freq_filters (im3d::image3d<T> const& f);

};// end functor



// INITIALIZING STATIC MEMBERS

template <typename T>
im3d::image3d<T> filtering<T>::time_filter;

template <typename T>
im3d::image3d<T> filtering<T>::cycl_freq_filter;

template <typename T>
im3d::image3d<T> filtering<T>::acycl_freq_filter;



// CONVOLUTION FUNCTIONS

/*!
 \brief This function performs an acyclic convolution between the previously
 initialized \ref filtering::acycl_freq_filter and a selected \ref im3d::image3d.

 It could be used outside of the functor \ref filtering or it could be passed directly
 to the functor using its member \ref filtering::changef
 or its constructor as it's shown
 in the following example. Note that when you call a static function to build static
 filters (for instance \ref filtering::build_gaussian_filter) and you want to perform
 a frequency convolution, it's necessary to give also the factor to convolve in order
 to compute correct dimensions of \ref filtering::cycl_freq_filter. So if you want to
 perform a new convolution with a new factor of different dimensions you have to
 update dimensions of frequency filters before doing it, calling
 \ref filtering::build_freq_filters static member.

 \code
    // initializing image3d
    im3d::im3d::image3d<double> cycl_res, acycl_res, f1(100,100,100,1,1,1), f2(50,30,45);

    ... importing datas on f1 and f2 ...

    // initializing filtering functor with cyclic convolution function
    conv::filtering<double> test_conv (conv::cyclic_fftw_convolution);

    // building gaussian filter and releted frequency filter of correct dimensions
    // to perform convolution with f1
    conv::filtering<double>::build_gaussian_filter (f1,5,2);

    // performing cyclic convolution with f1
    test_conv (cycl_res, f1);
    // now the result of convolution with f1 is stored on cycl_res

    // updating frequency filters to perform convolution with f2
    conv::filtering<double>::build_freq_filters (f2);

    // change function linked to the functor to perform an acyclic convolution
    test_conv.changef(conv::acyclic_fftw_convolution);

    // performing acyclic convolution with f2
    test_conv (acycl_res, f2);
    // now the result of convolution with f2 is stored on acycl_res
 \endcode

 \param result is the output of acyclic convolution
 \param factor is the factor of convolution to convolve with
 \ref filtering::acycl_freq_filter
 */
template <typename T>
void acyclic_fftw_convolution (im3d::image3d<T>& result, im3d::image3d<T> const& factor);

/*!
 \brief This function perform a cyclic convolution between the previously inizialized
 \ref filtering::cycl_freq_filter and a selcted \ref im3d::image3d.

 This function could be used outside of the functor \ref filtering or it could be
 passed directly to the functor using its member
 \ref filtering::changef or its constructor.
 Note that when you call a static function to build static
 filters (for instance \ref filtering::build_gaussian_filter) and you want to perform
 a frequency convolution it's necessary to give also the factor to convolve in order
 to compute correct dimensions of \ref filtering::cycl_freq_filter. So if you want to
 perform a new convolution with a new factor of different dimensions you have to
 update dimensions of frequency filters before doing it calling
 \ref filtering::build_freq_filters static member. Have a look at the exemple shown
 on \ref acyclic_fftw_convolution documentation.

 \param result is the output of cyclic convolution
 \param factor is the factor of convolution to convolve with
 \ref filtering::cycl_freq_filter
 */
template <typename T>
void cyclic_fftw_convolution (im3d::image3d<T>& result, im3d::image3d<T> const& factor);

/*!
 \brief This function perform a direct acyclic convolution between the previously
 inizialized \ref filtering::time_filter and a selcted \ref im3d::image3d.

 It could be used outside of the functor \ref filtering or it could be passed directly
 to the functor using its member \ref filtering::changef or its constructor.
 Before using it, you have to inizialze \ref filtering::time_filter using a chosen
 static function of the class \ref filtering (for instance
 \ref filtering::build_gaussian_filter). In this case it's not important to select the
 factor of convolution before performing it since \ref filtering::time_filter has got
 indipendent dimensions from the factor. Have a look at the exemple shown
 on \ref acyclic_fftw_convolution documentation to familiarize with the usage
 of this family of functions linked to various kinds of convolution.
 \param result is the output of direct acyclic convolution
 \param factor is the factor of convolution to convolve with
 \ref filtering::time_filter
 */
template <typename T>
void direct_convolution (im3d::image3d<T>& result, im3d::image3d<T> const& factor);

}// end namespace conv

#include "convolution_imp.hxx"

#endif // CONVOLUTION_HXX_INCLUDED
