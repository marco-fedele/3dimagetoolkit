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
 \file rsfe_splitbregman.hxx

  \brief header file of class \ref segm::rsfe_splitbregman
 */

#ifndef RSFE_SPLITBREGMAN_HXX_INCLUDED
#define RSFE_SPLITBREGMAN_HXX_INCLUDED

#include "image3d.hxx"
#include "interface.hxx"
#include "segmentation.hxx"
#include "poisson.hxx"
#include "convolution.hxx"
#include "exceptions.hxx"

#include <GetPot>
#include <sstream>


namespace segm
{

/*!
 \brief A class to apply a segmentation algorithm consisting on the application of
 the Split-Bregman method to minimize the Region Scalable Fitting Energy Functional.

 This algorithm divides an \ref im3d::image3d into two regions, hence it works good to
 segment images with two clearly defined set of intensity values.
 It works good also with inhomogeneous images thanks to its property of approximation
 of the intensity of the image to segment using convolution with a Gaussian kernel.\n
 More precisely this algorithm minimize the following functional, where \ref phi is
 the levelset function in which the contour of the image is stored and \f$u_0\f$ is
 the intensity function of the image:
 \f[
 \mathcal{F}(\phi(\mathbf{x}))
 \quad=\quad
    \nu~\left\|\nabla \phi \right\|_{g} \quad+\quad
    \left\langle \phi, \lambda_{1}e_{1} - \lambda_{2}e_{2} \right\rangle_{L^2(\Omega)}
 \f]
 \f[
 \mbox{where}\qquad e_i(\mathbf{x})
 \quad=\quad
    \int_{\Omega}
        K_\sigma(\mathbf{x}-\mathbf{y})\left| u_0 (\mathbf{x}) -
        f_i(\mathbf{y})\right|^2
    d \mathbf{y},
 \f]
 \f[
 \mbox{with}\qquad
 K_\sigma(\mathbf{x})=
\frac{1}{(\sqrt{2\pi}\ \sigma)^d}\ \exp \left\{ -\frac{|\mathbf{x}|^2}{2\sigma^2}\right\},
\quad \mathbf{x} \in \mathbf{R}^d,
\qquad \qquad
 f_{i}(\mathbf{x}) = \frac{K_{\sigma} \ast \Big[  M_{i}^{\varepsilon}
 \big(\phi(\mathbf{x})\big) u_{0}(\mathbf{x}) \Big] }
 {K_{\sigma}(\mathbf{x}) \ast M_{i}^{\varepsilon}\big(\phi(\mathbf{x})\big) },
 \f]
 \f[
 M_{1}^{\varepsilon}(\mathbf{x}) = \frac{1}{2}
 \left(
    1 + \frac{2}{\pi}\arctan\left(
    \frac{\phi(\mathbf{x})-(a_0 +\alpha)}{\varepsilon|b_0 -a_0|}
    \right)
 \right),
 \qquad \qquad
 M_{2}^{\varepsilon}(\mathbf{x}) = 1 - M_{1}^{\varepsilon}.
 \f]
 \f[
 \|f\|_{g}~=~
    \left\|~g(u_0)~f~\right\|_{L^1(\Omega)}, \qquad
g\big(u_0(\mathbf{x})\big)  ~=~
    \frac{1}{1+\beta\left|edge\_detector\big(u_0(\mathbf{x})\big)\right|^2}.
 \f]
 The behaviour of the algorithm changes with the modification of the parameters
 characterizing the functional. Using the default values is a good starting point, but
 it could be necessary to modify the parameters to obtain a better segmentation.
 To understand better how to modify each parameter see the documentation of each one:
 \ref sigma, \ref lamda, \ref lamda1, \ref lamda2, \ref nu, \ref beta, \ref dt,
 \ref ls_steps, \ref epsilon, \ref a0, \ref b0, \ref alpha.\n
 Another feature of the algorithm is to 'think globally' dividing the whole image into
 two regions. On the contrary, to focus only on a part of the image, you can use
 the automatic extraction of a connected component of the image setting some specific
 parameters (\ref auto_extract_conn_comp, \ref auto_extract_conn_comp_freq).\n
 Also other customization are allowed, see the documentation of all private variables
 for more details.
*/
template <typename T>
class rsfe_splitbregman : public segmentation<T>
{

private:

    /*!
     * \brief A type to define various kind of boundary condition
     */
    enum bc_type {NEU = 1, DIR};

    /*!
     * \brief A type to define various kind of edge detector
     */
    enum ed_type {GRAD = 1, GRAD_THRESHOLD, LBP};


    // GENERAL PARAMETERS
    /*!
     * \name General Parameters
     * \brief Variables related to numbers of iteration, frequency of showing
     * partial results and general setting
     */
    ///@{
    /*!
     * \brief Name of ".pot" file used by \ref apply to set some private attribute via
     * GetPot (default="data")
     *
     * This variable, if it is set to a value different from the default empty string,
     * allows \ref apply to call private member \ref set_param_from_getpot before
     * starting
     * iterations of the algorithm. In this way algorithm initializes all variables in
     * section [\\init] of file getpotfile and it updates all variables in section
     * [\\onthego]
     * every \ref showfrequency iterations if private variable \ref onthego is true.
     *
     * \note If the desired ".pot" file is not in your current directory you can
     * set this string with the path to reach it.
     *
     * \warning It is useless to set \ref onthego equal to true if this variable is
     * set equal to an empty string or to a string that doesn't match with any
     * existing file
     */
    std::string getpotfile;
    /*!
     * \brief It is true to produce a verbose output during execution of \ref apply,
     * otherwise false (default=true)
     */
    bool verbosity;
    /*!
     * \brief It indicates the maximum number of iterations allowed by the algorithm
     * (default=100)
     */
    uint maxiter;
    /*!
     * \brief Number that indicates how often \ref apply shows partial result
     * (0=never, default=0)
     */
    uint showfrequency;
    /*!
     * \brief Number that indicates how often \ref apply saves partial result
     * (0=never, default=0)
     */
    uint dumpfrequency;
    /*!
     * \brief Name of desired output file in which partial and final results are stored
     * (default="result")
     *
     * In case of not null \ref dumpfrequency, \ref apply will save partial results in
     * a .mhd-.raw file called <outputname>_<iteration number> adding also the string
     * "_final" to the last result. In case of null \ref dumpfrequency it will not
     * save any files.
     */
    std::string outputname;
    /*!
     * \brief Name of the logfile in which history of the algorithm is saved
     * (default="")
     *
     * If this string is not null, stderr is redirected to a text file called
     * logfilename in append modality.
     * In this way all important informations of the execution of of the algorithm
     * like initial private attributes, norm of errors at each iteration,
     * update of private attributes or extraction of a connected component of the
     * levelset are saved in a file for consultation. On the contrary, if this string
     * is null, the program will write all this informations on the screen.
     *
     * \note The algorithm execution will produce a verbose output only in case of a
     * true value of private attribute \ref verbosity
     */
    std::string logfilename;
    /*!
     * \brief If it is true (default), it allows user to update some private attributes
     * from file \ref getpotfile (section [\\onthego]) every \ref showfrequency
     * iterations.
     */
    bool onthego;
    /*!
     * \brief If is true, the levelset of the current iteration is saved in file
     * \ref outputname (default false)
     *
     * \note The only way to set this variable to true is through section [\\onthego]
     * of file \ref getpotfile, so it works only if \ref onthego is true and
     * \ref getpotfile links to an existing ".pot" file. After saving current levelset
     * its value is set to false again.
     */
    bool save_current;
    /*!
     * \brief Variable settable via \ref getpotfile (section [\\onthego]) to end
     * algorithm at the current iteration without reaching \ref maxiter
     *
     * \warning Settable only from \ref getpotfile during execution, if \ref onthego
     * is true.
     */
    bool end_now;

    ///@}
    // RSFE PARAMETERS
    /*!
     * \name RSFE Parameters
     * \brief All the parameters that appear in the Region Scalable Fitting Energy
     * Functional
     */
    ///@{
    /*!
     * \brief set the width of all the Gaussian kernels (default=\f$3\f$)
     *
     * A quite general rule is to set \f$\sigma\f$ an order of magnitude smaller than
     * the diameter of the object of interest in the image. A too high \f$\sigma\f$
     * will cause an aggregation between all near objects inside the image.
     * A \f$\sigma\f$ a bit bigger than the one of the rule speeds up the algorithm
     * in the very first iterations and avoids the problem of the creation of other
     * contours inside the one of interest. Hence, sometimes a good compromise could
     * be to adopt an hybrid strategy starting from a bigger \f$\sigma\f$ for a
     * few iterations and after decreasing it on-the-go.
     */
    T sigma;

    /*!
     * \brief it deals with the splitting of the functional with the Split-Bregman
     * method (default=\f$5 \times 10^{-2}\f$)
     *
     * It has an effect on the speed of movement of the contour. The greater it is,
     * the slower the algorithm becomes. Despite this seems to prove that a smaller
     * value is always better, the property of slowing down the algorithm increasing
     * \f$\lambda\f$ can be really useful in lots of practical applications.
     */
    T lamda;

    /*!
     * \brief \f$\lambda_1\f$ and \f$\lambda_2\f$ weigh the terms \f$e_1\f$ and
     * \f$e_2\f$ one in relation to the other (default=\f$10^{-5}\f$)
     *
     * It is suggested to set \f$\lambda_1 = \lambda_2 \f$. Indeed, setting
     * \f$\lambda_2 > \lambda_1\f$ will encourage the contour to develop towards the
     * outside of the object of interest overestimating it a bit also at the
     * convergence, while the opposite strategy underestimates the contour especially
     * during its evolution. The larger the difference between the two values is,
     * the more this behaviour is evident.
     * However, sometimes setting these parameters in order
     * to underestimate the contour during its evolution can help to avoid the problem
     * of aggregation between near objects.
     */
    T lamda1;

    /*!
     * \brief \f$\lambda_1\f$ and \f$\lambda_2\f$ weigh the terms \f$e_1\f$ and
     * \f$e_2\f$ one in relation to the other (default=\f$10^{-5}\f$)
     *
     * It is suggested to set \f$\lambda_1 = \lambda_2 \f$. Indeed, setting
     * \f$\lambda_2 > \lambda_1\f$ will encourage the contour to develop towards the
     * outside of the object of interest overestimating it a bit also at the
     * convergence, while the opposite strategy underestimates the contour especially
     * during its evolution. The larger the difference between the two values is,
     * the more this behaviour is evident.
     * However, sometimes setting these parameters in order
     * to underestimate the contour during its evolution can help to avoid the problem
     * of aggregation between near objects.
     */
    T lamda2;

    /*!
     * \brief it weighs the edge detector at the denominator of the function
     * \f$g\big(u_0(\mathbf{x})\big)\f$ (default = \f$100\f$)
     *
     * A general rule is to increase it if the edge detector used works good with
     * your image. Otherwise it is better to decrease it.
     */
    T beta;

    /*!
     * \brief desired time step of the time discretization (default=\f$1\f$).
     *
     * Useful to slow down the algorithm in some cases.
     */
    T dt;

    /*!
     * \brief The number of steps that \ref apply does at each iteration to solve
     * linear system (default=\f$5\f$)
     *
     * At each iteration \ref segmentation::phi becomes the solution of the unsteady
     * Poisson problem at \f$ t{n+1} = t_n + \f$\ref ls_steps \ref dt.
     * Hence this variable tells the algorithm how many steps of the resolution of
     * linear system it has to execute.
     * This also depends on the method chosen to solve it and on its convergence
     * properties. Generally the more steps you set, the more accurate is the solution
     * (but you need more time to solve it).
     * It is also suggested to use this parameter to slow down the algorithm when
     * necessary.
     */
    uint ls_steps;

    /*!
     * \brief the shape parameter that determines how sharp is the approximation of
     * the Heaviside (default=\f$10^{-3}\f$)
     *
     * A general rule is to set \f$\varepsilon\f$ smaller than two orders of magnitude
     * in respect of the image spacing.
     */
    T epsilon;

    /*!
     * \brief It is the lower bound that levelset \ref segmentation::phi can assume
     * (default=\f$0\f$)
     *
     * This is the minimal allowed value for \ref segmentation::phi. In fact, at
     * the end of each iteration of the algorithm, every lower value is forced to be
     * equal to \ref a0 in order to ensure that values of the levelset are in the range
     * [\ref a0, \ref b0].\n
     * In case of Dirichlet \ref bc , this is also the value
     * of the boundary condition of the Poisson problem.
     *
     * \warning Cannot set this attributes greater than the value of \ref b0
     */
    T a0;

    /*!
     * \brief It is the upper bound that levelset \ref segmentation::phi can assume
     * (default=\f$1\f$)
     *
     * This is the maximal allowed value for \ref segmentation::phi. In fact, at
     * the end of each iteration of the algorithm, every bigger value is forced to be
     * equal to \ref b0 in order to ensure that values of the levelset are in the range
     * [\ref a0, \ref b0].
     *
     * \warning Cannot set this attributes lower than the value of \ref a0
     */
    T b0;

    ///@}
    // CONNECTED COMPONENT
    /*!
     * \name Connected Component
     * \brief Variables related to the extraction of a connected component of the
     * partial result
     */
    ///@{
    /*!
     * \brief If it is set to true (default=false), automatic extraction of a connected
     * component is activate every \ref auto_extract_conn_comp_freq
     *
     * This variable, when it is true, allows \ref apply to extract a connected
     * component of the current levelset \ref segmentation::phi calling
     * public member \ref im3d::image3d::connected_component on the levelset
     * only if a selected pixel is white (a pixel is considered white if its value is
     * between \ref segmentation::alpha and \ref b0). The coordinates of this pixel
     * are settable by the user and they are stored in private attributes
     * \ref cc_init_pixel_x, \ref cc_init_pixel_y and \ref cc_init_pixel_z.
     * After extraction, the connected component is automatically set as current
     * levelset and the algorithm restarts from this new initial
     * \ref im3d::image3d.\n
     * By the way, this is not the only way to extract a connected component of
     * \ref segmentation::phi. If this variable is false, every \ref showfrequency
     * iterations, after the show of partial results, user can decide to extract it
     * answering a question. Only if he likes the output, he could decide to set it
     * as current levelset.\n
     * The difference between the two approaches is that the first is totally automatic
     * and it doesn't need user interaction with the show of partial results (it works
     * also if \ref showfrequency is 0), the second needs an 'active' user to works.
     */
    bool auto_extract_conn_comp;
    //! X coordinate of the starting pixel of the automatic extraction
    double cc_init_pixel_x;
    //! Y coordinate of the starting pixel of the automatic extraction
    double cc_init_pixel_y;
    //! Z coordinate of the starting pixel of the automatic extraction
    double cc_init_pixel_z;

    bool cc_binary_output;

    bool cc_init_variables;

    /*! Number that indicates, if \ref auto_extract_conn_comp is true, how often
     * \ref apply automatically extracts a connected_component
     */
    uint auto_extract_conn_comp_freq;

    ///@}
    // EXPERT
    /*!
     * \name Expert
     * \brief Parameters to set the edge detector that appears in the second addend
     * of the RSFE functional
     */
    ///@{
    /*!
     * \brief The kind of boundary conditions to give to the Poisson problem
     * depends on this value (default=NEU)
     *
     * If it is equal to NEU the problem is solved with Neumann condition and this
     * allows contours with holes on the boundary. If it is set to DIR it is solved
     * with Dirichlet condition imposing \ref a0 at the boundary and hence preventing
     * holes in it.
     */
    bc_type bc;
    /*!
     * \brief The tolerance of the chosen linear system method (default=0, not set)
     *
     * This variable by default is not used by the algorithm because it uses variable
     * \ref ls_steps to choose when finishing the iterations of the linear system
     * method.
     * Instead, if user set a value different from zero,
     * the iterative method to solve linear system continue also after \ref ls_steps
     * iterations till the
     * percentage L2-norm error is greater than this value. It is suggested not to use
     * a value too small, for instance a value of 0.01 it is more than enough.
     */
    T ls_tol;
    /*!
     * \brief The kind of edge detector used to compute an opportune norm of the
     * gradient of \ref segmentation::phi in the functional to be minimized
     *
     * The goal is to compute a function to be multiplied with the module of the
     * gradient in order to make his L1 norm smaller in regions with edges.
     * This function is of the form \f$ g=1\big/(1+\beta\xi) \f$, where \f$\xi\f$ is
     * the edge detector, hence it assumes higher values in regions with edges.
     * The choice of \f$\xi\f$ depends on this variable that can assume three values:
     * GRAD, GRAD_THRESHOLD and LBP.\n
     * GRAD is the default (and suggested) choice that sets
     * \f$\xi=\big|\nabla\phi\big|^2\f$.
     * Most of time this is a good choice because generally edges coincide with
     * regions with high gradient, however in some kind of images
     * (for instance images with really high noise)
     * gradient could be not null also in regions far from edges. To avoid this
     * problem the second values makes a threshold on the gradient using member
     * im3d::image3d::im_to_black_and_white. User can choose the preferred
     * value of the threshold during the beginning of the execution of the algorithm.
     * Note that this choice leads to a \f$\xi\f$ equal to 1 on the 'candidate' edges,
     * 0 otherwise, hence it produces a less efficient g if at the same time \ref beta
     * is not increased properly by the user.\n
     * Both of this two choices can be combined with a Gaussian blur of the images to
     * limit error linked with noise setting \ref edge_detector_sigma with a value
     * different from zero.\n
     * LBP is the last possibility user can choose, but it works only with 2D images.
     * The effect is similar to the usage of GRAD_THRESHOLD, but to compute the black
     * and white image is used member
     * \ref im3d::image3d::local_binary_pattern_edge_detector instead of the gradient.
     * With an opportune choice of the parameters this leads to an edge detector
     * less influenced by noise. Also in this case parameters are settable
     * during the beginning of the execution of the algorithm.
     */
    ed_type edge_detector;
    /*!
     * \brief Value of sigma used to compute a Gaussian blur on the \ref edge_detector
     * in case it is not equal to LBP (default=0, no Gaussian blur).
     *
     * In case of not null value and \ref edge_detector equal to GRAD or
     * GRAD_THRESHOLD, it allows \ref edge_detect to filter the image with a Gaussian
     * kernel using class \ref conv::filtering before computing the edge detector
     * function.
     */
    T edge_detector_sigma;
    /*!
     * \brief it weighs the first addend of the functional (default=1).
     *
     * No significant changes are noted on the contour evolution modifying its value.
     */
    T nu;
    /*!
     * \brief Private functor used by \ref apply to compute one step of the unsteady
     * Poisson problem
     *
     * This functor uses by default the \ref lapl::NeuGaussSeidel function (or the
     * \ref lapl::DirGaussSeidel function in case of \ref bc equal to DIR).
     * Expert user could implements himself a different function to solve the linear
     * system and he could pass it to this functor using public member
     * \ref set_onestep_poisson. See \ref lapl::unsteady_poisson_functor class
     * documentation for further information.
     */
    lapl::unsteady_poisson_functor<T> onestep_poisson;
    /*!
     * \brief Private functor used by \ref apply to compute convolution with
     * Gaussian kernel during execution of the algorithm
     *
     * Default function associated with this functor is
     * \ref conv::acyclic_fftw_convolution
     * because is the quickest way to compute a convolution.
     * Only if \ref sigma is really small, \ref conv::direct_convolution could be
     * also a good choice. Any other choice is deprecated unless user is comfortable
     * with Fourier Transform and Convolution Theory.
     */
    conv::filtering<T> cv;
    /*!
     * \brief Boolean variable (deafult=false) to change way of computing Gaussian
     * kernel used by \ref cv in every convolution
     *
     * See \ref conv::filtering::build_gaussian_filter for further information
     */
    bool gaussian_pixel_approach;
    /*!
     * \brief tollerance value, the algorithm ends if the error is lesser than its
     * value (default 1.e-5)
     *
     * Only if the image is really easy to segment this value could cause the end of
     * the algorithm, otherwise user could stop the algorithm if partial results
     * are good using \ref onthego parameter \ref end_now or setting \ref maxiter.
     * Despite changing this value is allowed, it is not suggested.
     */
    T tol;

    ///@}
    // NOT DIRECTLY SETTABLE
    /*!
     * \name Not Directly Settable
     * \brief These attributes are used my some members of the class, but they are not
     * directly settable by user using any specialized public members
     */
    ///@{
    int ix; //!< index of initial contour set by \ref initialize_contour_as_cube
    int iy; //!< index of initial contour set by \ref initialize_contour_as_cube
    int iz; //!< index of initial contour set by \ref initialize_contour_as_cube
    int fx; //!< index of initial contour set by \ref initialize_contour_as_cube
    int fy; //!< index of initial contour set by \ref initialize_contour_as_cube
    int fz; //!< index of initial contour set by \ref initialize_contour_as_cube
    uint dimx; //!< dimension in X direction of current analyzed image
    uint dimy; //!< dimension in Y direction of current analyzed image
    uint dimz; //!< dimension in Z direction of current analyzed image
    uint space_dim; //!< equal to 3 for 3D image, 2 for 2D image
    double hx; //!< spacing in X direction of current analyzed image
    double hy; //!< spacing in Y direction of current analyzed image
    double hz; //!< spacing in Z direction of current analyzed image
    /*!
     * \brief auxiliary variable to reinitialize variables after a connected
     * component extraction
     */
    bool init_variables;

    ///@}

    // PRIVATE MEMBERS USED BY APPLY

    /*!
     \brief private member used by \ref apply to update domain of the two regions at
     each iteration
     */
    void set_heaviside (im3d::image3d<T>& Heps) const;

    /*!
     \brief private member used by \ref apply to initialize \ref segmentation::phi
     with a cube
     */
    void initialize_phi_with_cube ();

    /*!
     \brief private member used by \ref apply to initialize \ref segmentation::phi
     with an \ref im3d::image3d chosen by user
     */
    void initialize_phi_with_init (im3d::image3d<T> const& init);

    /*!
     \brief private member used by \ref apply to set the \ref edge_detector function
     used by \ref shrink
     */
    void edge_detect (im3d::image3d<T>& res, im3d::image3d<T> const& f ) const;

    /*!
     \brief private member used by \ref apply to select after solving linear system
     only those values of levelset in the allowed range [\ref a0, \ref b0]
     */
    void cut_phi();

    /*!
     \brief private member used by \ref apply to apply shrink operator
     */
    void shrink (std::vector<im3d::image3d<T> >& res,
                 std::vector<im3d::image3d<T> > const& f1, im3d::image3d<T> const& f2) const;

    /*!
     \brief private member used by \ref apply to show partial result using private
     attribute \ref segmentation::levelset and \ref segmentation::image
     */
    void internal_show ();

    /*!
     \brief private member used by \ref apply to extract a connected component of
     \ref segmentation::phi at some iteration
     */
    void extract_connected_component ();

    /*!
     \brief private member used by \ref apply to update on the go those private
     attributes settable via GetPot using \ref update_param_onthego
     */
    void update_param_onthego ();

    /*!
     \brief private member used by \ref apply to set parameters from \ref getpotfile
     */
    void set_param_from_getpot (std::string const& section);


public:

    // CONSTRUCTOR
    /*!
     \brief constructor to set all private attributes to their default values
     */
    rsfe_splitbregman();



    // MEMBERS TO SET PRIVATE PARAMETERS

    // NOT CHANGABLE WITH GETPOTFILE

    /*!
     * \name Private Parameters
     * \brief Members to set some private attributes not changable via GetPot
     */
    ///@{
    /*!
     \brief to set private attribute \ref getpotfile
     \param name is the desired filename of GetPot file
     */
    inline void set_getpotfile (std::string const& name);

    /*!
     \brief to set private attribute \ref onthego
     \param onthego is the desired value
     */
    inline void set_onthego (bool const onthego);

    /*!
     \brief to set private attribute \ref verbosity
     \param v is the desired value
     */
    inline void set_verbosity (bool const v);

    /*!
     \brief to set private attribute \ref logfilename
     \param name is the desired filename of GetPot file
     */
    inline void set_logfilename (std::string const& name);

    /*!
     \brief to set private attribute \ref cv
     \param cv is the desired function to assign at the private functor
     */
    inline void set_cv (conv::filtering<T> const cv);

    /*!
     \brief to set private attribute \ref onestep_poisson
     \param osl is the desired function to assign at the private functor
     */
    inline void set_onestep_poisson (lapl::unsteady_poisson_functor<T> const osl);

    ///@} end name Private Parameters


    // CHANGABLE WITH GETPOTFILE

    /*!
     * \name GetPot Parameters
     * \brief Members to set all private attributes also settable via \ref getpotfile
     */
    ///@{
    /*!
     \brief to set private attribute \ref maxiter
     \param maxiter is the desired value
     */
    inline void set_maxiter (uint const& maxiter);

    /*!
     \brief to set private attribute \ref showfrequency
     \param sf is the desired value
     */
    inline void set_showfrequency (uint const& sf);

    /*!
     \brief to set private attribute \ref auto_extract_conn_comp_freq
     \param cc_freq is the desired value
     */
    inline void set_auto_extract_conn_comp_freq (uint const& cc_freq);

    /*!
     \brief to set private attribute \ref dumpfrequency
     \param dump is the desired value
     */
    inline void set_dumpfrequency (uint const& dump);

    /*!
     \brief to set private attribute \ref outputname
     \param name is the desired filename of GetPot file
     */
    inline void set_outputname (std::string const& name);

    /*!
     \brief to set private attribute \ref tol
     \param tol is the desired value
     */
    inline void set_tol (T const& tol);

    /*!
     \brief to set private attribute \ref cc_init_pixel_x
     \param cc_x is the desired value
     */
    inline void set_cc_init_pixel_x (double const& cc_x);

    /*!
     \brief to set private attribute \ref cc_init_pixel_y
     \param cc_y is the desired value
     */
    inline void set_cc_init_pixel_y (double const& cc_y);

    /*!
     \brief to set private attribute \ref cc_init_pixel_z
     \param cc_z is the desired value
     */
    inline void set_cc_init_pixel_z (double const& cc_z);

    /*!
     \brief to set private attribute \ref lamda
     \param lamda is the desired value
     */
    inline void set_lamda (T const& lamda);

    /*!
     \brief to set private attribute \ref lamda1
     \param lamda1 is the desired value
     */
    inline void set_lamda1 (T const& lamda1);

    /*!
     \brief to set private attribute \ref lamda2
     \param lamda2 is the desired value
     */
    inline void set_lamda2 (T const& lamda2);

    /*!
     \brief to set private attribute \ref sigma
     \param sigma is the desired value
     */
    inline void set_sigma (T const& sigma);

    /*!
     \brief to set private attribute \ref edge_detector_sigma
     \param eds is the desired value
     */
    inline void set_edge_detector_sigma (T const& eds);

    /*!
     \brief to set private attribute \ref nu
     \param nu is the desired value
     */
    inline void set_nu (T const& nu);

    /*!
     \brief to set private attribute \ref beta
     \param beta is the desired value
     */
    inline void set_beta (T const& beta);

    /*!
     \brief to set private attribute \ref dt
     \param dt is the desired value
     */
    inline void set_dt (T const& dt);

    /*!
     \brief to set private attribute \ref a0
     \param a0 is the desired value
     */
    inline void set_a0 (T const& a0);

    /*!
     \brief to set private attribute \ref b0
     \param b0 is the desired value
     */
    inline void set_b0 (T const& b0);

    /*!
     \brief to set private attribute \ref epsilon
     \param epsilon is the desired value
     */
    inline void set_epsilon (T const& epsilon);

    /*!
     \brief to set private attribute \ref bc
     \param bc is the desired value
     */
    inline void set_bc (bc_type const bc);

    /*!
     \brief to set private attribute \ref ls_tol
     \param tol is the desired value
     */
    inline void set_ls_tol (T const& tol);

    /*!
     \brief to set private attribute \ref ls_steps
     \param steps is the desired value
     */
    inline void set_ls_steps (uint const steps);

    /*!
     \brief to set private attribute \ref gaussian_pixel_approach
     \param gpa is the desired value
     */
    inline void set_gaussian_pixel_approach (bool const gpa);

    /*!
     \brief to set private attribute \ref bc
     \param ed is the desired value
     */
    inline void set_edge_detector (ed_type const ed);

    /*!
     \brief to set private attribute \ref auto_extract_conn_comp
     \param auto_cc is the desired value
     */
    inline void set_auto_extract_conn_comp (bool const auto_cc);

    /*!
     \brief to set private attribute \ref cc_binary_output
     \param cc_bo is the desired value
     */
    inline void set_cc_binary_output (bool const cc_bo);

    /*!
     \brief to set private attribute \ref cc_init_variables
     \param cc_bo is the desired value
     */
    inline void set_cc_init_variables (bool const cc_iv);

    /*!
     \brief to set private attribute \ref segmentation::alpha
     \param alpha is the desired value
     */
    inline void set_alpha (T const& alpha);

    ///@} end name GetPot Parameters


    // MEMBERS TO APPLY ALGORITHM

    /*!
     * \name Algorithm Execution
     * \brief Members to execute algorithm
     */
    ///@{
    /*!
     \brief Member to set inizial rectangle/cube.

     This member show image you'd like to segment and allows user to select
     coordinates of the two vertices of the desired inizial rectangle/cube.
     Clicking with mouse on the image user would be able to know coordinates and to
     write it when program asks him.
     Coordinates will be saved in private variables, so when user executes algorithm
     the desired initial rectangle/cube is built.

     \param image is image you'd like to segment using \ref apply.

     \warning use this member with the same image you'd like to segment and note that
     at the end of the algorithm coordinates will be set to the default values again.
     */
    void initialize_contour_as_cube (im3d::image3d<T> const& image);

    /*!
     \brief Member to apply the algorithm.

     \param myim is the image to segment.
     */
    void apply (im3d::image3d<T> const& myim);

    /*!
     \brief Member to apply the algorithm.

     \param myim is the image to segment
     \param init is the initial value with whom levelset \ref phi is initialized
     */
    void apply (im3d::image3d<T> const& myim, im3d::image3d<T> const& init);

    ///@} end name Algorithm Execution


};

}//end namespace segm

#include "rsfe_splitbregman_imp.hxx"

#endif // RSFE_SPLITBREGMAN_HXX_INCLUDED
