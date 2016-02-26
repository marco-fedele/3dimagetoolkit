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
  \file image_toolkit.cxx

  \brief Program that allows user to manipulate an image using most important
  members of \ref im3d::image3d, \ref conv::filtering and \ref segm::rsfe_splitbregman.
*/

#include "../header/names.hxx"
#include "../header/image3d.hxx"
#include "../header/interface.hxx"
#include "../header/convolution.hxx"
#include "../header/rsfe_splitbregman.hxx"


#include <string>
#include <vector>
#include <iostream>
#include <cmath>

typedef double MY_REAL;		



int main (int argc, char ** argv){
	
	if (argc==1){
		// showing filtering helper if user execute the program with a wrong number of options
		std::cout << std::endl << "*** Image Toolkit Help ***" << std::endl << std::endl;
		std::cout << "To use Image Toolkit you have to execute the program with 1 option:";
		std::cout << std::endl << std::endl;
		std::cout << "image_toolkit <name of image to filter>" << std::endl << std::endl;
		
		return 0;
	}
	else{
		string filename;
		filename = argv[1];
		std::cout << std::endl << "\tselected file:\t" << filename << std::endl;
		
		// importing choosen image using class interface
		im3d::interface<MY_REAL> myim_interface;
		
		// inizializing necessary image3d
		im3d::image3d<MY_REAL> myim, output_image; 
		
		{
			im3d::interface<MY_REAL> reader(filename);
			// converting chosen image in an image3d
			reader.convert2image3d(myim);
		}
		
		std::cout << "\tdimensions:\t" << myim.getdimx() << "\tx   " << myim.getdimy() << "   \tx   ";
		std::cout <<  myim.getdimz() << std::endl;

		std::cout << "\tspacing:\t" << myim.gethx() << "\tx   " << myim.gethy() << "\tx   ";
		std::cout <<  myim.gethz() << std::endl;

		std::cout << "\trange of intensity:\t[" << myim.min() << ", " << myim.max() << "]" << std::endl;

		// change range of intensity between 0 and 255 to facility visualization on a gray scale
		//myim.change_range_of_intensity(255);
		
		// converting myim in interface to avoid problem of visualization with some kind of file
		myim_interface.convertfromimage3d(myim);
		
		// initializing a filtering functor assigning it acyclic_fftw_convolution function
		conv::filtering<MY_REAL> filter(conv::acyclic_fftw_convolution);
		
		while(1){
			std::cout << std::endl << "\tChoose something to do with your image:";
			std::cout << std::endl << std::endl;
			std::cout << "\t1) Show" << std::endl;
			std::cout << "\t2) Create Cube" << std::endl;
			std::cout << "\t3) Crop" << std::endl;
			std::cout << "\t4) Change Resolution" << std::endl;
			std::cout << "\t5) Gaussian Blur" << std::endl;
			std::cout << "\t6) Sharpness Filter" << std::endl;
			std::cout << "\t7) Laplacian" << std::endl;
			std::cout << "\t8) Modulus of Gradient" << std::endl;
			std::cout << "\t9) Median Filter" << std::endl;
			std::cout << "\t10) Local Binary Pattern " << std::endl;
			std::cout << "\t11) Threshold" << std::endl;
			std::cout << "\t12) Connected Component" << std::endl;
			std::cout << "\t13) Range of Intensity" << std::endl;
			std::cout << "\t14) Histogram Equalization" << std::endl;
			std::cout << "\t15) Segmentation Algorithm" << std::endl;
			std::cout << "\t16) Delete Binary Image" << std::endl;
			std::cout << std::endl;
			std::cout << "\t0) Exit" << std::endl << std::endl << "\t";
			
			int choice;
			
			cin >> choice;
			std::cout << std::endl;
			
			char answer;
			
			switch (choice){
			case 1:
				{
					MY_REAL level;
					std::cout << "Showing original image..." << std::endl << std::endl;
					myim_interface.show_image();
					
					/*std::cout << "Choose a level between 0 and 255 to show:" << std::endl;
					std::cout << "\tlevel: ";
					cin >> level;
					
					myim_interface.show_contour(level);
					myim_interface.show_image_and_contour(level);
					*/
					std::cout << "Do you want to show this level with a different background image? (y/n): ";
					cin >> answer;
					
					if(answer=='y'){
						string background_file;
						std::cout << "Choose an image file as background: ";
						cin >> background_file;
						
						im3d::interface<MY_REAL> background_interface(background_file);
						level = ( myim.max()-myim.min() ) / 2.;
						myim_interface.show_contour_with_background_image(background_interface,level);
					}
					
					
				}				
				break;
					
				
				
			case 2:
				{
					uint start_i,start_j,start_k,end_i,end_j,end_k,tmp;	
					std::cout << "Choose a pair of pixels as vertixes of desired cube:" << std::endl;
					myim_interface.get_coordinates(start_i,start_j,start_k,end_i,end_j,end_k);
					
					if (start_i>end_i) { tmp=end_i; end_i=start_i; start_i=tmp; }
					if (start_j>end_j) { tmp=end_j; end_j=start_j; start_j=tmp; }
					if (start_k>end_k) { tmp=end_k; end_k=start_k; start_k=tmp; }
					
					output_image.setdim(myim.getdimx(),myim.getdimy(),myim.getdimz(),0);
					output_image.seth(myim.gethx(),myim.gethy(),myim.gethz());
					
					std::cout << "\tstart: " << start_i << "-" << start_j << "-" << start_k << std::endl;
					std::cout << "\tend: " << end_i << "-" << end_j << "-" << end_k << std::endl;
					std::cout << "\tdim: " << myim.getdimx() << "-" << myim.getdimy() << "-";
					std::cout << myim.getdimz() << std::endl;
					
					for (uint i= start_i; i<= end_i; ++i)
						for (uint j= start_j; j<= end_j; ++j)
							for (uint k= start_k; k<= end_k; ++k)
								output_image(i,j,k)=1.;
					
				}				
				break;
					
					
					
			case 3:
				{
					uint start_i,start_j,start_k,end_i,end_j,end_k;	
					std::cout << "Choose a pair of pixels to limit region to Crop:" << std::endl;
					myim_interface.get_coordinates(start_i,start_j,start_k,end_i,end_j,end_k);
					
					myim.crop(output_image,start_i,start_j,start_k,end_i,end_j,end_k);
				}				
				break;
					
					
					
			case 4:
				{
					bool increase;
					int ratio;
					std::cout << "Choose parameters of Change Resolution:" << std::endl;
					
					std::cout << "\tratio: ";
					cin >> ratio;
					
					std::cout << "\tincrease: ";
					cin >> increase;
					
					myim.change_resolution(output_image,ratio,increase);
				}
				break;
					
					
					
			case 5:
				{
					double sigma, radius;
					std::cout << "Choose parameters of Gaussian filter:" << std::endl;
					
					std::cout << "\tsigma: ";
					cin >> sigma;
					
					std::cout << "\tradius (suggested 3): ";
					cin >> radius;
					std::cout << std::endl;
					
					// building gaussian filter using static member of class filtering
					conv::filtering<MY_REAL>::build_gaussian_filter(myim, sigma,radius);
					
					// executing convolution
					filter(output_image,myim);
				}				
				break;
					
					
					
			case 6:
				{
					double central_node_coefficient;
					std::cout << "Choose parameter of High-Pass filter:" << std::endl;
					
					std::cout << "\tcentral_node_coeff (suggested 1): ";
					cin >> central_node_coefficient;
					std::cout << std::endl;
					
					// building high-pass filter using static member of class filtering
					conv::filtering<MY_REAL>::build_high_pass_filter(myim, central_node_coefficient);
					
					// executing convolution
					filter(output_image,myim);
				}				
				break;
					
					
					
			case 7:
				// building laplacian filter using static member of class filtering
				conv::filtering<MY_REAL>::build_laplacian_filter(myim);
					
				// executing convolution
				filter(output_image,myim);
					
				break;
					
					
					
			case 8: 
				{
					std::vector<im3d::image3d<MY_REAL> > image_gradient;
					
					if (myim.getdimz()==1)
						image_gradient.resize(2);
					else
						image_gradient.resize(3);
                    
					myim.grad(image_gradient);
                    
					output_image.setdim(myim.getdimx(),myim.getdimy(),myim.getdimz());
					output_image.seth(myim.gethx(),myim.gethy(),myim.gethz());

					vector_abs(output_image,image_gradient);
					
				}
				break;
					
				
					
			case 9:	
				{
					uint radius=0;
					std::cout << "Choose radius of the median filter" << std::endl;
					cin >> radius;
					myim.median_filter(output_image,radius);
				}
				break;	
					
					
					
			case 10:	
				{
					MY_REAL T;
					double t1, t2;
					
					std::cout << "Choose parameters of the local binary pattern" << std::endl;
					std::cout << "\tconstant T: ";
					cin >> T;
					std::cout << "\tthreshold t1: ";
					cin >> t1;
					std::cout << "\tthreshold t2 (>t1): ";
					cin >> t2;
					myim.local_binary_pattern_edge_detector(output_image,T,t1,t2);
				}
					
				break;
					
					
					
			case 11:
				{
					double threshold;
					bool negative;
					std::cout << "Choose parameters of Binary Threshold:" << std::endl;
					std::cout << "\tthreshold coefficient (between 0 and 1): ";
					cin >> threshold;
					std::cout << "\tnegative version (0=false, 1=true): ";
					cin >> negative;
					
					myim.im_to_black_and_white(output_image,threshold,negative);
				}				
				break;
					
					
					
			case 12:
				{
					double threshold;
					bool full_connected, binary_output;
					im3d::interface<MY_REAL> binary_interface;
					uint i,j,k;
					std::cout << "Choose parameters of Connected Component:" << std::endl;
					std::cout << "\tthreshold (between 0 and 1): ";
					cin >> threshold;
					std::cout << "\tfull connected (0=false, 1=true): ";
					cin >> full_connected;
					std::cout << "\tbinary output (0=false, 1=true): ";
					cin >> binary_output;
					std::cout << "\ta pixel of desired connected component: " << std::endl;
					
					myim.im_to_black_and_white(output_image,threshold);
					
					binary_interface.convertfromimage3d(output_image);
					
					binary_interface.get_coordinates(i,j,k);
					
					myim.connected_component(output_image,i,j,k,threshold,full_connected,binary_output);
				}
				break;
					
					
					
			case 13:
				{
					double lowerbound, upperbound, user_lowervalue(0.), user_uppervalue(0.);
					int black_background;
					std::cout << "\tactual range of intensity:\t[" << myim.min() << ", " << myim.max() << "]" << std::endl;
					std::cout << "Choose parameters of Range of Intensity:" << std::endl;
					std::cout << "\tlowerbound: ";
					cin >> lowerbound;
					std::cout << "\tupperbound: ";
					cin >> upperbound;
					std::cout << "\tvalues out of range (-2=user values, -1=lowerbound, 1=upperbound, 0=hybrid): ";
					cin >> black_background;
					if (black_background == -2)
					{
						std::cout << "\tuser lowervalue: ";
						cin >> user_lowervalue;
						std::cout << "\tuser uppervalue: ";
						cin >> user_uppervalue;
					}
					myim.select_range_of_intensity(output_image, lowerbound, upperbound, black_background, user_lowervalue, user_uppervalue);
				}
				break;
					
					
					
			case 14:
				{
					uint quantization;
					std::cout << "Choose parameter of Histogram Equalization:" << std::endl;
					std::cout << "quantization (suggested at least 256): ";
					cin >> quantization;
					
					myim.histogram_equalization(output_image,quantization);
				}		  
				break;
					
					
					
			case 15:
				{
					string getpot_file, logfile;
					bool onthego;
					int choice;
					segm::rsfe_splitbregman<MY_REAL> algo;
					
					std::cout << "Choose parameters of RSFE Split Bregman Segmentation Algorithm:";
					std::cout << std::endl;
					std::cout << "\tgetpot file (digit ! to use default): ";
					cin >> getpot_file;
					
					if (getpot_file!="!")
						algo.set_getpotfile(getpot_file);
					
					std::cout << "\tlogfilename (digit ! not to set a file): ";
					cin >> logfile;
					
					if (logfile!="!")
						algo.set_logfilename(logfile);
					
					std::cout << "\ton the go interactivity (0=false, 1=true): ";
					cin >> onthego;
					
					algo.set_onthego(onthego);
					
					std::cout << "\tset initial contour:" << std::endl;
					std::cout << "\t\t1) default" << std::endl;
					std::cout << "\t\t2) customized cube" << std::endl;
					std::cout << "\t\t3) from an image" << std::endl;
					
					cin >> choice;
					
					switch (choice) {
					case 1:
						algo.apply(myim);
						break;
							
					case 2:
						algo.initialize_contour_as_cube(myim);
						algo.apply(myim);
						break;
							
					case 3:
						{
							string init_file;
							
							std::cout << "Choose a file for the initialization (with path):\t";
							
							cin >> init_file;
							
							im3d::interface<MY_REAL> init_interface(init_file);
							im3d::image3d<MY_REAL> init_image;
							init_interface.convert2image3d(init_image);
							
							algo.apply(myim,init_image);
						}
						break;
							
					default:
						break;
					}
                    
					output_image = algo.getoutput();
				}		    
				break;
			
			
			case 16:
				{
					im3d::image3d<MY_REAL> aux;

					string filename_binaryimage;

					std::cout << "file with binary image to subtract (with path): " ;
					std::cin >> filename_binaryimage;
					{
						im3d::interface<MY_REAL> reader(filename_binaryimage);
						// converting chosen image in an image3d
						reader.convert2image3d(aux);
					}
					MY_REAL min, max;
					min = myim.min();
					max = myim.max();
					output_image = myim;
					output_image.change_range_of_intensity(max-min);
					aux *= max-min;
					output_image -= aux;
					aux = output_image;
					aux.select_range_of_intensity(output_image,min,max);
				}
				
				break;


			default:
				return 0;
					
				break;
					
					
			}// end initial menu switch case
			
			
			if(choice!=1){
				
				// change range of intensity between 0 and max-min to facility visualization
				// on a gray scale
				// output_image.change_range_of_intensity(output_image.max()-output_image.min());
				
				// create output_interface
				im3d::interface<MY_REAL> output_interface(output_image);
				
				std::cout << "Showing original image..." << std::endl;
				myim_interface.show_image();
				
				std::cout << "Showing output image..." << std::endl;
				if ( output_image.getdimz()!=1 && (choice==2 || choice==12 || choice==15) )
					output_interface.show_image_and_contour((output_image.max()-output_image.min())/2);
				else
					output_interface.show_image();
				
				// allowing user to choose if he would like to save result of filtering or not
				std::cout << "Do you want to save this filtered image? (y/n): ";
				cin >> answer;
				
				if (answer=='y'){
					string output_filename;
					std::cout << "Choose name of output file (full path): ";
					cin >> output_filename;
					
					if (output_image.getdimz()==1){
						std::cout << "Choose format (j=.jpg, m=.mhd): ";
						cin >> answer;
					}
					
					if (answer=='j')
						output_interface.write(output_filename,".jpg");
					else
						output_interface.write(output_filename,".mhd");
				}
				
				// allowing user to choose if he would like to continue to modify output image
				std::cout << "Do you want to continue to modify output image? (y/n): ";
				cin >> answer;
				
				if (answer=='y'){
					myim = output_image;
					myim_interface.convertfromimage3d(myim);
				}
				
			}
			
		}//end while
		
	}// end initial else
	
	
	return 0;
}
