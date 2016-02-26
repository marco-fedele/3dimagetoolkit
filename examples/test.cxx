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
  \file test.cxx

  \brief a file to debug the library
*/

//#include "../header/names.hxx"
#include "../header/image3d.hxx"
//#include "../header/interface.hxx"
//#include "../header/convolution.hxx"
//#include "../header/poisson.hxx"
//#include "../header/rsfe_splitbregman.hxx"


//#include <string>
// #include <vector>
// #include <iostream>
// #include <cmath>

// typedef double MY_REAL;		

// using namespace std;
// using namespace im3d;
// using namespace segm;
// using namespace conv;
// using namespace lapl;

void show(im3d::image3d<double> const &);
void enumerate(im3d::image3d<double> &);
void chess(im3d::image3d<double> &);

int main (int argc, char ** argv)
{
	// write here the various classes of the library

	// COSTRUTTORI
	im3d::image3d<double> image1;
	im3d::image3d<double> image2(3,3,3);
	im3d::image3d<double> image3(image2);

	std::cout << "SET IMAGE1 TO 4X4X4" << std::endl;
	image1.setdimx(4); image1.setdimy(4); image1.setdimz(4);

	std::cout << "IMAGE1 DIMX X DIMY X DIMZ:\t";
	std::cout << image1.getdimx() << " X "; 
	std::cout << image1.getdimy() << " X ";
	std::cout << image1.getdimz() << " X ";
	std::cout << std::endl;

	std::cout << "GET IMAGE1 SPACING" << std::endl;
	std::cout << "IMAGE1 HX X HY X HZ:\t";
	std::cout << image1.gethx() << " X "; 
	std::cout << image1.gethy() << " X ";
	std::cout << image1.gethz() << " X ";
	std::cout << std::endl;

	std::cout << "SET IMAGE1 TO 3X3X3 WITH VALUES 1.5" << std::endl;
	image1.setdim(3,3,3,1.5);
	
	std::cout << "IMAGE1 DIMX X DIMY X DIMZ:\t";
	std::cout << image1.getdimx() << " X "; 
	std::cout << image1.getdimy() << " X ";
	std::cout << image1.getdimz() << " X ";
	std::cout << std::endl;

	std::cout << "SET SPACING OF IMAGE1 EQUAL TO 2" << std::endl;
	image1.sethx(2.);
	image1.sethy(2.);
	image1.sethz(2.);

	std::cout << "GET IMAGE1 SPACING" << std::endl;
	std::cout << "IMAGE1 HX X HY X HZ:\t";
	std::cout << image1.gethx() << " X "; 
	std::cout << image1.gethy() << " X ";
	std::cout << image1.gethz() << " X ";
	std::cout << std::endl;

	std::cout << "SET UNITARY SPACING TO IMAGE1" << std::endl;
	image1.seth(1.,1.,1.);
	std::cout << "GET IMAGE1 SPACING" << std::endl;
	std::cout << "IMAGE1 HX X HY X HZ:\t";
	std::cout << image1.gethx() << " X "; 
	std::cout << image1.gethy() << " X ";
	std::cout << image1.gethz() << " X ";
	std::cout << std::endl;

	// OPERATOR () CONST
	std::cout << "image1(1,1,1) --> 1.5 " << std::endl;
	std::cout << image1(1,1,1) << std::endl;

	// OPERATOR ()
	double a(image1(1,1,1));
	std::cout << "a=image1(1,1,1) --> 1.5 " << std::endl;
	std::cout << a << std::endl;
	
	// OPERATOR =
	image2=image1;
	image3=3.;

	std::cout << "AFTER OPERATOR=" << std::endl;
	std::cout << "SHOWING IMAGE1 = 1.5" << std::endl;
	show(image1);
	std::cout << "SHOWING IMAGE2 = 1.5" << std::endl;
	show(image2);
	std::cout << "SHOWING IMAGE3 = 3." << std::endl;
	show(image3);

	// OPERATOR +=
	image1+=image2;
	image3+=1.5;

	std::cout << "AFTER OPERATOR+=" << std::endl;
	std::cout << "SHOWING IMAGE1 = 3." << std::endl;
	show(image1);
	std::cout << "SHOWING IMAGE2 = 1.5" << std::endl;
	show(image2);
	std::cout << "SHOWING IMAGE3 = 4.5" << std::endl;
	show(image3);

	// OPERATOR -=
	image1-=image2;
	image3-=1.5;

	std::cout << "AFTER OPERATOR-=" << std::endl;
	std::cout << "SHOWING IMAGE1 = 1.5." << std::endl;
	show(image1);
	std::cout << "SHOWING IMAGE2 = 1.5" << std::endl;
	show(image2);
	std::cout << "SHOWING IMAGE3 = 3." << std::endl;
	show(image3);

	// OPERATOR *=
	image1*=image2;
	image3*=1.5;

	std::cout << "AFTER OPERATOR*=" << std::endl;
	std::cout << "SHOWING IMAGE1 = 2.25" << std::endl;
	show(image1);
	std::cout << "SHOWING IMAGE2 = 1.5" << std::endl;
	show(image2);
	std::cout << "SHOWING IMAGE3 = 4.5" << std::endl;
	show(image3);

	// OPERATOR /=
	image1/=image2;
	image3/=1.5;

	std::cout << "AFTER OPERATOR/=" << std::endl;
	std::cout << "SHOWING IMAGE1 = 1.5" << std::endl;
	show(image1);
	std::cout << "SHOWING IMAGE2 = 1.5" << std::endl;
	show(image2);
	std::cout << "SHOWING IMAGE3 = 3." << std::endl;
	show(image3);

	// FRIEND FUNCTION +,-,*,/
	std::cout << "SHOWING RESULT OF + = 3." << std::endl;
	image3= image1 + image2;
	show(image3);
	image3 = image1 + 1.5;
	show(image3);

	std::cout << "SHOWING RESULT OF - = 0." << std::endl;
	image3= image1 - image2;
	show(image3);
	image3 = image1 - 1.5;
	show(image3);

	std::cout << "SHOWING RESULT OF * = 2.25" << std::endl;
	image3= image1 * image2;
	show(image3);
	image3 = image1 * 1.5;
	show(image3);

	std::cout << "SHOWING RESULT OF / = 1." << std::endl;
	image3= image1 / image2;
	show(image3);
	image3 = image1 / 1.5;
	show(image3);


	image3=1.;
	// SCALAR PROD = SCALAR PROD L2
	std::cout << "SCALARPROD IMAGE1*1=40.5" << std::endl;
	std::cout << im3d::scalarprod(image1,image3) << std::endl;
	std::cout << "SCALARPROD_L2 IMAGE1*1=40.5" << std::endl;
	std::cout << im3d::scalarprod_L2(image1,image3) << std::endl;

	enumerate(image3);
	std::cout << "IMAGE3 = ENUMERATE " << std::endl;
	show(image3);
	uint start(1), end(2);
	double dstart(-0.5), dend(0.5);
	// CROP
	std::cout << " UINT CROP IMAGE3 --> LAST NUMBER" << std::endl;
	image3.crop(image1,start,start,start,end,end,end);
	show(image1);
	std::cout << " DOUBLE CROP IMAGE3 --> FIRST NUMBER" << std::endl;
	image3.crop(image1,dstart,dstart,dstart,dend,dend,dend);
	show(image1);

	// CHANGE RESOLUTION
	std::cout << "CHANGE RESOLUTION --> CORNER OF FIRST AND LAST MATRIX" << std::endl;
	image3.change_resolution(image1,2);
	show(image1);

	// CHANGE RESOLUTION
	std::cout << "CHANGE RESOLUTION --> CORNER OF FIRST AND LAST MATRIX" << std::endl;
	image3.change_resolution(image1,2,true);
	show(image1);
	
	// GRAD
	std::vector<im3d::image3d<double> > gradient;//(3);
	// gradient[0]=image1;
	// gradient[1]=image1;
	// gradient[2]=image1;
	std::cout << "GRAD OF IMAGE3" << std::endl;
	image3.grad(gradient);
	std::cout << "GRAD-COMPONENT 1 --> 1" << std::endl;
	show(gradient[0]);
	std::cout << "GRAD-COMPONENT 2 --> 3" << std::endl;
	show(gradient[1]);
	std::cout << "GRAD-COMPONENT 3 --> 9" << std::endl;
	show(gradient[2]);

	// DIV
	std::cout << "DIV OF GRADIENT --> 0" << std::endl;
	im3d::div(image1,gradient);
	show(image1);

	// NORME
	std::cout <<"NORM1 --> 351\t" << image3.norm1() << std::endl;
	std::cout <<"NORMl1 --> 351\t" << image3.normL1() << std::endl;
	std::cout <<"???\t" << image3.norm2() << std::endl;
	std::cout <<"???\t" << image3.normL2() << std::endl;
	std::cout <<"26\t" << image3.norminf() << std::endl;
	std::cout <<"26\t" << image3.max() << std::endl;
	std::cout <<"0\t" << image3.min() << std::endl;

	// HISTOGRAM EQ
	std::cout << "HISTOGRAM EQUALIZATION --> ENUMERATE" << std::endl;
	image3.histogram_equalization(image1,27);
	show(image1);

	// IM TO BW
	std::cout << "IM_TO_BLACK_AND_WHITE --> " << std::endl;
	image3.im_to_black_and_white(image1,0.5);
	show(image1);
	image3.im_to_black_and_white(image1,0.5,true);
	show(image1);

	show(image3);
	// CHANGE RANGE
	std::cout << "CHANGE RANGE --> " << std::endl;
	image1=image3;
	image1.change_range_of_intensity(20,10);
	show(image1);

	// SELECT RANGE
	std::cout << "SELECT RANGE --> " << std::endl;
	image3.select_range_of_intensity(image1,10,20,-1);
	show(image1);
	image3.select_range_of_intensity(image1,10,20,0);
	show(image1);
	image3.select_range_of_intensity(image1,10,20,1);
	show(image1);

	chess(image2);
	im3d::image3d<double> image4(4,4,4),image5(4,4,4);
	enumerate(image4);
	// CONNECTED COMPONENT
	std::cout << "CONNECTED COMPONENT --> " << std::endl;
	image4.connected_component(image5,3,3,3,0.5,false);
	show(image5);
	image4.connected_component(image5,3,3,3,0.5,true);
	show(image5);
	std::cout << "CONNECTED COMPONENT --> " << std::endl;
	image2.connected_component(image1,3,3,3,0.5,false);
	show(image1);
	image2.connected_component(image1,3,3,3,0.5,true);
	show(image1);

	// MEDIAN FILTER
	image4=1.;
	std::cout << "MEDIAN FILTER --> " << std::endl;
	image4.median_filter(image5,1);
	show(image5);

	
  return 0;	
}


void show(im3d::image3d<double> const & image)
{
	std::cout << std::endl;
	
	for (uint j=0; j< image.getdimy(); ++j)
		{
			for (uint k=0; k< image.getdimz(); ++k)
				{
					
					for (uint i=0; i< image.getdimx(); ++i)
						{
							std::cout << "\t" << image(i,j,k);
						}
					std::cout << "\t\t";
				}
			std::cout << std::endl;
		}

	std::cout << std::endl;

	return;
}

void enumerate (im3d::image3d<double> & result)
{
	for (uint k=0; k< result.getdimz(); ++k)
		for (uint j=0; j< result.getdimy(); ++j)
			for (uint i=0; i< result.getdimx(); ++i)
				result(i,j,k) = i + j*result.getdimx() + k*result.getdimx()*result.getdimy();

	return;
}

void chess (im3d::image3d<double> & result)
{
	for (uint k=0; k< result.getdimz(); ++k)
		for (uint j=0; j< result.getdimy(); ++j)
			for (uint i=0; i< result.getdimx(); ++i)
				if( (i + j*result.getdimx() +
				     k*result.getdimx()*result.getdimy() ) %2 == 0 )
					result(i,j,k)=1.;
				else
					result(i,j,k)=0.;

	return;
}
