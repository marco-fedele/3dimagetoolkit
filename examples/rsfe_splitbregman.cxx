/*
=========================================================================================
3D Image Toolkit

Copyright 2013:
Marco Fedele, Luca Barbarotta, Francesco Cremonesi.
All rights reserved.
Read details of the license in the file license.txt provided with the library.
=========================================================================================
*/

/*!
 \file rsfe_splitbregman.cxx

 \brief File to create an executable usefull to apply \ref segm::rsfe_splitbregman algorithm
 to a \ref im3d::image3d allowing user to set all parameters of the algorithm from file
 data.
 */

//#include "../header/names.hxx"
//#include "../header/image3d.hxx"
//#include "../header/interface.hxx"
#include "../header/rsfe_splitbregman.hxx"


#include <string>
#include <iostream>


typedef double MY_REAL;



int main (int argc, char** argv)
{
    // reading
    string filename, getpotfile = "./data", logfilename = "!", initname = "!";

    if (argc > 1)
    {
        filename = argv[1];
    }
    else
    {
        // showing rsfe_splitbregman helper if user execute the program with a wrong number of options
        std::cout << std::endl << "\t*** rsfe_splitbregman helper ***" << std::endl << std::endl;
        std::cout << "\tTo use rsfe_splitbregman you have to execute the program with at least 1 option:";
        std::cout << std::endl << std::endl;
        std::cout << "\trsfe_splitbregman <name of image to segment>\n\t\t\t[<getpotfile> <logfilename>";
        std::cout << " <init phi filename>]" << std::endl << std::endl;
        std::cout << "\t<getpotfile> is the \".pot\" to import parameters of algorithm." << std::endl;
        std::cout << "\tDefault value is \"./data\"." << std::endl << std::endl;
        std::cout << "\t<logfilename> is name of logfile where history is saved. " << std::endl;
        std::cout << "\tIf it is not set terminal is used." << std::endl << std::endl;
        std::cout << "\t<init phi filename> is name of image with whom you'd like to initilize phi." << std::endl;
        std::cout << "\tIf it is not set default cube is used." << std::endl << std::endl;

        return 0;
    }


    if (argc > 2)
    {
        getpotfile = argv[2];
    }

    if (argc > 3)
    {
        logfilename = argv[3];
    }

    if (argc > 4)
    {
        initname = argv[4];
    }

    std::cout << "you have selected: " << filename << std::endl;

    im3d::image3d<MY_REAL> myim, init;

    {
        im3d::interface<MY_REAL> myim_inteface (filename);
        myim_inteface.convert2image3d (myim);
    }

    myim.change_range_of_intensity (255);

    if (initname != "!")
    {
        im3d::interface<MY_REAL> init_inteface (initname);
        init_inteface.convert2image3d (init);
    }

    segm::rsfe_splitbregman<MY_REAL> algo;

    algo.set_getpotfile (getpotfile);
    algo.set_onthego (true);

    if (logfilename != "!")
    {
        algo.set_logfilename (logfilename);
    }

    std::cout << "Applying segmentation algorithm ..." << std::endl;

    if (initname == "!")
    {
        algo.apply (myim);
    }
    else
    {
        algo.apply (myim, init);
    }


    return 0;
}
