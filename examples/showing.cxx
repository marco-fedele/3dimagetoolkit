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
  \file showing.cxx

  \brief File to create an executable to show a selected image in some diferent way using
  some members of class \ref im3d::interface.
*/

#include "../header/names.hxx"
#include "../header/image3d.hxx"
#include "../header/interface.hxx"


#include <string>



using namespace std;
using namespace im3d;


int main (int argc, char** argv)
{
    //Lettura filename
    string filename, backgroundimage_filename, outputname = "./results/result_";
    uint number_of_show = 1;
    double level_value;


    if (argc > 1)
    {
        filename = argv[1];
    }

    else
    {
        // showing show helper if user execute the program with a wrong number of options
        cout << endl << "\t*** show helper ***" << endl << endl;
        cout << "\tTo use show you have to execute the program with at least 1 option:";
        cout << endl << endl;
        cout << "\tshow <name of image to show> [<name of background image> <number of show> <level value>]" << endl << endl;
        cout << "\t<number of show> default value is 1. The following option are possible:" << endl;
        cout << "\t1. show a level of your image with background image" << endl;
        cout << "\t2. show your image with planes" << endl;
        cout << "\t3. show your image with planes and its level" << endl;
        cout << "\t4. show only a level" << endl;
        cout << "\t5. show your background image" << endl << endl;
        cout << "\t<level value> default value is min(image)+[max(image)-min(image)]/2." << endl << endl;
        cout << "\tIf you execute the program without specifying a background image," << endl;
        cout << "\tit will be shown only the type 3 show." << endl << endl;

        return 0;
    }

    if (argc > 2)
    {
        backgroundimage_filename = argv[2];
    }

    if (argc > 3)
    {
        number_of_show = atoi (argv[3]);
    }

    image3d<int> myim;

    interface<int> shower (filename);

    cout << "\tend reading file" << endl;


    if (argc == 2)
    {
        cout << "\tshowing your image with planes ..." << endl;
        shower.show_image();
    }
    else
    {
        interface<int> background_interface (backgroundimage_filename);

        if (argc > 4)
        {
            level_value = atof (argv[4]);
        }

        cout << endl << "\tshowing a level of your image with background image ..." << endl;
        shower.show_contour_with_background_image (background_interface, level_value);

        if (number_of_show >= 2)
        {

            shower.convert2image3d (myim);

            cout << "\tend convert2image3d" << endl;

            level_value = myim.max() - myim.min();
            level_value /= 2;
            level_value += myim.min();

            cout << endl << "\tshowing your image with planes and its level ..." << endl;
            shower.show_image_and_contour (level_value);
        }

        if (number_of_show >= 3)
        {
            cout << "\tshowing your image with planes and its level ..." << endl;
            shower.show_image_and_contour (level_value);
        }

        if (number_of_show >= 4)
        {
            cout << "\tshowing only a level of your image ..." << endl;
            shower.show_contour (level_value);
        }

        if (number_of_show >= 5)
        {
            cout << "\tshowing background image with planes ..." << endl;
            background_interface.show_image();
        }

    }

    cout << endl;

    return 0;
}
