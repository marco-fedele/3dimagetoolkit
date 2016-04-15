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
  \file  interface_imp.hxx

   \brief file with the implementation of all the methods of class \ref im3d::interface
*/

#ifndef INTERFACE_IMP_HXX_INCLUDED
#define INTERFACE_IMP_HXX_INCLUDED



template <typename S>
im3d::interface<S>::interface ()
{
    setcolour();
    setopacity();
    this->imageptr = vtkImageData::New();
}



template <typename S>
im3d::interface<S>::interface (std::string const& imagename)
{
    setcolour();
    setopacity();

    // Read a specified .mhd file.
    if ( imagename.compare (imagename.size() - 4, 4, ".mhd") == 0 )
    {
        vtkSmartPointer<vtkMetaImageReader> metaReader = vtkMetaImageReader::New();
        metaReader->SetFileName (imagename.c_str() );
        metaReader->Update();
        this->imageptr = metaReader->GetOutput();
    }

    // Read a specified .dcm file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".dcm") == 0 )
    {
        vtkSmartPointer<vtkDICOMImageReader> dicom_reader = vtkDICOMImageReader::New();
        dicom_reader->SetFileName (imagename.c_str() );
        dicom_reader->Update();
        this->imageptr = dicom_reader->GetOutput();
    }

    // Read a specified .jpg file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".jpg") == 0 ||
              imagename.compare (imagename.size() - 4, 4, ".JPG") == 0 ||
              imagename.compare (imagename.size() - 5, 5, ".jpeg") == 0 ||
              imagename.compare (imagename.size() - 5, 5, ".JPEG") == 0 )
    {
        vtkSmartPointer<vtkJPEGReader> jpeg_reader = vtkJPEGReader::New();
        jpeg_reader->SetFileName (imagename.c_str() );
        jpeg_reader->Update();
        this->imageptr = jpeg_reader->GetOutput();
    }

    // Read a specified .bmp file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".bmp") == 0 )
    {
        vtkSmartPointer<vtkBMPReader> bmp_reader = vtkBMPReader::New();
        bmp_reader->SetFileName (imagename.c_str() );
        bmp_reader->Update();
        this->imageptr = bmp_reader->GetOutput();
    }

    // Read all the DICOM files in the specified directory.
    else
    {
        vtkSmartPointer<vtkDICOMImageReader> dicom_reader = vtkDICOMImageReader::New();
        dicom_reader->SetDirectoryName (imagename.c_str() );
        dicom_reader->Update();
        this->imageptr = dicom_reader->GetOutput();
    }

    return;
}



template <typename S>
im3d::interface<S>::interface (image3d<S> const& myim)
{
    setcolour();
    setopacity();
    this->imageptr = vtkImageData::New();
    this->convertfromimage3d (myim);
}



template <typename S>
void im3d::interface<S>::convertfromimage3d (image3d<S> const& myim)
{
    this->imageptr->SetSpacing (static_cast<double> (myim.gethx() ),
                                static_cast<double> (myim.gethy() ),
                                static_cast<double> (myim.gethz() ) );
    this->imageptr->SetExtent (0, myim.getdimx() - 1, 0, myim.getdimy() - 1, 0, myim.getdimz() - 1);
#if VTK_MAJOR_VERSION <= 5
    this->imageptr->SetNumberOfScalarComponents (1);
    this->imageptr->SetScalarTypeToDouble();
#else
    this->imageptr->AllocateScalars (VTK_DOUBLE, 1);
#endif
    for (uint i = 0; i < myim.getdimx(); ++i)
        for (uint j = 0; j < myim.getdimy(); ++j)
            for (uint k = 0; k < myim.getdimz(); ++k)
                this->imageptr->SetScalarComponentFromDouble
                (i, j, k, 0, static_cast<double> (myim (i, j, k) ) );
    return;
}



template <typename S>
void im3d::interface<S>::convert2image3d (image3d<S>& myim) const
{
    int dims[3];
    double h[3];

    this->imageptr->GetDimensions (dims);

    myim.setdim (dims[0], dims[1], dims[2]);

    this->imageptr->GetSpacing (h);

    myim.sethx ( static_cast<double> (h[0]) );
    myim.sethy ( static_cast<double> (h[1]) );
    myim.sethz ( static_cast<double> (h[2]) );

    #pragma omp parallel for
    for (uint i = 0; i < myim.getdimx(); ++i)
        for (uint j = 0; j < myim.getdimy(); ++j)
            for (uint k = 0; k < myim.getdimz(); ++k)
                myim (i, j, k) =
                    static_cast<S> (this->imageptr->GetScalarComponentAsDouble (i, j, k, 0) );

    return;
}




template <typename S>
void im3d::interface<S>::show_contour (double const& level)
{
    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    setopacity();
    setcolour();
    addcontour2renderer (renderer, level);

    startshowing (renderer, renWinInt);

    deleterendering (renderer, renWinInt);

    return;
}



template <typename S>
void im3d::interface<S>::show_image ()
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];


    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);

    vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetZ = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
    imagePlaneWidgetZ->SetInputData (this->imageptr);
#else
    imagePlaneWidgetZ->SetInput (this->imageptr);
#endif
    imagePlaneWidgetZ->SetInteractor (renWinInt);
    imagePlaneWidgetZ->SetPlaneOrientationToZAxes();
    imagePlaneWidgetZ->DisplayTextOn();
    imagePlaneWidgetZ->UseContinuousCursorOn();
    imagePlaneWidgetZ->On();

    if (dimz != 1)
    {
        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetX = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetX->SetInputData (this->imageptr);
#else
        imagePlaneWidgetX->SetInput (this->imageptr);
#endif
        imagePlaneWidgetX->SetInteractor (renWinInt);
        imagePlaneWidgetX->SetPlaneOrientationToXAxes();
        imagePlaneWidgetX->DisplayTextOn();
        imagePlaneWidgetX->UseContinuousCursorOn();
        imagePlaneWidgetX->On();

        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetY = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetY->SetInputData (this->imageptr);
#else
        imagePlaneWidgetY->SetInput (this->imageptr);
#endif
        imagePlaneWidgetY->SetInteractor (renWinInt);
        imagePlaneWidgetY->SetPlaneOrientationToYAxes();
        imagePlaneWidgetY->DisplayTextOn();
        imagePlaneWidgetY->UseContinuousCursorOn();
        imagePlaneWidgetY->On();
    }

    interface::startshowing (renderer, renWinInt);

    setopacity();
    setcolour();

    deleterendering (renderer, renWinInt);
}



template <typename S>
void im3d::interface<S>::get_coordinates (double& x, double& y, double& z)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    std::cout << "Showing your image ..." << std::endl;
    std::cout << "Click with mouse on the image to know coordinates." << std::endl;
    std::cout << "When you have chosen a pixel write ";
    std::cout << "coordinates (positive) in a piece of paper and close the image.";
    std::cout << std::endl << std::endl;

    this->show_image();

    // 2d case
    if (dimz == 1)
    {
        std::cout << "Write coordinates of the pixel (x,y):" << std::endl;

        do
        {
            std::cout << "\tx: ";
            cin >> x;
        }
        while (x < 0);

        do
        {
            std::cout << "\ty: ";
            cin >> y;
        }
        while (y < 0);

        z = 0;
    }
    // 3d case
    else
    {
        std::cout << "Write coordinates of the pixel (x,y,z):" << std::endl;

        do
        {
            std::cout << "\tx: ";
            cin >> x;
        }
        while (x < 0);

        do
        {
            std::cout << "\ty: ";
            cin >> y;
        }
        while (y < 0);

        do
        {
            std::cout << "\tz: ";
            cin >> z;
        }
        while (z < 0);
    }


    return;
}



template <typename S>
void im3d::interface<S>::get_coordinates (uint& i, uint& j, uint& k)
{
    double x, y, z;
    int dims[3];
    double h[3];

    this->imageptr->GetDimensions (dims);
    this->imageptr->GetSpacing (h);


    this->get_coordinates (x, y, z);

    i = floor (x / h[0]);
    j = floor (y / h[1]);

    if (dims[2] == 1)
    {
        k = 0;
    }
    else
    {
        k = floor (z / h[2]);
    }

    if ( static_cast<int> (i) >= dims[0] )
    {
        i = dims[0] - 1;
    }
    if ( static_cast<int> (j) >= dims[1] )
    {
        j = dims[1] - 1;
    }
    if ( static_cast<int> (k) >= dims[2] )
    {
        k = dims[2] - 1;
    }

    return;
}



template <typename S>
void im3d::interface<S>::get_coordinates (double& x1, double& y1, double& z1,
                                          double& x2, double& y2, double& z2)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    std::cout << "Showing your image ..." << std::endl;
    std::cout << "Click with mouse on the image to know coordinates." << std::endl;
    std::cout << "When you have chosen a pair of pixels write ";
    std::cout << "coordinates (positive) in a piece of paper and close the image.";
    std::cout << std::endl << std::endl;

    this->show_image();

    // 2d case
    if (dimz == 1)
    {
        std::cout << "Write coordinates of the first pixel (x1,y1):" << std::endl;

        do
        {
            std::cout << "\tx1: ";
            cin >> x1;
        }
        while (x1 < 0);

        do
        {
            std::cout << "\ty1: ";
            cin >> y1;
        }
        while (y1 < 0);

        std::cout << "Write coordinates of the second pixel (x2,y2):" << std::endl;

        do
        {
            std::cout << "\tx2: ";
            cin >> x2;
        }
        while (x2 < 0);

        do
        {
            std::cout << "\ty2: ";
            cin >> y2;
        }
        while (y2 < 0);

        z1 = z2 = 0;
    }
    // 3d case
    else
    {
        std::cout << "Write coordinates of the first pixel (x1,y1,z1):" << std::endl;

        do
        {
            std::cout << "\tx1: ";
            cin >> x1;
        }
        while (x1 < 0);

        do
        {
            std::cout << "\ty1: ";
            cin >> y1;
        }
        while (y1 < 0);

        do
        {
            std::cout << "\tz1: ";
            cin >> z1;
        }
        while (z1 < 0);

        std::cout << "Write coordinates of the second pixel (x2,y2,z2):" << std::endl;

        do
        {
            std::cout << "\tx2: ";
            cin >> x2;
        }
        while (x2 < 0);

        do
        {
            std::cout << "\ty2: ";
            cin >> y2;
        }
        while (y2 < 0);

        do
        {
            std::cout << "\tz2: ";
            cin >> z2;
        }
        while (z2 < 0);
    }

    return;
}



template <typename S>
void im3d::interface<S>::get_coordinates (uint& i1, uint& j1, uint& k1,
                                          uint& i2, uint& j2, uint& k2)
{
    double x1, y1, z1, x2, y2, z2;
    int dims[3];
    double h[3];

    this->imageptr->GetDimensions (dims);
    this->imageptr->GetSpacing (h);

    this->get_coordinates (x1, y1, z1, x2, y2, z2);

    i1 = floor (x1 / h[0]);
    i2 = floor (x2 / h[0]);
    j1 = floor (y1 / h[1]);
    j2 = floor (y2 / h[1]);

    if (dims[2] == 1)
    {
        k1 = k2 = 0;
    }
    else
    {
        k1 = floor (z1 / h[2]);
        k2 = floor (z2 / h[2]);
    }

    if ( static_cast<int> (i1) >= dims[0] )
    {
        i1 = dims[0] - 1;
    }
    if ( static_cast<int> (j1) >= dims[1] )
    {
        j1 = dims[1] - 1;
    }
    if ( static_cast<int> (k1) >= dims[2] )
    {
        k1 = dims[2] - 1;
    }

    if ( static_cast<int> (i2) >= dims[0] )
    {
        i2 = dims[0] - 1;
    }
    if ( static_cast<int> (j2) >= dims[1] )
    {
        j2 = dims[1] - 1;
    }
    if ( static_cast<int> (k2) >= dims[2] )
    {
        k2 = dims[2] - 1;
    }

    return;
}



template <typename S>
void im3d::interface<S>::show_image_and_contour (double const& level)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    setopacity();
    setcolour();
    addcontour2renderer (renderer, level);

    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);

    vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetZ = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
    imagePlaneWidgetZ->SetInputData (this->imageptr);
#else
    imagePlaneWidgetZ->SetInput (this->imageptr);
#endif
    imagePlaneWidgetZ->SetInteractor (renWinInt);
    imagePlaneWidgetZ->SetPlaneOrientationToZAxes();
    imagePlaneWidgetZ->DisplayTextOn();
    imagePlaneWidgetZ->UseContinuousCursorOn();
    imagePlaneWidgetZ->On();

    if (dimz != 1)
    {
        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetX = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetX->SetInputData (this->imageptr);
#else
        imagePlaneWidgetX->SetInput (this->imageptr);
#endif
        imagePlaneWidgetX->SetInteractor (renWinInt);
        imagePlaneWidgetX->SetPlaneOrientationToXAxes();
        imagePlaneWidgetX->DisplayTextOn();
        imagePlaneWidgetX->UseContinuousCursorOn();
        imagePlaneWidgetX->On();

        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetY = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetY->SetInputData (this->imageptr);
#else
        imagePlaneWidgetY->SetInput (this->imageptr);
#endif
        imagePlaneWidgetY->SetInteractor (renWinInt);
        imagePlaneWidgetY->SetPlaneOrientationToYAxes();
        imagePlaneWidgetY->DisplayTextOn();
        imagePlaneWidgetY->UseContinuousCursorOn();
        imagePlaneWidgetY->On();
    }

    interface::startshowing (renderer, renWinInt);

    setopacity();
    setcolour();

    deleterendering (renderer, renWinInt);

    return;
}



template <typename S>
void im3d::interface<S>::show_contour_with_background_image (interface& backgroundimage,
                                                             double const& level)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    setopacity();
    setcolour();
    addcontour2renderer (renderer, level);

    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);

    vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetZ = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
    imagePlaneWidgetZ->SetInputData (backgroundimage.getimageptr() );
#else
    imagePlaneWidgetZ->SetInput (backgroundimage.getimageptr() );
#endif
    imagePlaneWidgetZ->SetInteractor (renWinInt);
    imagePlaneWidgetZ->SetPlaneOrientationToZAxes();
    imagePlaneWidgetZ->DisplayTextOn();
    imagePlaneWidgetZ->UseContinuousCursorOn();
    imagePlaneWidgetZ->On();

    if (dimz != 1)
    {
        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetX = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetX->SetInputData (backgroundimage.getimageptr() );
#else
        imagePlaneWidgetX->SetInput (backgroundimage.getimageptr() );
#endif
        imagePlaneWidgetX->SetInteractor (renWinInt);
        imagePlaneWidgetX->SetPlaneOrientationToXAxes();
        imagePlaneWidgetX->DisplayTextOn();
        imagePlaneWidgetX->UseContinuousCursorOn();
        imagePlaneWidgetX->On();

        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetY = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetY->SetInputData (backgroundimage.getimageptr() );
#else
        imagePlaneWidgetY->SetInput (backgroundimage.getimageptr() );
#endif
        imagePlaneWidgetY->SetInteractor (renWinInt);
        imagePlaneWidgetY->SetPlaneOrientationToYAxes();
        imagePlaneWidgetY->DisplayTextOn();
        imagePlaneWidgetY->UseContinuousCursorOn();
        imagePlaneWidgetY->On();
    }

    interface::startshowing (renderer, renWinInt);

    setopacity();
    setcolour();

    deleterendering (renderer, renWinInt);

    return;
}



template <typename S>
void im3d::interface<S>::write (std::string const& imagename, std::string const& extension) const
{
    if (extension == ".mhd")
    {
        vtkSmartPointer<vtkMetaImageWriter> metaWriter = vtkMetaImageWriter::New();
        metaWriter->SetFileName ( (imagename + extension).c_str() );
        metaWriter->SetRAWFileName ( (imagename + ".raw").c_str() );
#if VTK_VERSION_MAJOR >= 6
        metaWriter->SetInputData (this->imageptr);
#else
        metaWriter->SetInput (this->imageptr);
#endif
        metaWriter->SetCompression (false);
        metaWriter->Write();
    }
    else if (extension == ".jpg")
    {
        vtkSmartPointer<vtkImageCast> castFilter = vtkSmartPointer<vtkImageCast>::New();
        castFilter->SetOutputScalarTypeToUnsignedChar ();
#if VTK_VERSION_MAJOR >= 6
        castFilter->SetInputData (this->imageptr);
#else
        castFilter->SetInput (this->imageptr);
#endif
        castFilter->Update();

        vtkSmartPointer<vtkJPEGWriter> writer = vtkSmartPointer<vtkJPEGWriter>::New();
        writer->SetFileName ( (imagename + extension).c_str() );
        writer->SetInputConnection (castFilter->GetOutputPort() );
        writer->Write();
    }
    else
    {
        std::cout << "ERROR: " << imagename << " has an extension " << extension;
        std::cout << " not allowed in im3d::interface<S>::write" << std::endl;
    }

    return;
}



template <typename S>
void im3d::interface<S>::addcontour2renderer (vtkSmartPointer<vtkRenderer>& renderer,
                                              double const& level)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    vtkSmartPointer<vtkPolyData> surface;

    // 3d case
    if (dimz != 1)
    {

        vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkMarchingCubes::New();

#if VTK_VERSION_MAJOR >= 6
        marchingCubes->SetInputData (this->imageptr);
#else
        marchingCubes->SetInput (this->imageptr);
#endif
        marchingCubes->SetValue (0, level);
        marchingCubes->Update();

        surface = marchingCubes->GetOutput();
    }

    // 2d case
    else
    {
        setcolour (1., 0., 0.);
        setopacity (1.);

        vtkSmartPointer<vtkContourFilter> contour = vtkContourFilter::New();

#if VTK_VERSION_MAJOR >= 6
        contour->SetInputData (this->imageptr);
#else
        contour->SetInput (this->imageptr);
#endif
        contour->SetValue (0, level);
        contour->Update();

        surface = contour->GetOutput();
    }

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkPolyDataMapper::New();
#if VTK_VERSION_MAJOR >= 6
    mapper->SetInputData (surface);
#else
    mapper->SetInput (surface);
#endif
    mapper->ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> actor = vtkActor::New();

    actor->SetMapper (mapper);
    actor->GetProperty()->SetOpacity (opacity);
    actor->GetProperty()->SetColor (colour[0], colour[1], colour[2]);

    renderer->AddActor (actor);

    return;
}



template <typename S>
void im3d::interface<S>::addcontour2renderer (double const& level,
                                              vtkSmartPointer<vtkRenderer>& renderer)
{
    addcontour2renderer (renderer, level);

    return;
}



template <typename S>
void im3d::interface<S>::addplanes2renderwindowinteractor
(vtkSmartPointer<vtkRenderer>& renderer,
 vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt) const
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);

    vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetZ = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
    imagePlaneWidgetZ->SetInputData (this->imageptr);
#else
    imagePlaneWidgetZ->SetInput (this->imageptr);
#endif
    imagePlaneWidgetZ->SetInteractor (renWinInt);
    imagePlaneWidgetZ->SetPlaneOrientationToZAxes();
    imagePlaneWidgetZ->DisplayTextOn();
    imagePlaneWidgetZ->UseContinuousCursorOn();
    imagePlaneWidgetZ->On();

    if (dimz != 1)
    {
        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetX = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetX->SetInputData (this->imageptr);
#else
        imagePlaneWidgetX->SetInput (this->imageptr);
#endif
        imagePlaneWidgetX->SetInteractor (renWinInt);
        imagePlaneWidgetX->SetPlaneOrientationToXAxes();
        imagePlaneWidgetX->DisplayTextOn();
        imagePlaneWidgetX->UseContinuousCursorOn();
        imagePlaneWidgetX->On();

        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetY = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetY->SetInputData (this->imageptr);
#else
        imagePlaneWidgetY->SetInput (this->imageptr);
#endif
        imagePlaneWidgetY->SetInteractor (renWinInt);
        imagePlaneWidgetY->SetPlaneOrientationToYAxes();
        imagePlaneWidgetY->DisplayTextOn();
        imagePlaneWidgetY->UseContinuousCursorOn();
        imagePlaneWidgetY->On();
    }

    return;
}



template <typename S>
void im3d::interface<S>::newrendering
(vtkSmartPointer<vtkRenderer>& renderer,
 vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt)
{
    renWinInt = vtkRenderWindowInteractor::New();
    renderer = vtkRenderer::New();

    //renWinInt.SetInteractorStyle(vtkInteractorStyleTrackballCamera::New());

    return;
}



template <typename S>
void im3d::interface<S>::deleterendering
(vtkSmartPointer<vtkRenderer>& renderer,
 vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt)
{
    renderer->Delete();
    renWinInt->Delete();

    return;
}



template <typename S>
void im3d::interface<S>::reset (vtkSmartPointer<vtkRenderer>& renderer,
                                vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt)
{
    deleterendering (renderer, renWinInt);
    newrendering (renderer, renWinInt);

    return;
}



template <typename S>
void im3d::interface<S>::startshowing (vtkSmartPointer<vtkRenderer>& renderer,
                                       vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt,
                                       int const& windowwidth, int const& windowheight)
{
    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);
    renWin->SetSize (windowwidth, windowheight);
    renderer->SetBackground (0.25, 0.25, 0.75);

    vtkSmartPointer<vtkInteractorStyleTrackballCamera> trackStyle = vtkInteractorStyleTrackballCamera::New();
    renWinInt->SetInteractorStyle (trackStyle);

    renWinInt->Initialize();
    renWin->Render();
    renWinInt->Start();

    return;
}


template <typename S>
void im3d::interface<S>::startshowing (int const& windowwidth, int const& windowheight,
                                       vtkSmartPointer<vtkRenderer>& renderer,
                                       vtkSmartPointer<vtkRenderWindowInteractor>& renWinInt)
{
    startshowing (renderer, renWinInt, windowwidth, windowheight);

    return;
}



#endif // INTERFACE_IMP_HXX_INCLUDED
