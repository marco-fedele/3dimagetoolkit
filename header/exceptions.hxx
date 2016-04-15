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
 \file  exceptions.hxx

 \brief file with some methods useful in debug mode to capture floating point exceptions
 */

#ifndef EXCEPTIONS_HXX_INCLUDED
#define EXCEPTIONS_HXX_INCLUDED

#include <fenv.h>
#include <exception>

// flag to toggle automatic abort at Floating Point Exceptions
#ifdef FPE_ABORT
struct FpeTrap
{
    //tells the compiler to automatically call this function at the beginning of execution
    static void __attribute__ ( (constructor) )
    trapfpe()
    {
        feenableexcept (FE_INVALID | FE_DIVBYZERO /*| FE_OVERFLOW | FE_UNDERFLOW*/);
    }
};
#endif

/* Definition of Exceptions */
struct BadFPOper : public std::exception
{
    virtual char const* what() const throw()
    {
        return "Error::Invalid FP Operation: Possibly 0/0 or sqrt of negative number. \
Program is not stopping. Compile with the option -DFPE_ABORT to automatically \
abort program at Floating Point Operation Exceptions ";
    }
};

struct BadDivision : public std::exception
{
    virtual char const* what() const throw()
    {
        return "Error::Division by Zero. Program is not stopping. Compile with the option \
-DFPE_ABORT to automatically abort program at Floating Point Operation Exceptions";
    }
};


void test_fpe_exception() throw (BadFPOper, BadDivision)
{
    int set_excepts = fetestexcept (FE_INVALID | FE_DIVBYZERO);
    feclearexcept (FE_INVALID | FE_DIVBYZERO);
    if (set_excepts & FE_INVALID)
    {
        throw BadFPOper();
    }
    if (set_excepts & FE_DIVBYZERO)
    {
        throw BadDivision();
    }
}

#endif

