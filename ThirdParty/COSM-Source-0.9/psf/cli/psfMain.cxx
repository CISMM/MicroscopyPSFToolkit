/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/

#include "psf/psf.h"
#include "wu/wuHeader.h"
#include "xml/psfXml.h"
#include <iostream>
#include <blitz/timer.h>
#include <tclap/CmdLine.h>

#define XSTR(s) STR(s)
#define STR(s) #s
#define VERSION XSTR(COSM_VERSION)

using namespace blitz;
using namespace TCLAP;
using namespace cosm;

template<typename T>
void evaluate( const std::string& file ) {
    PsfUser user;
    blitz::Timer timer;
    PsfXml<T> psfXml;
    if ( !psfXml.open(file) )
    {
        std::cout <<"PSF generation failed - file cannot be opened"<< std::endl;;
        return;
    }
    Psf<T> psf(psfXml.model(), psfXml.type(), psfXml.eval(), &user);  
    psf.parameters( psfXml);
    timer.start();
    psf.evaluate();
    timer.stop();
    Array<T,3> psfArray = psf.psf();
    wuDataWrite(psfArray, psfXml.file());
    std::cout <<"PSF generation took "<< timer.elapsedSeconds() << " seconds "<< std::endl;;
}

int main(
    int argc,
    char* argv[]
) {
    // Define command line object
    CmdLine cmdLine("COSM PSF Generation",' ', VERSION);

    // Define the argument options
    ValueArg<std::string> xmlArg("f","file", "PSF XML file", true, "psf.xml", "psf xml");
    SwitchArg doubleArg("d", "double", "Use double precision", false);
    SwitchArg longdoubleArg("l", "long", "Use long double precision", false);

    // add arguments to command line options
    cmdLine.add(xmlArg);
    cmdLine.add(doubleArg);
    cmdLine.add(longdoubleArg);
                                                                                
    // parse command line
    cmdLine.parse(argc, argv);

    if ( doubleArg.getValue() && longdoubleArg.getValue() ) 
    {
        evaluate<long double>(xmlArg.getValue());
    }
    else if ( doubleArg.getValue() ) 
    {
        evaluate<double>(xmlArg.getValue());
    } 
        else 
    {
        evaluate<float>(xmlArg.getValue());
    }

    return 0;
}
