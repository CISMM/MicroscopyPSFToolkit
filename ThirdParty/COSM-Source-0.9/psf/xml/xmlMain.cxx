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

#include "psfXml.h"
#include <iostream>
#include <tclap/CmdLine.h>

using namespace cosm;
using namespace TCLAP;

//#define REAL double
#define REAL float

int main(
    int argc,
    char* argv[]
) {
    // Define command line object
    CmdLine cmdLine("COSM PSF Generation",' ', "0.1");

    // Define the argument options
    ValueArg<std::string> inputXmlArg("i","input", "Input XML file", true, "in.xml", "input xml");
    ValueArg<std::string> outputXmlArg("o","output", "Output XML file", true, "out.xml", "output xml");

    // add arguments to command line options
    cmdLine.add(inputXmlArg);
    cmdLine.add(outputXmlArg);
                                                                                
    // parse command line
    cmdLine.parse(argc, argv);

    PsfXml<float> psfXml;
    psfXml.open(inputXmlArg.getValue());
    psfXml.save(outputXmlArg.getValue());

    return 0;
}
