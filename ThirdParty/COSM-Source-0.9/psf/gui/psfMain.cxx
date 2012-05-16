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

#include "psfMain.h"
#include "wu/wuHeader.h"
#include <iostream>
#include <blitz/timer.h>
#include "blitz/arrayManip.h"
#include <tclap/CmdLine.h>
#include "psfGUI.h"
#include <FL/Fl.H>

#define XSTR(s) STR(s)
#define STR(s) #s
#define VERSION XSTR(COSM_VERSION)

using namespace blitz;
using namespace TCLAP;
using namespace cosm;

template<typename T>
void evaluate( const PsfXml<T>& psfXml, PsfUser* user ) {
    Psf<T> psf(psfXml.model(), psfXml.type(), psfXml.eval(), user);  
    psf.parameters( psfXml);
    blitz::Timer timer;
    timer.start();
    std::cout <<"Enter PSFMain::evaluate()"<<std::endl;
    psf.evaluate(MAGNITUDE|REAL|IMAGINARY);


    std::cout <<"PSFMain psf.evaluate()  completed"<< std::endl;
    if ( psfXml.sum() )
    {
        psf.sumY();
    }
    timer.stop();

    if ( psfXml.type() == DIC || psfXml.type() == DIC_2D )
    {
	  psf.rotateXY(psfXml.rotation());
	}
	if ( psfXml.realAndImag() ) 
	{
 
	  Array<T,3> psfArrayReal = psf.psfReal();
	  Array<T,3> psfArrayImag = psf.psfImag();

	  size_t i = psfXml.file().rfind('.', psfXml.file().length());
	  std::string psfprefix;
	  std::string psfsuffix;
      if ( i != std::string::npos )
      {
         psfprefix = psfXml.file().substr(0,i);
		 psfsuffix = psfXml.file().substr(i,  psfXml.file().length());
      }
      else
      {
         psfprefix = psfXml.file().substr(0,psfXml.file().length());
      }
	  TinyVector<int,3> shift((psf.psf().shape()+1)/2);
	  if ( psfXml.centered() )
	  {
		  Array<T,3> psfArrayRealShift = circularShift(psfArrayReal, shift);
		  Array<T,3> psfArrayImagShift = circularShift(psfArrayImag, shift);
          wuDataWrite(psfArrayRealShift, psfprefix + "_re" + psfsuffix);
	      wuDataWrite(psfArrayImagShift, psfprefix + "_im" + psfsuffix);
	  }
	  else
	  {
          wuDataWrite(psfArrayReal, psfprefix + "_re" + psfsuffix);
	      wuDataWrite(psfArrayImag, psfprefix + "_im" + psfsuffix);
      }
    }
    std::cout <<"filename is "<< psfXml.file().c_str()<< std::endl;
	Array<T,3> psfArray = psf.psf();
	if ( psfXml.centered() )
    {
        TinyVector<int,3> shift((psf.psf().shape()+1)/2);
		std::cout << "Centered "<< shift << std::endl;
        Array<T,3> psfArrayShift = circularShift(psfArray, shift);
		wuDataWrite(psfArrayShift, psfXml.file());
	}
	else
	{
	    wuDataWrite(psfArray, psfXml.file());
	}

    std::cout <<"PSF generation took "<< timer.elapsedSeconds() << " seconds "<< std::endl;
}

// ----------------------------------------------------------------------------
int
main (int argc, char* argv[])
{
    // Define command line object
    CmdLine cmdLine("COSM PSF Generation",' ', VERSION);

    // Define the argument options
    SwitchArg cliArg("c", "cli", "Use command line interface", false);
    ValueArg<std::string> xmlArg("f","file", "PSF XML file", false, "psf.xml", "psf xml");
    SwitchArg doubleArg("d", "double", "Use double precision", false);
    SwitchArg longdoubleArg("l", "long", "Use long double precision", false);

    //Version 0.8.11
    ValueArg<double> IntervalArg("i","interval", "Psf Thickness Interval", false, 0, "thickness interval");
    ValueArg<int> TotalNumArg("t","total", "Psf Total Number", false, 1, "total number");
    SwitchArg SumArg("s", "sum", "Sum over Y", false);
    SwitchArg UncenteredArg("u", "uncentered", "Uncentered psf file", false);
	SwitchArg realAndImagArg("g", "gri", "Generate Real and Imaginary", false);

    // add arguments to command line options
    cmdLine.add(cliArg);
    cmdLine.add(xmlArg);
    cmdLine.add(doubleArg);
    cmdLine.add(longdoubleArg);

    //Version 0.8.11
    cmdLine.add(IntervalArg);
    cmdLine.add(TotalNumArg);
    cmdLine.add(SumArg);
    cmdLine.add(UncenteredArg);
                                                                                
    // parse command line
    cmdLine.parse(argc, argv);

    PsfGUI gui;
    std::cout <<"gui()"<< std::endl;
    if ( cliArg.getValue() )
    {
        PsfUser user;
        if ( doubleArg.getValue() && longdoubleArg.getValue() ) 
        {
            PsfXml<long double> psfXml;
            if ( psfXml.open(xmlArg.getValue()) )
            {
                std::cout <<"Enter main()::long double"<< std::endl;
                evaluate<long double>(psfXml, &user);
            }
        }
        else if ( doubleArg.getValue() ) 
        {
            PsfXml<double> psfXml;
            if ( psfXml.open(xmlArg.getValue()) )
            {
                std::cout <<"Enter main()::double"<< std::endl;
                evaluate<double>(psfXml, &user);
            }
        } 
            else 
        {
            PsfXml<float> psfXml;
            if ( psfXml.open(xmlArg.getValue()) )
            {
                std::cout <<"Enter main()::float"<< std::endl;

                //Version 0.8.11
                //Read filename
                std::string psfname = psfXml.file();
                if (psfname == ""){
                	std::cout <<"No PSF Output Filename specified"<< std::endl;
                    return 1;
                }
                std::string psfsuffix;
                std::string psfprefix;
                size_t i = psfname.rfind('.', psfname.length());
                if ( i != std::string::npos )
                {
                    psfprefix = psfname.substr(0,i);
                    psfsuffix = psfname.substr(i, psfname.length());
                }

                //Read parameters from CLI input: sum & centered
                psfXml.sum(SumArg.isSet());
                psfXml.centered(!UncenteredArg.isSet());
				psfXml.realAndImag(realAndImagArg.isSet());

                //Read initial thicknesss
                float thickness = psfXml.ts();

                //processing evaluate for single or multiple psf(s)
                for ( int i = 0; i < TotalNumArg.getValue(); i++ )
                {
                    std::ostringstream psfindex;
                    psfindex << psfprefix << i;
                    std::string psffile = (TotalNumArg.getValue() == 1) ? psfname : (psfindex.str() + psfsuffix);
                    std::string xmlfile = (TotalNumArg.getValue() == 1) ? (psfprefix + ".xml") : (psfindex.str() + ".xml");
                    psfXml.file(psffile);
                    psfXml.ts(thickness);
                    psfXml.save(xmlfile);
                    evaluate<float>(psfXml, &user);
                    thickness += IntervalArg.getValue();
                }
            }
        }
    }
    else
    {
        Fl::scheme("plastic");
        gui.show(argc, argv);
        int fl_ret = Fl::run();
        return fl_ret;
    }

}
