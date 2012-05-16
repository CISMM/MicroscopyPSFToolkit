#include "estGUI.h"
#include "estAlgo.h"
#include <FL/Fl.H>

#include "wu/wuHeader.h"
#include "wu/wuImage.h"
#include <string>
#include <iostream>
#include <tclap/CmdLine.h>
#include <blitz/timer.h>

#define XSTR(s) STR(s)
#define STR(s) #s
#define VERSION XSTR(COSM_VERSION)

using namespace blitz;
using namespace TCLAP;
using namespace cosm;

const int N = 3;

// ----------------------------------------------------------------------------
int
main (int argc, char* argv[])
{

    // Define command line object
    CmdLine cmdLine("COSM Estimation", ' ', VERSION);

    SwitchArg cliArg("c", "cli", "Command Line Interface", false);

    // Define the argument options
    ValueArg<std::string> estArg("t", "est", "Estimate filename prefix", false, "est", "estimate prefix");
    ValueArg<std::string> imgArg("i", "img", "Image filename prefix", false, "img", "image prefix");
    ValueArg<std::string> psfArg("p", "psf", "PSF filename prefix", false, "psf", "psf prefix");
    ValueArg<std::string> otfArg("f", "otf", "OTF filename prefix", false, "otf", "otf prefix");
    ValueArg<std::string> phaArg("y", "pha", "Phantom filename prefix", false, "pha", "phantom prefix");
    ValueArg<std::string> suffixArg("x", "suffix", "Filename suffix", false, "wu", "suffix"); //Version 0.8.11::DEBUG (".wu" -> "wu")
    ValueArg<int> updateArg("u", "update", "Number of iterations between statistics updates", false, 100, "update");
    ValueArg<int> writeArg("w", "write", "Number of iterations between write updates", false, 100, "write update");
    ValueArg<double> llsArg("l", "lls", "LLS estimate", false, 1E-4, "pval");
    ValueArg<double> mapArg("m", "map", "MAP estimate", false, 1E-4, "alpha");
    ValueArg<int> emArg("e", "em", "EM estimate", false, 100, "iterations");
    ValueArg<int> jvcArg("j", "jvc", "JVC estimate", false, 100, "iterations");
    SwitchArg svArg("v", "sv", "EM-SV estimate", false);
    SwitchArg osArg("o", "os", "EM-OS estimate", false);
    ValueArg<int> numStrataArg("n", "num", "Number of strata", false, 1, "stratanumber");
    ValueArg<int> startStrataArg("s", "start", "Start of strata", false, 1, "stratastart");
    ValueArg<int> sizeStrataArg("k", "size", "Size of strata", false, 1, "stratasize");
    SwitchArg doubleArg("d", "double", "Use double Z dimension", false);
    ValueArg<double> intensityArg("I", "intensity", "Intensity penalty, 0 =< penalty =< 1", false, 0, "penalty");
    ValueArg<double> roughnessArg("R", "roughness", "Roughness penalty, 0 =< penalty =< 1", false, 0, "penalty");
    ValueArg<unsigned short> errorArg("E", "err", "Error estimate", false, 0x1f, "error");

    //Option of centerdPsf choice
    ValueArg<bool> centeredPsfArg("C", "center", "Centered Psf file", false, true, "centered psf");     //Version 0.8.11::Improve-5(CLI)
   
    // Add arguments to command line options
    cmdLine.add(cliArg);
    cmdLine.add(estArg);
    cmdLine.add(imgArg);
    cmdLine.add(psfArg);
    cmdLine.add(otfArg);
    cmdLine.add(phaArg);
    cmdLine.add(suffixArg);

/*
    vector<Arg*> xorArg(4);
    xorArg[0] = &llsArg;
    xorArg[1] = &mapArg;
    xorArg[2] = &emArg;
    xorArg[3] = &jvcArg;
    cmdLine.xorAdd(xorArg);
*/
    cmdLine.add(llsArg);
    cmdLine.add(mapArg);
    cmdLine.add(emArg);
    cmdLine.add(jvcArg);

    cmdLine.add(svArg);
    cmdLine.add(osArg);
    cmdLine.add(doubleArg);
    cmdLine.add(intensityArg);
    cmdLine.add(roughnessArg);
    cmdLine.add(updateArg);
    cmdLine.add(writeArg);
    cmdLine.add(errorArg);

    cmdLine.add(centeredPsfArg);    //Version 0.8.11::Improve-5(CLI)
	cmdLine.add(numStrataArg);      //Version 0.8.11::DEBUG(cmdLine.Add Strata)
    cmdLine.add(startStrataArg);    //Version 0.8.11::DEBUG(cmdLine.Add Strata)
	cmdLine.add(sizeStrataArg);     //Version 0.8.11::DEBUG(cmdLine.Add Strata)
        
    // Parse the command line
    cmdLine.parse(argc, argv);

    if ( cliArg.getValue() )
    {

        std::string estname = estArg.getValue();
        std::string psfname = psfArg.getValue();
        std::string otfname = otfArg.getValue();
        std::string imgname = imgArg.getValue();
        std::string phaname = phaArg.getValue();
        std::string suffix = suffixArg.getValue();

	    std::ostringstream suffix2; //Version 0.8.11::DEBUG (".wu" -> "wu")
	    suffix2 << "." << suffix;   //Version 0.8.11::DEBUG (".wu" -> "wu")

        // read in img file
        WUImage imgData;
        if ( imgData.ReadData(imgname+suffix2.str()) == false ) //Version 0.8.11::DEBUG (".wu" -> "wu")
        {
            std::cout << "Reading image data failed"<< std::endl;
            return -1;
        };  

        // read phantom file
        bool usePhantom = false;
        WUImage phantomData;
        if ( phaArg.isSet() ) 
        {
            if ( phantomData.ReadData(phaname+suffix2.str()) == false ) //Version 0.8.11::DEBUG (".wu" -> "wu")
            {
                std::cout << "Reading image phantom failed" << std::endl;
                return -1;
            }
            usePhantom = true;
        }

        // read in psf file
        WUImage psfData;
        if ( !svArg.isSet() && !osArg.isSet() ) 
        {
            if (psfData.ReadData(psfname+suffix2.str()) == false ) //Version 0.8.11::DEBUG (".wu" -> "wu")
            {
                std::cout << "Reading psf failed" << std::endl;
                return -1;
            }
        }
        // command line parameters
        int iterations = jvcArg.isSet() ? jvcArg.getValue() : emArg.getValue();
        int update = updateArg.getValue();
        int writeUpdate = writeArg.getValue();
        unsigned short err = errorArg.getValue();
        double value = llsArg.isSet() ? llsArg.getValue() : mapArg.getValue();
        bool intensity = intensityArg.isSet();
        bool roughness = roughnessArg.isSet();
		if ( jvcArg.isSet() || emArg.isSet() || svArg.isSet() || osArg.isSet() ) //Version 0.8.11::DEBUG (prevent 'value' override)
		{                                                                        //Version 0.8.11::DEBUG (prevent 'value' override)
			value = intensity ? intensityArg.getValue() : 0;
			value = roughness ? roughnessArg.getValue() : 0;
		}                                                                        //Version 0.8.11::DEBUG (prevent 'value' override)

        int numberOfStrata = numStrataArg.getValue();
        int startOfStrata = startStrataArg.getValue();
        int sizeOfStrata = sizeStrataArg.getValue();

        bool centeredPsf = centeredPsfArg.getValue(); //Version 0.8.11::Improve-5(CLI)


        // determin the algorithm
        int algo = 0;
        if ( llsArg.isSet() )
        {
            algo = EST_LLS;
        } 
        else if ( mapArg.isSet() )
        {
            algo = EST_MAP;
        }
        else if ( jvcArg.isSet() )
        {
            algo = EST_JVC;
        }
        else if ( emArg.isSet() )
        {
		// emArg.isSet() means that whether either EM, SV or OS be checked,
	    // and jvcArg.isset() above only means that whether JVC be checked.
            if ( svArg.isSet() || osArg.isSet() )    //Version 0.8.11::DEBUG (algo 133->131)
            {                                        //Version 0.8.11::DEBUG (algo 133->131)
                if ( svArg.isSet() )                 //Version 0.8.11::DEBUG (algo 133->131)
                {                                    //Version 0.8.11::DEBUG (algo 133->131)
                    algo = EST_EMSV;                 //Version 0.8.11::DEBUG (algo 133->131)
                }                                    //Version 0.8.11::DEBUG (algo 133->131)
                else if ( osArg.isSet() )            //Version 0.8.11::DEBUG (algo 133->131)
                {                                    //Version 0.8.11::DEBUG (algo 133->131)
                    algo = EST_EMOS;                 //Version 0.8.11::DEBUG (algo 133->131)
                }                                    //Version 0.8.11::DEBUG (algo 133->131)
                if ( otfArg.isSet() )                //Version 0.8.11::DEBUG (algo 133->131)
                {                                    //Version 0.8.11::DEBUG (algo 133->131)
                    algo |= EST_IO;                  //Version 0.8.11::DEBUG (algo 133->131)
                }                                    //Version 0.8.11::DEBUG (algo 133->131)
			}else {                                  //Version 0.8.11::DEBUG (algo 133->131)
				algo = EST_EM;                       //Version 0.8.11::DEBUG (algo 133->131)
			}                                       
  
			if ( doubleArg.isSet() )                 //Version 0.8.11::DEBUG (algo 133->131)
            {                                        //Version 0.8.11::DEBUG (algo 133->131)
                algo |= EST_DOUBLE;                  //Version 0.8.11::DEBUG (algo 133->131)
            }                                        //Version 0.8.11::DEBUG (algo 133->131)
        }

        // convert to same data type and do estimation
        if ( imgData.IsFloat() )
        {
            Array<float,N> img = imgData.GetFloatArray();
            psfData.ConvertToFloat();
            Array<float,N> psf = psfData.GetFloatArray();
            if ( !larger(img.shape(), psf.shape()) )
            {
                std::cout<<"PSF dimensions larger than image"<< std::endl;
                return -1;
            }
            phantomData.ConvertToFloat();
            Array<float, N> phantom = phantomData.GetFloatArray();
            EstimatePenalty<float,N>* penalty = 
                (intensity ? (cosm::EstimatePenalty<float,N>*)(new cosm::IntensityPenalty<float,N>(value)) :
                (roughness ? (cosm::EstimatePenalty<float,N>*)(new cosm::RoughnessPenalty<float,N>(value)) : NULL));

            Array<float,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, writeUpdate, numberOfStrata,       //Version 0.8.11::Improve-5(CLI)
            		startOfStrata, sizeOfStrata, (float)value, err, psfname, suffix, otfname, suffix, estname, suffix, centeredPsf, penalty);  //Version 0.8.11::Improve-5(CLI)
					wuDataWrite(est, estname + suffix2.str()); //Version 0.8.11::DEBUG (".wu" -> "wu")
        } 
        else if ( imgData.IsDouble() )
        {
            Array<double,N> img = imgData.GetDoubleArray();
            psfData.ConvertToDouble();
            Array<double,N> psf = psfData.GetDoubleArray();
            if ( !larger(img.shape(), psf.shape()) )
            {
                std::cout<<"PSF dimensions larger than image"<< std::endl;
                return -1;
            }
            phantomData.ConvertToDouble();
            Array<double, N> phantom = phantomData.GetDoubleArray();
            EstimatePenalty<double,N>* penalty = 
                (intensity ? (cosm::EstimatePenalty<double,N>*)(new cosm::IntensityPenalty<float,N>(value)) :
                (roughness ? (cosm::EstimatePenalty<double,N>*)(new cosm::RoughnessPenalty<float,N>(value)) : NULL));
            Array<double,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, writeUpdate, numberOfStrata,       //Version 0.8.11::Improve-5(CLI)
            		startOfStrata, sizeOfStrata, (double)value, err, psfname, suffix, otfname, suffix, estname, suffix, centeredPsf, penalty);  //Version 0.8.11::Improve-5(CLI)
            wuDataWrite(est, estname + suffix2.str()); //Version 0.8.11::DEBUG (".wu" -> "wu")
        }
        else if ( imgData.IsLongDouble() )
        {
            Array<long double,N> img = imgData.GetLongDoubleArray();
            psfData.ConvertToLongDouble();
            Array<long double,N> psf = psfData.GetLongDoubleArray();
            if ( !larger(img.shape(), psf.shape()) )
            {
                std::cout<<"PSF dimensions larger than image"<< std::endl;
                return -1;
            }
            phantomData.ConvertToLongDouble();
            Array<long double, N> phantom = phantomData.GetLongDoubleArray();
            EstimatePenalty<long double,N>* penalty = 
                (intensity ? (cosm::EstimatePenalty<long double,N>*)(new cosm::IntensityPenalty<float,N>(value)) :
                (roughness ? (cosm::EstimatePenalty<long double,N>*)(new cosm::RoughnessPenalty<float,N>(value)) : NULL));
            Array<long double,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, writeUpdate, numberOfStrata,        //Version 0.8.11::Improve-5(CLI)
            		startOfStrata, sizeOfStrata, (long double) value, err, psfname, suffix, otfname, suffix, estname, suffix, centeredPsf, penalty);  //Version 0.8.11::Improve-5(CLI)
            wuDataWrite(est, estname + suffix2.str()); //Version 0.8.11::DEBUG (".wu" -> "wu")
        }
    }
    else 
    { 
        EstGUI gui;
        Fl::scheme("plastic");
        gui.show(argc, argv);
        int fl_ret = Fl::run();
        return fl_ret;
    }
}
