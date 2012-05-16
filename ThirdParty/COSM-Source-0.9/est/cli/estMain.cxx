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

#include "est/estimateLLS.h"
#include "est/estimateMAP.h"
#include "est/estimateEM.h"
#include "est/estimateEM2.h"
#include "est/estimateJVC.h"
#include "est/estimateEMSV.h"
#include "est/estimateEMSV2.h"
#include "est/estimateEMSVbig.h"
#include "est/estimateEMSV2big.h"
#include "est/estimateEMOS.h"
#include "est/estimateEMOS2.h"
#include "est/estimateEMOSbig.h"
#include "est/estimateEMOS2big.h"
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

//        [<io><double>][<dept variant><iterative>][<algorithm variation>]
enum {
    EST_ITER   = 0x010,
    EST_VAR    = 0x020,
    EST_DOUBLE = 0x100,
    EST_IO     = 0x200,

    EST_LLS    = 0x001,
    EST_MAP    = 0x002,

    EST_JVC    = 0x011,
    EST_EM     = 0x012,
    EST_EM2    = 0x112,

    EST_EMSV   = 0x031,
    EST_EMSV2  = 0x131,
    EST_EMSVB  = 0x231,
    EST_EMSV2B = 0x331,
    EST_EMOS   = 0x032,
    EST_EMOS2  = 0x132,
    EST_EMOSB  = 0x232,
    EST_EMOS2B = 0x332
};

template<typename T, int N>
class ImageIO : public EstimateIO<T,N> {
public:
    ImageIO() {};
    ~ImageIO() {};
    void ReadData( Array<T,N>& img, const std::string& name ) 
    { wuDataRead(img, name); };
    void WriteData( Array<T,N>& img, const std::string& name ) 
    { wuDataWrite(img, name); };
    void ReadData( Array<std::complex<T>,N>& img, const std::string& name ) 
    { wuDataRead(img, name); };
    void WriteData( Array<std::complex<T>,N>& img, const std::string& name ) 
    { wuDataWrite(img, name); };
};

template<typename T, int N>
class TestUser : public EstimateUser<T,N> {

  public:

    TestUser( unsigned short errorMask ) : EstimateUser<T,N>(), hasPhantom_(false), errorEstimate_(MAXIMUM_ERROR), errorMask_(errorMask) {};
    TestUser(Array<T,N>& phantom, unsigned short errorMask) : EstimateUser<T,N>(), hasPhantom_(true), errorEstimate_(MAXIMUM_ERROR), phantom_(phantom), errorMask_(errorMask) {};
    void setPhantom(Array<T,N>& phantom) { phantom_.resize(phantom.extent()); phantom_ = phantom; hasPhantom_ = true; }

    ~TestUser() {};

    virtual void update( EstimateIterative<T,N>& estimate )
    {
        cout <<"iterations: "<< estimate.getIterationsDone() <<", error: "; 

        if ( errorMask_ & 0x1 ) {
	    cout <<
	        errorEstimate_.error( 
		    MAXIMUM_ERROR,
	            estimate.results(), 
	            hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
	         ) << " ";
		
	}
        if ( errorMask_ & 0x2 ) {
	    cout <<
	        errorEstimate_.error( 
		    MEAN_ERROR,
	            estimate.results(), 
	            hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
	         ) << " ";
	}
        if ( errorMask_ & 0x4 ) {
	    cout <<
	        errorEstimate_.error( 
		    MEAN_SQUARE_ERROR,
	            estimate.results(), 
	            hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
	         ) << " ";
	}

	cout <<endl;
    }

  private:

    bool hasPhantom_;
    ErrorEstimate<T,N> errorEstimate_;
    Array<T, N> phantom_;
    unsigned short errorMask_;
};

template<typename T, int N>
Array<T,N> performEstimation(
    int algo,
    Array<T,N>& img, 
    Array<T,N>& psf, 
    Array<T,N>& phantom, 
    bool usePhantom,
    int iterations, 
    int update, 
    int numberOfStrata, 
    int startOfStrata, 
    int sizerOfStrata, 
    T value, 
    unsigned short err,
    const std::string& psfname,
    const std::string& otfname,
    const std::string& suffix
);

int main( 
    int argc, 
    char* argv[] 
) {

    // Define command line object
    CmdLine cmdLine("COSM Estimation", ' ', VERSION);

    // Define the argument options
    ValueArg<std::string> estArg("t", "est", "Estimate filename prefix", true, "est", "estimate prefix");
    ValueArg<std::string> imgArg("i", "img", "Image filename prefix", true, "img", "image prefix");
    ValueArg<std::string> psfArg("p", "psf", "PSF filename prefix", true, "psf", "psf prefix");
    ValueArg<std::string> otfArg("f", "otf", "OTF filename prefix", false, "otf", "otf prefix");
    ValueArg<std::string> phaArg("a", "pha", "Phantom filename prefix", false, "pha", "phantom prefix");
    ValueArg<std::string> suffixArg("s", "suffix", "Filename suffix", false, ".wu", "suffix");
    ValueArg<int> updateArg("u", "update", "Number of iterations between statistics updates", false, 100, "update");
    ValueArg<double> llsArg("l", "lls", "LLS estimate", true, 1E-4, "pval");
    ValueArg<double> mapArg("m", "map", "MAP estimate", true, 1E-4, "alpha");
    ValueArg<int> emArg("e", "em", "EM estimate", true, 100, "iterations");
    ValueArg<int> jvcArg("j", "jvc", "JVC estimate", true, 100, "iterations");
    SwitchArg svArg("v", "sv", "EM-SV estimate", false);
    SwitchArg osArg("o", "os", "EM-OS estimate", false);
    ValueArg<int> numStrataArg("n", "num", "Number of strata", false, 1, "stratanumber");
    ValueArg<int> startStrataArg("s", "start", "Start of strata", false, 1, "stratastart");
    ValueArg<int> sizeStrataArg("k", "size", "Size of strata", false, 1, "stratasize");
    SwitchArg doubleArg("d", "double", "Use double Z dimension", false);
    ValueArg<unsigned short> errorArg("r", "err", "Error estimate", false, 0x1, "error");
    ValueArg<double> alphaArg("a", "alpha", "Intensity Penalty parameter 0 =< alpha =< 1", false, 0, "alpha");
   
    // Add arguments to command line options
    cmdLine.add(estArg);
    cmdLine.add(imgArg);
    cmdLine.add(psfArg);
    cmdLine.add(otfArg);
    cmdLine.add(phaArg);
    cmdLine.add(suffixArg);

    vector<Arg*> xorArg(4);
    xorArg[0] = &llsArg;
    xorArg[1] = &mapArg;
    xorArg[2] = &emArg;
    xorArg[3] = &jvcArg;
    cmdLine.xorAdd(xorArg);

    cmdLine.add(svArg);
    cmdLine.add(osArg);
    cmdLine.add(doubleArg);
    cmdLine.add(updateArg);
    cmdLine.add(errorArg);
    cmdLine.add(alphaArg);
	
    // Parse the command line
    cmdLine.parse(argc, argv);
 
    std::string estname = estArg.getValue();
    std::string psfname = psfArg.getValue();
    std::string otfname = otfArg.getValue();
    std::string imgname = imgArg.getValue();
    std::string phaname = phaArg.getValue();
    std::string suffix = suffixArg.getValue();

    // read in img file
    WUImage imgData;
    if ( imgData.ReadData(imgname+suffix) == false )
    {
        std::cout << "Reading image data failed"<< std::endl;
        return -1;
    };  

    // read phantom file
    bool usePhantom = false;
    WUImage phantomData;
    if ( phaArg.isSet() ) 
    {
        if ( phantomData.ReadData(phaname+suffix) == false )
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
        if (psfData.ReadData(psfname+suffix) == false ) 
        {
            std::cout << "Reading psf failed" << std::endl;
            return -1;
        }
    }

    // command line parameters
    int iterations = jvcArg.isSet() ? jvcArg.getValue() : emArg.getValue();
    int update = updateArg.getValue();
    unsigned short err = errorArg.getValue();
    double value = llsArg.isSet() ? llsArg.getValue() : mapArg.getValue();
    value = (emArg.isSet() || jvcArg.isSet() || svArg.isSet() || osArg.isSet()) ? alphaArg.getValue() : 0;

    int numberOfStrata = numStrataArg.getValue();
    int startOfStrata = startStrataArg.getValue();
    int sizeOfStrata = sizeStrataArg.getValue();

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
        algo = EST_EM;

        if ( doubleArg.isSet() )
        {
            algo |= EST_DOUBLE;
        }
        if ( svArg.isSet() || osArg.isSet() )
        {
            if ( svArg.isSet() )
            {
                algo |= EST_EMSV;
            }
            else if ( osArg.isSet() )
            {
                algo |= EST_EMOS;
            }
            if ( otfArg.isSet() )
            {
                algo |= EST_IO;
            }
        }
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
        Array<float,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, numberOfStrata, startOfStrata, sizeOfStrata, (float)value, err, psfname, otfname, suffix);
        wuDataWrite(est, estname + suffix);
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
        Array<double,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, numberOfStrata, startOfStrata, sizeOfStrata, (double)value, err, psfname, otfname, suffix);
        wuDataWrite(est, estname + suffix);
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
        Array<long double,N> est = performEstimation(algo, img, psf, phantom, usePhantom, iterations, update, numberOfStrata, startOfStrata, sizeOfStrata, (long double)value, err, psfname, otfname, suffix);
        wuDataWrite(est, estname + suffix);
    }
}


template<typename T, int N>
Array<T,N> performEstimation(
    int algo,
    Array<T,N>& img, 
    Array<T,N>& psf, 
    Array<T,N>& phantom, 
    bool usePhantom,
    int iterations, 
    int update, 
    int numberOfStrata, 
    int startOfStrata, 
    int sizeOfStrata, 
    T value, 
    unsigned short err,
    const std::string& psfname,
    const std::string& otfname,
    const std::string& suffix
) {

    // compute strata and psfs
    Array<RectDomain<N>, 1> strata(numberOfStrata);
    Array<Array<T,N>, 1> psfs(numberOfStrata+1);
    if ( algo & EST_VAR )
    {
        for ( int m = 0; m <= numberOfStrata+1; m++ )
        {
            TinyVector<int, N> l(img.lbound());
            TinyVector<int,N> u(img.ubound());
            if( m == 0 )
            {
                u(0) = startOfStrata - 1;
            }
            else if ( m <= numberOfStrata )
            {
                l(0) = startOfStrata + sizeOfStrata * (m-1);
                u(0) = startOfStrata + sizeOfStrata * m - 1;
            }
            else
            {
                l(0) = startOfStrata + sizeOfStrata * (m-1);
            }
            strata(m).setlbound(l);
            strata(m).setubound(u);
            std::cout <<"Strata: "<<m<<" "<<l<<", "<<u<<std::endl;

            if ( m <= numberOfStrata )
            {
                if ( (algo & EST_IO) != EST_IO )
                {
                    std::ostringstream psfm;
                    psfm << psfname << m << suffix;
                    std::cout <<psfm.str()<<std::endl;
                    cosm::wuDataRead(psfs(m),psfm.str());
                }
            }
        }
    }

    // create user object
    TestUser<T,N>* user = new TestUser<T,N>(err);
    if ( usePhantom )
    {
        user->setPhantom(phantom);
    }

    // create image io object for variant algorithms
    ImageIO<T, N>* imageIO = new ImageIO<T,N>();

    // start estimation;
    Estimate<T, N>* estimate = NULL;

    switch( algo )
    {
        case EST_LLS: 
        {
	    estimate = new EstimateLLS<T,N>(psf,img, value);
            break;
        }
        case EST_MAP:
        {
	    estimate = new EstimateMAP<T,N>(psf,img, value);
            break;
        }
        case EST_JVC:
        {
	    estimate = new EstimateJVC<T,N>(psf, img, iterations, user, update, value);
            break;
        }
        case EST_EM:
        {
	    estimate = new EstimateEM<T,N>(psf, img, iterations, user, update, value);
            break;
        }
        case EST_EM2:
        {
	    estimate = new EstimateEM2<T,N>(psf, img, iterations, user, update, value);
            break;
        }
        case EST_EMSV:
        {
	    estimate = new EstimateEMSV<T, N>(strata, psfs, img, iterations, user, update);
            break;
        }
        case EST_EMSV2:
        {
	    estimate = new EstimateEMSV2<T, N>(strata, psfs, img, iterations, user, update);
            break;
        }
        case EST_EMOS:
        {
	    estimate = new EstimateEMOS<T, N>(strata, psfs, img, iterations, user, update);
            break;
        }
        case EST_EMOS2:
        {
	    estimate = new EstimateEMOS2<T, N>(strata, psfs, img, iterations, user, update);
            break;
        }
        case EST_EMSVB:
        {
	    estimate = new EstimateEMSVbig<T, N>(strata, psfname, otfname, suffix, img, iterations, imageIO, user, update);
            break;
        }
        case EST_EMSV2B:
        {
	    estimate = new EstimateEMSV2big<T, N>(strata, psfname, otfname, suffix, img, iterations, imageIO, user, update);
            break;
        }
        case EST_EMOSB:
        {
	    estimate = new EstimateEMOSbig<T, N>(strata, psfname, otfname, suffix, img, iterations, imageIO, user, update);
            break;
        }
        case EST_EMOS2B:
        {
	    estimate = new EstimateEMOS2big<T, N>(strata, psfname, otfname, suffix, img, iterations, imageIO, user, update);
            break;
        }
    }
    Timer timer;
    timer.start();
    estimate->run();
    timer.stop();
    cout <<"Estimation completed in "<< timer.elapsedSeconds() << endl;
    Array<T,N> res = estimate->results();

    delete estimate;
    delete imageIO;
    delete user;

    return res;
}
