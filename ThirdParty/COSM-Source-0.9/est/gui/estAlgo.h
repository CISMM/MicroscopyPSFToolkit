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
 
#ifndef _EST_ALGO_H
#define _EST_ALGO_H
 
#include "est/estimateIO.h"
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
#include "fltk/Cartesian.h"
#include "fltk/Table.h"
#include <FL/Fl.H>
#include <FL/Fl_Progress.H>

#include <blitz/array.h>
#include <blitz/timer.h>
#include <sstream>

extern bool gStop;

//        [<io><double>][<dept variant><iterative>][<algorithm variation>]
enum {
    EST_LIN    = 0x001,
    EST_ITER   = 0x010,
    EST_VAR    = 0x020,
    EST_DOUBLE = 0x100,
    EST_IO     = 0x200,
 
    EST_LLS    = 0x003,
    EST_MAP    = 0x005,
 
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
class ImageIO : public cosm::EstimateIO<T,N> {
public:
    ImageIO() {};
    virtual ~ImageIO() {};
    void ReadData( blitz::Array<T,N>& img, const std::string& name )
    { cosm::wuDataRead(img, name); };
    void WriteData( blitz::Array<T,N>& img, const std::string& name )
    { cosm::wuDataWrite(img, name); };
    void ReadData( blitz::Array<std::complex<T>,N>& img, const std::string& name )
    { cosm::wuDataRead(img, name); };
    void WriteData( blitz::Array<std::complex<T>,N>& img, const std::string& name )
    { cosm::wuDataWrite(img, name); };
};


template<typename T, int N>
class ErrorObserver : public cosm::EstimateObserver<T,N> {
  public:
    ErrorObserver( 
        unsigned int period, 
        unsigned short errorMask,
		Table* table 
    ) : cosm::EstimateObserver<T,N>(period), 
        errorMask_(errorMask), 
        errorEstimate_(cosm::MAXIMUM_ERROR),
        hasPhantom_(false),
		table_(table)
    {  for ( int i = 0; i < 5; i++ ) line_[i] = 0; };

    ~ErrorObserver() {};

    void setPhantom(
        blitz::Array<T,N>& phantom
    ) {
        phantom_.resize(phantom.extent()); 
        phantom_ = phantom; 
        hasPhantom_ = true; 
    }

  protected:
    virtual void update( cosm::EstimateIterative<T,N>& estimate )
    {
	    std::ostringstream ostr[6];
		ostr[0] << estimate.getIterationsDone();
		double x = (double) estimate.getIterationsDone();
        cout <<"iterations: "<< estimate.getIterationsDone() << ", error: "; 

        if ( errorMask_ & 0x1 ) {
            T error = errorEstimate_.error( 
                    cosm::MAXIMUM_ERROR,
                    estimate.results(), 
                    hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
			) ;
			cout << error << " ";
			ostr[1] << error;
			if ( table_ != NULL ) line_[0] = new Ca_LinePoint(line_[0], x, double(error), 0, FL_BLACK);
                
        }
        if ( errorMask_ & 0x2 ) {
            T error = errorEstimate_.error( 
                    cosm::MEAN_ERROR,
                    estimate.results(), 
                    hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
			);
			cout << error << " ";
			ostr[2] << error;
			if ( table_ != NULL ) line_[1] = new Ca_LinePoint(line_[1], x, double(error), 0, FL_BLUE);
        }
        if ( errorMask_ & 0x4 ) {
            T error = errorEstimate_.error( 
                    cosm::MEAN_SQUARE_ERROR,
                    estimate.results(), 
                    hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
			);
			cout << error << " ";
			ostr[3] << error;
			if ( table_ != NULL ) line_[2] = new Ca_LinePoint(line_[2], x, double(error), 0, FL_GREEN);
        }
        if ( errorMask_ & 0x8 ) {
            T error = errorEstimate_.error( 
                    cosm::LOG_LIKELYHOOD_ERROR,
                    estimate.results(), 
                    hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
			);
			cout << error << " ";
			ostr[4] << error;
			if ( table_ != NULL ) line_[3] = new Ca_LinePoint(line_[3], x, double(error), 0, FL_RED);
        }
        if ( errorMask_ & 0x10 ) {
            T error = errorEstimate_.error( 
                    cosm::I_DIVERGENCE_ERROR,
                    estimate.results(), 
                    hasPhantom_ ? phantom_ : estimate.getOldEstimate() 
			);
			cout << error << " ";
			ostr[5] << error;
			if ( table_ != NULL ) line_[4] = new Ca_LinePoint(line_[4], x, double(error), 0, FL_YELLOW);
        }
        cout <<endl;
		if ( table_ != NULL )
		{
		    table_->addRow(6, 
			    ostr[0].str().c_str(),
			    ostr[1].str().c_str(), 
				ostr[2].str().c_str(), 
				ostr[3].str().c_str(), 
				ostr[4].str().c_str(),
				ostr[5].str().c_str());
			Fl::wait(0.0);
		}
    }

  private:

    unsigned short errorMask_;
    cosm::ErrorEstimate<T,N> errorEstimate_;
    blitz::Array<T, N> phantom_;
    bool hasPhantom_;
	Table* table_;
	Ca_LinePoint* line_[5];
};

template<typename T, int N>
class OutputObserver : public cosm::EstimateObserver<T,N> {
  public:
    OutputObserver(
        unsigned int period, 
        const std::string& prefix,
        const std::string& suffix
    ) : cosm::EstimateObserver<T,N>(period), prefix_(prefix), suffix_(suffix) { };

    ~OutputObserver() { };

  protected:
    virtual void update( cosm::EstimateIterative<T,N>& estimate )
    { 
        char number[9];
        sprintf(number, "%06d", estimate.getIterationsDone());
        std::string filename = prefix_ + "_" + number + "." + suffix_; 
        cosm::wuDataWrite(estimate.results(), filename); 
    };
 
  private:
    std::string prefix_;
    std::string suffix_;
};

template<typename T, int N>
class ProgressObserver : public cosm::EstimateObserver<T,N> {
  public:   
    ProgressObserver( 
        unsigned int period, 
        Fl_Progress* progress
    ) : cosm::EstimateObserver<T,N>(period), progress_(progress) { };

    ~ProgressObserver() { };

  protected:
    virtual void update( cosm::EstimateIterative<T,N>& estimate )
    { 
        // update progress bar
        if ( progress_ )
        {
            progress_->value(float(100.0)*float(estimate.getIterationsDone())/float(estimate.getIterations()));
            Fl::wait(0.0);
        }
    };
  protected:
    Fl_Progress* progress_;
};
   

template<typename T, int N>
class CtrlObserver : public cosm::EstimateObserver<T,N> {
  public:
    CtrlObserver( unsigned int period ) : cosm::EstimateObserver<T,N>(period) {};

    ~CtrlObserver() { };

  protected:
    void update( cosm::EstimateIterative<T,N>& estimate ) 
    {
        if ( gStop )
        {
	    estimate.abort();
        }
        Fl::wait(0.0);
    }
};


template<typename T, int N>
blitz::Array<T,N> performEstimation(
    int algo,
    blitz::Array<T,N>& img,
    blitz::Array<T,N>& psf1,
    blitz::Array<T,N>& phantom,
    bool usePhantom,
    int iterations,
    int update,
    int writeUpdate,
    int numberOfStrata,
    int startOfStrata,
    int sizeOfStrata,
    T value,
    unsigned short err,
    const std::string& psfprefix,
    const std::string& psfsuffix,
    const std::string& otfprefix,
    const std::string& otfsuffix,
    const std::string& estprefix,
    const std::string& estsuffix,
	bool centeredPsf = true,
    cosm::EstimatePenalty<T,N>* penalty = 0,
    Fl_Progress* prog = 0,
	Table* table = 0
) {
    // compute strata and psfs
	TinyVector<int, N> shift = (psf1.length()+1)/2;
	std::cout <<"shift: "<< shift << std::endl;
    blitz::Array<RectDomain<N>, 1> strata(numberOfStrata+2);
    blitz::Array<blitz::Array<T,N>, 1> psfs(numberOfStrata+1);
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
                    psfm << psfprefix << m << "." << psfsuffix;
                    std::cout << "Reading: "<< psfm.str() << std::endl;
                    cosm::wuDataRead(psf1, psfm.str());
                    TinyVector<int, N> shift = (psf1.length()+1)/2;
                    std::cout <<"shift: "<< shift << std::endl;
					psfs(m).resize(psf1.extent()); 
					psfs(m) = (!centeredPsf ? psf1 : circularShift(psf1, shift));
                }
            }
        }
    }
	blitz::Array<T,N> psf = (!centeredPsf ? psf1 : circularShift(psf1, shift));
	
    // create observer objects
    std::vector<cosm::EstimateObserver<T,N>* > observers;
    cosm::EstimateObserver<T,N>* observer = 0;
    observer = new OutputObserver<T,N>( writeUpdate,  estprefix, estsuffix);
    observers.push_back(observer);
    observer = new ErrorObserver<T,N>(update, (cosm::ErrorType)err, table);
    if ( usePhantom )
    {
        ((ErrorObserver<T,N>*)observer)->setPhantom(phantom);
    }
    observers.push_back(observer);

    if ( prog != NULL )
    {
        observer = new ProgressObserver<T,N>(update, prog);
        observers.push_back(observer);
        observer = new CtrlObserver<T,N>(update);
        observers.push_back(observer);
    }
 
    // create image io object for variant algorithms
    ImageIO<T, N>* imageIO = new ImageIO<T,N>();
 
    // start estimation;
    cosm::Estimate<T, N>* estimate = NULL;
 
    switch( algo )
    {
        case EST_LLS:
        {
            estimate = new cosm::EstimateLLS<T,N>(psf,img, value);
            break;
        }
        case EST_MAP:
        {
            estimate = new cosm::EstimateMAP<T,N>(psf,img, value);
            break;
        }
        case EST_JVC:
        {
            estimate = new cosm::EstimateJVC<T,N>(psf, img, iterations);            break;
        }
        case EST_EM:
        {
            estimate = new cosm::EstimateEM<T,N>(psf, img, iterations, penalty);
            break;
        }
        case EST_EM2:
        {
            estimate = new cosm::EstimateEM2<T,N>(psf, img, iterations, penalty);            break;
        }
        case EST_EMSV:
        {
            estimate = new cosm::EstimateEMSV<T, N>(strata, psfs, img, iterations);
            break;
        }
        case EST_EMSV2:
        {
            estimate = new cosm::EstimateEMSV2<T, N>(strata, psfs, img, iterations);
            break;
        }
        case EST_EMOS:
        {
            estimate = new cosm::EstimateEMOS<T, N>(strata, psfs, img, iterations);
            break;
        }
        case EST_EMOS2:
        {
            estimate = new cosm::EstimateEMOS2<T, N>(strata, psfs, img, iterations);
            break;
        }
        case EST_EMSVB:
        {
            estimate = new cosm::EstimateEMSVbig<T, N>(strata, psfprefix, otfprefix, psfsuffix, img, iterations, imageIO);
            break;
        }
        case EST_EMSV2B:
        {
            estimate = new cosm::EstimateEMSV2big<T, N>(strata, psfprefix, otfprefix, psfsuffix, img, iterations, imageIO);
            break;
        }
        case EST_EMOSB:
        {
            estimate = new cosm::EstimateEMOSbig<T, N>(strata, psfprefix, otfprefix, psfsuffix, img, iterations, imageIO );
            break;
        }
        case EST_EMOS2B:
        {
            estimate = new cosm::EstimateEMOS2big<T, N>(strata, psfprefix,  otfprefix, psfsuffix, img, iterations, imageIO);
            break;
        }
    }
    if ( estimate == 0)
    {
        cout <<"Unknown estimate : "<< hex << algo << dec << endl;
        blitz::Array<T,N> res;
        return res;
    }
    if( algo & EST_ITER )
    {
        for ( unsigned int i = 0; i < observers.size(); i++ )
        {
            ((cosm::EstimateIterative<T,N>*)estimate)->registerObserver(observers[i]);
        }
    }
    blitz::Timer timer;
    timer.start();
    estimate->run();
    timer.stop();
    cout <<"Estimation completed in "<< timer.elapsedSeconds() << endl;
    blitz::Array<T,N> res = estimate->results();

    delete estimate;
    delete imageIO;
    return res;
}


#endif // _EST_ALGO_H
