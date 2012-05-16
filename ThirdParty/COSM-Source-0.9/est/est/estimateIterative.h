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

#ifndef _ESTIMATE_ITERATIVE_H
#define _ESTIMATE_ITERATIVE_H

#include "estimate.h"
#include "estimateObserver.h"
#include "estimatePenalty.h"
#include <blitz/tinyvec-et.h>
#include <vector>

namespace cosm {

template<typename T, int N>
class EstimateIterative : public Estimate<T,N> {

  public:

    // class constructor
    EstimateIterative(
	Array<T,N>& psf,
	Array<T,N>& img,
	int iterations,
        EstimatePenalty<T,N>* penalty = NULL
    ) : Estimate<T,N>(psf,img), 
	iterations_(iterations), 
	iterationsDone_(0), 
	penalty_(penalty),
	running_(true)
    {   
        old_.resize(this->est_.extent()); 
        this->est_ = 1;
    };

    // class constructor
    EstimateIterative(
	Array<T,N>& img,
	int iterations,
        EstimatePenalty<T,N>* penalty = NULL
    ) : Estimate<T,N>(img), 
	iterations_(iterations), 
	iterationsDone_(0), 
	penalty_(penalty),
	running_(true)
    {   
        old_.resize(this->est_.extent()); 
        this->est_ = 1;
    };

    // emEstimate class destructor
    virtual ~EstimateIterative() { };
 
    // Function to register observers that will be notified each iteration
    void registerObserver( EstimateObserver<T,N>* observer ) 
    { observers_.push_back(observer); };

    // Function to run the estimation Algorithm for a number of iterations
    virtual int run()
    {
        while ( running_ && (iterationsDone_ < iterations_) )
        {
            for ( unsigned int i = 0; i < observers_.size(); i++ )
            {
                EstimateObserver<T,N>* observer = observers_[i];
                if ( observer != NULL )
                {
                    observer->notify(*this);
                }
            } 
            iterate();
            iterationsDone_++;
        }
        return iterationsDone_;
    };

    // Function to specify the psf
    virtual void setPSF(
        Array<T,N>& psf
    ) { 
	Estimate<T,N>::setPSF(psf); 
	iterationsDone_ = 0;
        this->est_ = 1;
    };

    // Function to specify the raw image
    virtual void setImage(
        Array<T,N>& img
    ) { 
	Estimate<T,N>::setImage(img); 
        old_.resize(this->est_.extent()); 
	iterationsDone_ = 0;
        this->est_ = 1;
    };

    // Function to set the initial condition of estimate
    void setInitialContitions(
        Array<T,N>& initial
    ) { this->est_ = initial; };

    // Function to set the number of iterations to run
    void setIterations( 
	int iterations 
    ) { iterations_ = iterations; };

    // Function to get the number of iterations to be run
    int getIterations( ) 
    {   return iterations_; };

    // Function to get the number of iterations that have been run
    int getIterationsDone( ) 
    {   return iterationsDone_; };

    // Function to abort em estimation
    void abort() 
    {   running_ = false; };
    
    // Function to resume em estimation
    void resume()
    {   running_ = true; run(); };

    Array<T,N>& getOldEstimate()
    {   return old_; };

  protected:

    // Function for single iteration of the em algorithm
    virtual void iterate() = 0;

  private:

    // not allowed
    EstimateIterative( EstimateIterative<T,N>& );
    void operator=( EstimateIterative<T,N>& );
 
  protected:
    
    int iterations_;		// number of iterations to perform
    int iterationsDone_;	// number of iterations already done
    EstimatePenalty<T,N>* penalty_; // estimation penalty function
    bool running_;		// continue running flag
    Array<T,N> old_; 		// old estimate for error calculations
    std::vector< EstimateObserver<T,N>* > observers_; // observer objects
};

}

#endif  // _ESTIMATE_ITERATIVE_H
