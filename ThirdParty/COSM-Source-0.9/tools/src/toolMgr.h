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

#ifndef TOOL_MGR_H
#define TOOL_MGR_H

#include "wu/wuHeader.h"
#include "blitz/arrayManip.h"
#include "blitz/fftwInterface.h"
#include <blitz/array.h>

#include <string>

namespace cosm {

class wuHeader;

class ToolMgr {

public:
    ToolMgr( 
        const std::string& inName, 
        const std::string& outName = "" 
    );
    ToolMgr( 
        const std::string& inName, 
        const std::string& refName, 
        const std::string& outName 
    );
    ~ToolMgr();

    bool info( wuHeader& header );

    bool addHeader( 
        int xDim, 
        int yDim, 
        int zDim, 
        wuHeader::type type, 
        bool endian 
    );

    bool removeHeader( int );

    bool convert( wuHeader::type type, int );

    bool transformation( 
        int xDim, 
        int yDim, 
        int zDim, 
        int xShift, 
        int yShift, 
        int zShift, 
        double value  
    );

    bool resample( 
        int xDim, 
        int yDim, 
        int zDim, 
        bool up 
    );

    bool compare( 
        double& maximum, 
        double& mean, 
        double& meanSquare 
    ); 

    bool convolve(
    		unsigned short size = 0,
    		bool centeredPsf = true
    		);

    bool shift(
        int xShift,
        int yShift,
        int zShift,
        bool circular
    );
 
    bool scale(
        bool sumFlag = false,
        bool maxFlag = false,
        double value = 1.0
    );

    bool create( 
        int xDim, 
        int yDim, 
        int zDim, 
        wuHeader::type type,
        double value 
    );

    bool ellipsoid( 
        int xCenter, 
        int yCenter, 
        int zCenter, 
        int xRadius, 
        int yRadius, 
        int zRadius, 
        double value  
    );

    bool box( 
        int xCenter, 
        int yCenter, 
        int zCenter, 
        int xRadius, 
        int yRadius, 
        int zRadius, 
        double value  
    );

    bool point( 
        int xCenter, 
        int yCenter, 
        int zCenter, 
        double value  
    );

    bool variant(
        int number,
        int start,
        int size,
    	bool centeredPsf = true
    );

protected:

    template<typename T, int  N>
    Array<T,N> variantImage(
        int numberOfStrata,
        int startOfStrata,
        int sizeOfStrata,
        bool centeredPsf
    );
    
private:
     
    ToolMgr( const ToolMgr& );                    // not allowed
    const ToolMgr& operator()( const ToolMgr& );  // not allowed

private:

    std::string inName_;
    std::string outName_;
    std::string refName_;

};


template<typename T, int  N>
Array<T,N> ToolMgr::variantImage(
    int numberOfStrata,
    int startOfStrata,
    int sizeOfStrata,
    bool centeredPsf
) {
    // read the object image
    Array<T, N> img;
    wuDataRead(img, inName_);

    // find the psf prefix and suffix
    std::string psfname = refName_;
    std::string psfsuffix;
    std::string psfprefix;
    size_t i = psfname.rfind('.', psfname.length());
    if ( i != std::string::npos )
    {
        psfprefix = psfname.substr(0,i);
        psfsuffix = psfname.substr(i, psfname.length());
    }

    Array<RectDomain<N>, 1> strata(numberOfStrata+2);
    Array<Array<T,N>, 1> psfs(numberOfStrata+1);

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
            std::ostringstream psfm;
            psfm << psfprefix << m << psfsuffix;
            std::cout <<psfm.str()<<std::endl;
            wuDataRead(psfs(m), psfm.str());
        }
    }

    // compute the interpolation constants
    Array<T,1> a;
    a.resize(img.extent(0));
    a = 0;
    int index = strata(0).lbound(int(0));
    for ( int m = 0; m < strata.length(0); m++ ) 
    {
        int l = strata(m).lbound(int(0));
        int u = strata(m).ubound(int(0));
        int size = u - l + 1;
        for ( int j = 0; j < size; j++ ) 
        { 
            a(index) = T(1) - T(j)/T(size);
            index++;
        }
    }
    std::cout <<"a: "<< a << std::endl;

    TinyVector<int,N> extent(img.extent());
    Array<T,N> s;              
    Array<T,N> s1;             
    s.resize(extent);
    s1.resize(extent);
    extent(N-1) = extent(N-1)/2+1;

    // resize the psfs and compute the Fourier transform (OTF)
    Array<Array<std::complex<T>,N>, 1> psfsF;  
    psfsF.resize(psfs.length(0));
    for ( int m = 0; m < psfs.length(0); m++ ) 
    {
        psfsF(m).resize(extent);
        Array<T,N> psfResized(img.extent()); 
		//Version 0.8.11::External-3
        //Process Psf file based on the choice of user that whether it's centered
        if (centeredPsf){
        	padAround(psfs(m), psfResized);
        }else{
        	padCenter(psfs(m), psfResized);
        }
        psfsF(m) = forwardFFT(psfResized);
    }

    // compute the image
    fftwInterface<T,N> fftw;
    Array<std::complex<T>,N> imgF;  
    Array<std::complex<T>,N> sF;  

    imgF.resize(extent);
    sF.resize(extent);

    int m;
    // get convolution of img with interpolation of psfs 
    // (multiplication and in Fourier domain)
    imgF = 0;
    int size = strata.length(0) - 2;
    for ( m = 0; m < size; m++ ) 
    {
        s = 0;
        if ( m == 0 )
        {
            s(strata(m)) = img(strata(m));
        }
        else if ( m == size - 1 )
        {
            s(strata(m+2)) = img(strata(m+2));
        }
        s(strata(m+1)) = img(strata(m+1));

        s1 = s;
        multiplyStratum(strata(m+1), s, a, true);
        multiplyStratum(strata(m+1), s1, a, false);

        fftw.plan(s, sF);
        fftw.execute();
        imgF += sF * psfsF(m);

        fftw.plan(s1, sF);
        fftw.execute();
        imgF += sF * psfsF(m+1);
    }
    // convert back to time domain
    fftw.plan(imgF, img);
    fftw.execute();
    img /= img.size();

	//Version 0.8.11::External-3
    //Process Psf file based on the choice of user that whether it's centered
    if (centeredPsf){
		TinyVector<int,3> shift((img.shape()+1)/2);
		Array<T,3> res = circularShift(img, shift);

		return res;
    }else{
        return img;

    }
}

};

#endif // TOOL_MGR_H
