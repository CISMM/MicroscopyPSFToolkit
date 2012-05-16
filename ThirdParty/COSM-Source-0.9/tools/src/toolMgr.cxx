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

#include "toolMgr.h"
#include "wu/wuHeader.h"
#include "wu/wuImage.h"
#include "est/estimate.h"
#include "blitz/arrayManip.h"
#include "blitz/convolution.h"
#include "blitz/fftwInterface.h"
#include <blitz/array.h>

using namespace cosm;
using namespace blitz;

ToolMgr::ToolMgr( 
    const std::string& inName, 
    const std::string& outName 
) : inName_(inName), outName_(outName), refName_("") {
}

ToolMgr::ToolMgr( 
    const std::string& inName, 
    const std::string& refName, 
    const std::string& outName 
) : inName_(inName), outName_(outName), refName_(refName) {
}

ToolMgr::~ToolMgr() {}

bool 
ToolMgr::info( 
    wuHeader& header 
) {
    if ( !inName_.empty() )
    {
        header.read(inName_);
        return true;
    }
    return false;
}

bool 
ToolMgr::addHeader( 
    int xDim, 
    int yDim, 
    int zDim, 
    wuHeader::type type, 
    bool endian 
) {
    wuHeader header;
    header.columns(xDim);
    header.rows(yDim);
    header.sections(zDim);
    header.dataType(type);
    header.swap(endian);

    WUImage image;
    image.header(header);
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }
    if ( !image.ReadData(inName_, false) )
    {
        return false;
    }
    if ( !image.WriteData(outName_) )
    {
        return false;
    }
    return true;
}

bool 
ToolMgr::removeHeader( 
	int numberOfFile
) {
    WUImage image;
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }

    if ( numberOfFile <= 1){
		if ( image.ReadData(inName_) && image.WriteData(outName_, false) )
		{
			return true;
		}
    }else{
		
		//Version 0.8.11::Improve-4
		//Support processing of multiple input files

        // find the I/O files' prefix and suffix
        std::string inputfile = inName_;
        std::string inputprefix;
        std::string inputsuffix;

        std::string outfile = outName_;
        std::string outprefix;
        std::string outsuffix;

        size_t i = inputfile.rfind('.', inputfile.length());
        if ( i != std::string::npos )
        {
        	inputprefix = inputfile.substr(0,i);
        	inputsuffix = inputfile.substr(i, inputfile.length());
        }

        size_t j = outfile.rfind('.', outfile.length());
        if ( j != std::string::npos )
        {
        	outprefix = outfile.substr(0,j);
        	outsuffix = outfile.substr(j, outfile.length());
        }

        for ( int m = 0; m < numberOfFile; m++){
            std::ostringstream inputindex;
            std::ostringstream outindex;
            inputindex << inputprefix << m;
            outindex << outprefix << m;

            std::string subInputfile = inputindex.str() + inputsuffix;
            std::string subOutfile = outindex.str() + outsuffix;
    		if ( !image.ReadData(subInputfile) )
    		{
    			return false;
    		}
    		if ( !image.WriteData(subOutfile, false) )
    		{
    			return false;
    		}
    	}
        return true;
    }
    return false;
}

bool 
ToolMgr::convert( 
    wuHeader::type type,
    int numberOfFile
) {
    WUImage image;

	if ( inName_.empty() || outName_.empty() )
	{
		return false;
	}

	if ( numberOfFile <= 1){
		if ( !image.ReadData(inName_) )
		{
			return false;
		}
		switch (type )
		{
			case wuHeader::FLOAT:
			{
				blitz::Array<float,3> res = image.GetFloatArray();
				return wuDataWrite(res, outName_) == 0;
			}
			case wuHeader::DOUBLE:
			{
				blitz::Array<double,3> res = image.GetDoubleArray();
				return wuDataWrite(res, outName_) == 0;
			}
			case wuHeader::LONG_DOUBLE:
			{
				blitz::Array<long double,3> res = image.GetLongDoubleArray();
				return wuDataWrite(res, outName_) == 0;
			}
			case wuHeader::USHORT:
			{
				blitz::Array<unsigned short,3> res = image.GetUshortArray();
				return wuDataWrite(res, outName_) == 0;
			}
			default:
				return false;
		}

    } else{
	    //Version 0.8.11::Improve-4
		//Support processing of multiple input files

        // find the I/O files' prefix and suffix
        std::string inputfile = inName_;
        std::string inputprefix;
        std::string inputsuffix;

        std::string outfile = outName_;
        std::string outprefix;
        std::string outsuffix;

        size_t i = inputfile.rfind('.', inputfile.length());
        if ( i != std::string::npos )
        {
        	inputprefix = inputfile.substr(0,i);
        	inputsuffix = inputfile.substr(i, inputfile.length());
        }

        size_t j = outfile.rfind('.', outfile.length());
        if ( j != std::string::npos )
        {
        	outprefix = outfile.substr(0,j);
        	outsuffix = outfile.substr(j, outfile.length());
        }

        for ( int m = 0; m < numberOfFile; m++){
            std::ostringstream inputindex;
            std::ostringstream outindex;
            inputindex << inputprefix << m;
            outindex << outprefix << m;

            std::string subInputfile = inputindex.str() + inputsuffix;
            std::string subOutfile = outindex.str() + outsuffix;
    		if ( !image.ReadData(subInputfile) )
    		{
    			return false;
    		}
    		switch (type )
    		{
    			case wuHeader::FLOAT:
    			{
    				blitz::Array<float,3> res = image.GetFloatArray();
    				if(!(wuDataWrite(res, subOutfile) == 0)){
    					return false;
    				}
    			}
    			break;

    			case wuHeader::DOUBLE:
    			{
    				blitz::Array<double,3> res = image.GetDoubleArray();
    				if(!(wuDataWrite(res, subOutfile) == 0)){
    					return false;
    				}
    			}
    			break;

    			case wuHeader::LONG_DOUBLE:
    			{
    				blitz::Array<long double,3> res = image.GetLongDoubleArray();
    				if(!(wuDataWrite(res, subOutfile) == 0)){
    					return false;
    				}
    			}
    			break;

    			case wuHeader::USHORT:
    			{
    				blitz::Array<unsigned short,3> res = image.GetUshortArray();
    				if(!(wuDataWrite(res, subOutfile) == 0)){
    					return false;
    				}
    			}
    			break;

    			default:
    				return false;
    		}
        }
        return true;
    }
	return false;
} 

bool
ToolMgr::transformation( 
    int xDim, 
    int yDim, 
    int zDim, 
    int xShift, 
    int yShift, 
    int zShift, 
    double value  
){
    WUImage data;
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }

    if ( !data.ReadData(inName_) )
    {
        return false;
    }
    blitz::TinyVector<int,3> shift(zShift,yShift,xShift);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(zDim,yDim,zDim);
        transform(image, res, shift, (float)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(zDim,yDim,zDim);
        transform(image, res, shift, value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(zDim,yDim,xDim);
        transform(image, res, shift, (long double)value);
        return wuDataWrite(res, outName_) == 0;
    }

    else if ( data.IsUshort() )
    {
        blitz::Array<unsigned short, 3> image = data.GetUshortArray();
        blitz::Array<unsigned short, 3> res(zDim,yDim,xDim);
        transform(image, res, shift, (unsigned short)value);
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool 
ToolMgr::resample( 
    int xDim, 
    int yDim, 
    int zDim, 
    bool up 
){
    WUImage data;
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }
    if ( !data.ReadData(inName_) )
    {
        return false;
    }
    blitz::TinyVector<int,3> factor(zDim,yDim,xDim);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res;
        if ( !up ) downsample(image, res, factor);
        if (  up ) upsample(image, res, factor);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res;
        if ( !up ) downsample(image, res, factor);
        if (  up ) upsample(image, res, factor);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res;
        if ( !up ) downsample(image, res, factor);
        if (  up ) upsample(image, res, factor);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsUshort() )
    {
        blitz::Array<unsigned short, 3> image = data.GetUshortArray();
        blitz::Array<unsigned short, 3> res;
        if ( !up ) downsample(image, res, factor);
        if (  up ) upsample(image, res, factor);
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool 
ToolMgr::compare( 
    double& maximum,
    double& mean,
    double& meanSquare
){
    if ( inName_.empty() || refName_.empty() )
    {
        return false;
    }

    // create the reference
    WUImage ref;
    if ( !ref.ReadData(refName_) )
    {
        return false;
    }

    // create the data 
    WUImage data;
    if ( !data.ReadData(inName_) )
    {
        return false;
    }

    if ( data.IsFloat() && ref.IsFloat() )
    {
        ErrorEstimate<float,3> estimate;
        blitz::Array<float,3> dataImg = data.GetFloatArray();
        blitz::Array<float,3> refImg = ref.GetFloatArray();
        maximum = (double)estimate.error(MAXIMUM_ERROR, dataImg, refImg);
        mean = (double)estimate.error(MEAN_ERROR, dataImg, refImg);
        meanSquare = (double)estimate.error(MEAN_SQUARE_ERROR, dataImg, refImg);
        return true;
    }
    else if ( data.IsDouble() && ref.IsDouble() )
    {
        ErrorEstimate<double,3> estimate;
        blitz::Array<double,3> dataImg = data.GetDoubleArray();
        blitz::Array<double,3> refImg = ref.GetDoubleArray();
        maximum = estimate.error(MAXIMUM_ERROR, dataImg, refImg);
        mean = estimate.error(MEAN_ERROR, dataImg, refImg);
        meanSquare = estimate.error(MEAN_SQUARE_ERROR, dataImg, refImg);
        return true;
    }
    else if ( data.IsLongDouble() & ref.IsLongDouble() )
    {
        ErrorEstimate<long double,3> estimate;
        blitz::Array<long double,3> dataImg = data.GetLongDoubleArray();
        blitz::Array<long double,3> refImg = ref.GetLongDoubleArray();
        maximum = (double)estimate.error(MAXIMUM_ERROR, dataImg, refImg);
        mean = (double)estimate.error(MEAN_ERROR, dataImg, refImg);
        meanSquare = (double)estimate.error(MEAN_SQUARE_ERROR, dataImg, refImg);
        return true;
    }
    return false;
} 

bool 
ToolMgr::convolve( 
    unsigned short size,
	bool centeredPsf
){
    if ( inName_.empty() || outName_.empty() || refName_.empty() )
    {
        return false;
    }
    WUImage img1;
    WUImage img2;
    if ( !img1.ReadData(inName_) ) 
    {
        return false;
    }
    if ( !img2.ReadData(refName_) )
    {
        return false;
    }

    if ( img1.IsFloat() && img2.IsFloat() )
    {
        blitz::Array<float, 3> res = convolutionWithFFT(img1.GetFloatArray(),img2.GetFloatArray(), size, centeredPsf);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( img1.IsDouble() && img2.IsDouble() )
    {
        blitz::Array<double, 3> res = convolutionWithFFT(img1.GetDoubleArray(),img2.GetDoubleArray(), size, centeredPsf);
        return wuDataWrite(res, outName_) == 0;


    }
    else if ( img1.IsLongDouble() && img2.IsLongDouble() )
    {
        blitz::Array<long double, 3> res = convolutionWithFFT(img1.GetLongDoubleArray(),img2.GetLongDoubleArray(), size, centeredPsf);
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool 
ToolMgr::create( 
    int xDim, 
    int yDim, 
    int zDim, 
    wuHeader::type type,
    double value  
){
    WUImage data;
    if ( outName_.empty() )
    {
        return false;
    }
    switch( type )
    {
        case wuHeader::FLOAT: data.Create(xDim,yDim,zDim, (float)value); break;
        case wuHeader::DOUBLE: data.Create(xDim,yDim,zDim, (double)value); break;
        case wuHeader::LONG_DOUBLE: data.Create(xDim,yDim,zDim, (long double)value); break;
        case wuHeader::USHORT: data.Create(xDim,yDim,zDim, (unsigned short)value); break;
        default:
            return false;
    }
    if ( !data.WriteData(outName_) )
    {
        return false;
    }
    return true;
}

bool 
ToolMgr::ellipsoid( 
    int xCenter, 
    int yCenter, 
    int zCenter, 
    int xRadius, 
    int yRadius, 
    int zRadius, 
    double value  
){
    WUImage data;
    if ( inName_.empty() &&  outName_.empty() )
    {
        return false; 
    }
    if ( !data.ReadData(inName_) )
    {
        return false; 
    }

    blitz::TinyVector<int,3> center(zCenter,yCenter,xCenter);
    blitz::TinyVector<int,3> radius(zRadius, yRadius, xRadius);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(image.extent());
        addEllipsoid(image, res, center, radius, (float)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(image.extent());
        addEllipsoid(image, res, center, radius, (double)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(image.extent());
        addEllipsoid(image, res, center, radius, (long double)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsUshort() )
    {
        blitz::Array<unsigned short, 3> image = data.GetUshortArray();
        blitz::Array<unsigned short, 3> res(image.extent());
        addEllipsoid(image, res, center, radius, (unsigned short)value);
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool 
ToolMgr::box( 
    int xCenter, 
    int yCenter, 
    int zCenter, 
    int xRadius, 
    int yRadius, 
    int zRadius, 
    double value  
){
    WUImage data;
    if ( inName_.empty() &&  outName_.empty() )
    {
        return false; 
    }
    if ( !data.ReadData(inName_) )
    {
        return false; 
    }

    blitz::TinyVector<int,3> center(zCenter,yCenter,xCenter);
    blitz::TinyVector<int,3> radius(zRadius, yRadius, xRadius);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(image.extent());
        addBox(image, res, center, radius, (float)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(image.extent());
        addBox(image, res, center, radius, (double)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(image.extent());
        addBox(image, res, center, radius, (long double)value);
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsUshort() )
    {
        blitz::Array<unsigned short, 3> image = data.GetUshortArray();
        blitz::Array<unsigned short, 3> res(image.extent());
        addBox(image, res, center, radius, (unsigned short)value);
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool 
ToolMgr::point( 
   int xCenter, 
   int yCenter, 
   int zCenter, 
   double value  
){

    WUImage data;
    if ( inName_.empty() &&  outName_.empty() )
    {
        return false; 
    }
    if ( !data.ReadData(inName_) )
    {
        return false; 
    }
    blitz::TinyVector<int,3> center(zCenter,yCenter,xCenter);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(image.extent());
        res = image;
        res(center) = (float)value;
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(image.extent());
        res = image;
        res(center) = (double)value;
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(image.extent());
        res = image;
        res(center) = (long double)value;
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsUshort() )
    {
        blitz::Array<unsigned short, 3> image = data.GetUshortArray();
        blitz::Array<unsigned short, 3> res(image.extent());
        res = image;
        res(center) = (short)value;
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}

bool
ToolMgr::shift(
    int xShift,
    int yShift,
    int zShift,
    bool circular
){
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }

    WUImage data;
    if ( !data.ReadData(inName_) )
    {
        return false;
    }
    blitz::TinyVector<int,3> shift(zShift,yShift,xShift);
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(image.shape());
        if ( circular )
        {
            res = circularShift(image, shift);
        }
        else
        {
            res = normalShift(image, shift);
        }
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() )
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(image.shape());
        if ( circular )
        {
            res = circularShift(image, shift);
        }
        else
        {
            res = normalShift(image, shift);
        }
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() )
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(image.shape());
        if ( circular )
        {
            res = circularShift(image, shift);
        }
        else
        {
            res = normalShift(image, shift);
        }
        return wuDataWrite(res, outName_) == 0;
    }
    return false;
}


bool 
ToolMgr::scale(
    bool sumFlag,
    bool maxFlag,
    double value
){
    if ( inName_.empty() || outName_.empty() )
    {
        return false;
    }

    WUImage data;
    if ( !data.ReadData(inName_) )
    {
        return false;
    }
    if ( data.IsFloat() ) 
    {
        blitz::Array<float, 3> image = data.GetFloatArray();
        blitz::Array<float, 3> res(image.shape());
        if ( sumFlag ) image /= (sum)(image);
        else if ( maxFlag ) image /= (max)(image);
        else image /= value;
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsDouble() ) 
    {
        blitz::Array<double, 3> image = data.GetDoubleArray();
        blitz::Array<double, 3> res(image.shape());
        if ( sumFlag ) image /= (sum)(image);
        else if ( maxFlag ) image /= (max)(image);
        else image /= value;
        return wuDataWrite(res, outName_) == 0;
    }
    else if ( data.IsLongDouble() ) 
    {
        blitz::Array<long double, 3> image = data.GetLongDoubleArray();
        blitz::Array<long double, 3> res(image.shape());
        if ( sumFlag ) image /= (sum)(image);
        else if ( maxFlag ) image /= (max)(image);
        else image /= value;
        return wuDataWrite(res, outName_) == 0;
    }
    
    return false;
}

bool
ToolMgr::variant(
    int number,
    int start,
    int size,
    bool centeredPsf
) {
    if ( inName_.empty() || outName_.empty() || refName_.empty() )
    {
        return false;
    }
    WUImage object;
    if ( !object.ReadData(inName_) ) 
    {
        return false;
    }
    std::string psfname = refName_;
    std::string psfsuffix;
    std::string psfprefix;
    size_t i = psfname.rfind('.', psfname.length());
    if ( i != std::string::npos )
    {
        psfprefix = psfname.substr(0,i);
        psfsuffix = psfname.substr(i, psfname.length());
    }

    // read in the first image for type check
    WUImage psf;
    if ( !psf.ReadData(psfprefix + "0" + psfsuffix) ) 
    {
        return false;
    }

    if ( object.IsFloat() && psf.IsFloat() )
    {
        blitz::Array<float, 3> res = variantImage<float,3>(number, start, size, centeredPsf);
        wuDataWrite(res, outName_);
        return true;
    }
    else if ( object.IsDouble() && psf.IsDouble() )
    {
        blitz::Array<double, 3> res = variantImage<double,3>(number, start, size, centeredPsf);
        wuDataWrite(res, outName_);
        return true;
    }
    else if ( object.IsLongDouble() && psf.IsLongDouble() )
    {
        blitz::Array<long double, 3> res = variantImage<long double, 3>(number, start, size, centeredPsf);
        wuDataWrite(res, outName_);
        return true;
    }
    return false;
}
