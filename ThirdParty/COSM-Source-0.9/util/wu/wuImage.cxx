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

#include "wu/wuImage.h"
 
using namespace blitz;
using namespace cosm; 

void 
WUImage::header( 
    const wuHeader& val 
) { 
    header_ = val; 
}

void WUImage::Create(
    int x,
    int y,
    int z,
    unsigned char value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::BYTE); 
    ushortArray_.resize(z,y,x);
    ushortArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    short value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::SHORT); 
    ushortArray_.resize(z,y,x);
    ushortArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    unsigned short value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::USHORT); 
    ushortArray_.resize(z,y,x);
    ushortArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    int value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::INT); 
    ushortArray_.resize(z,y,x);
    ushortArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    float value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::FLOAT); 
    floatArray_.resize(z,y,x);
    floatArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    double value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::DOUBLE); 
    doubleArray_.resize(z,y,x);
    doubleArray_ = value;
}

void WUImage::Create(
    int x,
    int y,
    int z,
    long double value 
){
    header_.columns(x);
    header_.rows(y);
    header_.sections(z);
    header_.dataType(cosm::wuHeader::LONG_DOUBLE); 
    longDoubleArray_.resize(z,y,x);
    longDoubleArray_ = value;
}

bool WUImage::ReadData( 
    const std::string& file,
    bool headerFlag 
) {
    int res = -1;
    if ( headerFlag )
    {
        header_.read(file);
    }
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
	    if ( headerFlag )
            {
                res = wuDataRead(byteArray_, file);
            } else {
                res = wuDataRead(byteArray_, file, header_);
            }
            break;
        case cosm::wuHeader::SHORT: 
	    if ( headerFlag )
            {
                res = wuDataRead(shortArray_, file);
            } else {
                res = wuDataRead(shortArray_, file, header_);
            }
            break;
        case cosm::wuHeader::USHORT: 
	    if ( headerFlag )
            {
                res = wuDataRead(ushortArray_, file);
            } else {
                res = wuDataRead(ushortArray_, file, header_);
            }
            break;
        case cosm::wuHeader::INT: 
	    if ( headerFlag )
            {
                res = wuDataRead(intArray_, file);
            } else {
                res = wuDataRead(intArray_, file, header_);
            }
            break;
        case cosm::wuHeader::FLOAT: 
	    if ( headerFlag )
            {
                res = wuDataRead(floatArray_, file);
            } else {
                res = wuDataRead(floatArray_, file, header_);
            }
            break;
        case cosm::wuHeader::DOUBLE: 
	    if ( headerFlag )
            {
                res = wuDataRead(doubleArray_, file);
            } else {
                res = wuDataRead(doubleArray_, file, header_);
            }
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
	    if ( headerFlag )
            {
                res = wuDataRead(longDoubleArray_, file);
            } else {
                res = wuDataRead(longDoubleArray_, file, header_);
            }
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res == -1 ? false : true;
}

bool WUImage::WriteData( 
    const std::string& file,
    bool headerFlag 
) {
    int res = -1;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res = wuDataWrite(byteArray_, file, headerFlag);
            break;
        case cosm::wuHeader::SHORT: 
            res = wuDataWrite(shortArray_, file, headerFlag);
            break;
        case cosm::wuHeader::USHORT: 
            res = wuDataWrite(ushortArray_, file, headerFlag);
            break;
        case cosm::wuHeader::INT: 
            res = wuDataWrite(intArray_, file, headerFlag);
            break;
        case cosm::wuHeader::FLOAT: 
            res = wuDataWrite(floatArray_, file, headerFlag);
            break;
        case cosm::wuHeader::DOUBLE: 
            res = wuDataWrite(doubleArray_, file, headerFlag);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res = wuDataWrite(longDoubleArray_, file, headerFlag);
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res == -1 ? false : true;
}

Array<unsigned char,3> 
WUImage::GetByteArray( void )  
{
    Array<unsigned char,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = byteArray_;
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<unsigned char>(byteArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<unsigned char>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<unsigned char>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<unsigned char>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<unsigned char>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res =  cast<unsigned char>(longDoubleArray_);
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}

Array<short,3> 
WUImage::GetShortArray( void )  
{
    Array<short,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<short>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = shortArray_;
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<short>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<short>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<short>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<short>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res =  cast<short>(longDoubleArray_);
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}

Array<unsigned short,3> 
WUImage::GetUshortArray( void )  
{
    Array<unsigned short,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<unsigned short>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<unsigned short>(shortArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = ushortArray_;
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<unsigned short>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<unsigned short>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<unsigned short>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res =  cast<unsigned short>(longDoubleArray_);
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}

Array<int,3> 
WUImage::GetIntArray( void )  
{
    Array<int,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<int>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<int>(shortArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<int>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = intArray_;
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<int>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<int>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res =  cast<int>(longDoubleArray_);
            break;
       default:
           std::cout << "unknown data type" << std::endl;
    }
    return res;
}

Array<float,3> 
WUImage::GetFloatArray( void )  
{
    Array<float,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<float>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<float>(shortArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<float>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<float>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = floatArray_;
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<float>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res =  cast<float>(longDoubleArray_);
            break;
       default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}

Array<double,3> 
WUImage::GetDoubleArray( void )
{
    Array<double,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<double>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<double>(shortArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<double>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<double>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<double>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = doubleArray_;
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res = cast<double>(longDoubleArray_);
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}
   
Array<long double,3> 
WUImage::GetLongDoubleArray( void )  
{
    Array<long double,3> res;
    switch ( header_.dataType() ) 
    {
        case cosm::wuHeader::BYTE: 
            res.resize(byteArray_.extent());
            res = cast<long double>(byteArray_);
            break;
        case cosm::wuHeader::SHORT: 
            res.resize(shortArray_.extent());
            res = cast<long double>(shortArray_);
            break;
        case cosm::wuHeader::USHORT: 
            res.resize(ushortArray_.extent());
            res = cast<long double>(ushortArray_);
            break;
        case cosm::wuHeader::INT: 
            res.resize(intArray_.extent());
            res = cast<long double>(intArray_);
            break;
        case cosm::wuHeader::FLOAT: 
            res.resize(floatArray_.extent());
            res = cast<long double>(floatArray_);
            break;
        case cosm::wuHeader::DOUBLE: 
            res.resize(doubleArray_.extent());
            res = cast<long double>(doubleArray_);
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            res.resize(longDoubleArray_.extent());
            res = longDoubleArray_;
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    }
    return res;
}

void 
WUImage::ConvertToByte( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            break;
        case cosm::wuHeader::SHORT: 
            byteArray_.resize(shortArray_.extent());
            byteArray_ = cast<unsigned char>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            byteArray_.resize(ushortArray_.extent());
            byteArray_ = cast<unsigned char>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            byteArray_.resize(intArray_.extent());
            byteArray_ = cast<unsigned char>(intArray_);
            intArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            byteArray_.resize(floatArray_.extent());
            byteArray_ = cast<unsigned char>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            byteArray_.resize(doubleArray_.extent());
            byteArray_ = cast<unsigned char>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            byteArray_.resize(longDoubleArray_.extent());
            byteArray_ = cast<unsigned char>(longDoubleArray_);
            longDoubleArray_.free();
            break;
       default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::BYTE); 
}

void 
WUImage::ConvertToShort( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            shortArray_.resize(byteArray_.extent());
            shortArray_ = cast<short>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            break;
        case cosm::wuHeader::USHORT: 
            shortArray_.resize(ushortArray_.extent());
            shortArray_ = cast<short>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            shortArray_.resize(intArray_.extent());
            shortArray_ = cast<short>(intArray_);
            intArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            shortArray_.resize(floatArray_.extent());
            shortArray_ = cast<short>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            shortArray_.resize(doubleArray_.extent());
            shortArray_ = cast<short>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            shortArray_.resize(longDoubleArray_.extent());
            shortArray_ = cast<short>(longDoubleArray_);
            longDoubleArray_.free();
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::SHORT); 
}

void 
WUImage::ConvertToUshort( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            ushortArray_.resize(byteArray_.extent());
            ushortArray_ = cast<unsigned short>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            ushortArray_.resize(shortArray_.extent());
            ushortArray_ = cast<unsigned short>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            break;
        case cosm::wuHeader::INT: 
            ushortArray_.resize(intArray_.extent());
            ushortArray_ = cast<unsigned short>(intArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            ushortArray_.resize(floatArray_.extent());
            ushortArray_ = cast<unsigned short>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            ushortArray_.resize(doubleArray_.extent());
            ushortArray_ = cast<unsigned short>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            ushortArray_.resize(longDoubleArray_.extent());
            ushortArray_ = cast<unsigned short>(longDoubleArray_);
            longDoubleArray_.free();
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::USHORT); 
}

void 
WUImage::ConvertToInt( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            intArray_.resize(byteArray_.extent());
            intArray_ = cast<int>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            intArray_.resize(shortArray_.extent());
            intArray_ = cast<int>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            intArray_.resize(ushortArray_.extent());
            intArray_ = cast<int>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            break;
        case cosm::wuHeader::FLOAT: 
            intArray_.resize(floatArray_.extent());
            intArray_ = cast<int>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            intArray_.resize(doubleArray_.extent());
            intArray_ = cast<int>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            intArray_.resize(longDoubleArray_.extent());
            intArray_ = cast<int>(longDoubleArray_);
            longDoubleArray_.free();
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::INT); 
}

void 
WUImage::ConvertToFloat( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            floatArray_.resize(byteArray_.extent());
            floatArray_ = cast<float>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            floatArray_.resize(shortArray_.extent());
            floatArray_ = cast<float>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            floatArray_.resize(ushortArray_.extent());
            floatArray_ = cast<float>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            floatArray_.resize(intArray_.extent());
            floatArray_ = cast<float>(intArray_);
            intArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            break;
        case cosm::wuHeader::DOUBLE: 
            floatArray_.resize(doubleArray_.extent());
            floatArray_ = cast<float>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            floatArray_.resize(longDoubleArray_.extent());
            floatArray_ = cast<float>(longDoubleArray_);
            longDoubleArray_.free();
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::FLOAT); 
}
 
void
WUImage::ConvertToDouble( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            doubleArray_.resize(byteArray_.extent());
            doubleArray_ = cast<double>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            doubleArray_.resize(shortArray_.extent());
            doubleArray_ = cast<double>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            doubleArray_.resize(ushortArray_.extent());
            doubleArray_ = cast<double>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            doubleArray_.resize(intArray_.extent());
            doubleArray_ = cast<double>(intArray_);
            intArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            doubleArray_.resize(doubleArray_.extent());
            doubleArray_ = cast<double>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            doubleArray_.resize(longDoubleArray_.extent());
            doubleArray_ = cast<double>(longDoubleArray_);
            longDoubleArray_.free();
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::DOUBLE); 
}
 
void
WUImage::ConvertToLongDouble( void )
{
    switch ( header_.dataType() )
    {
        case cosm::wuHeader::BYTE: 
            longDoubleArray_.resize(byteArray_.extent());
            longDoubleArray_ = cast<long double>(byteArray_);
            byteArray_.free();
            break;
        case cosm::wuHeader::SHORT: 
            longDoubleArray_.resize(shortArray_.extent());
            longDoubleArray_ = cast<long double>(shortArray_);
            shortArray_.free();
            break;
        case cosm::wuHeader::USHORT: 
            longDoubleArray_.resize(ushortArray_.extent());
            longDoubleArray_ = cast<long double>(ushortArray_);
            ushortArray_.free();
            break;
        case cosm::wuHeader::INT: 
            longDoubleArray_.resize(intArray_.extent());
            longDoubleArray_ = cast<long double>(intArray_);
            intArray_.free();
            break;
        case cosm::wuHeader::FLOAT: 
            floatArray_.resize(floatArray_.extent());
            longDoubleArray_ = cast<long double>(floatArray_);
            floatArray_.free();
            break;
        case cosm::wuHeader::DOUBLE: 
            longDoubleArray_.resize(doubleArray_.extent());
            longDoubleArray_ = cast<long double>(doubleArray_);
            doubleArray_.free();
            break;
        case cosm::wuHeader::LONG_DOUBLE: 
            break;
        default:
            std::cout << "unknown data type" << std::endl;
    } 
    header_.dataType(cosm::wuHeader::LONG_DOUBLE); 
} 

unsigned char* 
WUImage::GetByteData(
    void
) 
{
    return byteArray_.data();
}

short* 
WUImage::GetShortData(
    void
) 
{
    return shortArray_.data();
}

unsigned short* 
WUImage::GetUshortData(
    void
) 
{
    return ushortArray_.data();
}

int*
WUImage::GetIntData(
    void
) 
{
    return intArray_.data();
}

float* 
WUImage::GetFloatData(
    void
) 
{
    return floatArray_.data();
}

double* 
WUImage::GetDoubleData(
    void
){ 
    return doubleArray_.data();
}

long double* 
WUImage::GetLongDoubleData(
    void
){ 
    return longDoubleArray_.data();
}
