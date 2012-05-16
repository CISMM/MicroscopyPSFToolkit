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
 
#ifndef _WU_IMAGE_H
#define _WU_IMAGE_H
 
#include "wu/wuHeader.h"

namespace cosm {

class WUImage {

public:
    WUImage() {};
    ~WUImage() {};
 
    void Create( int x, int y, int z, unsigned char value );
    void Create( int x, int y, int z, short value );
    void Create( int x, int y, int z, unsigned short value );
    void Create( int x, int y, int z, int value );
    void Create( int x, int y, int z, float value );
    void Create( int x, int y, int z, double value );
    void Create( int x, int y, int z, long double value );

    bool ReadData( const std::string& file, bool headerFlag = true );
    bool WriteData( const std::string& file, bool headerFlag = true );

    void ConvertToByte( void );
    void ConvertToShort( void );
    void ConvertToUshort( void );
    void ConvertToInt( void );
    void ConvertToFloat( void );
    void ConvertToDouble( void );
    void ConvertToLongDouble( void );

    bool IsByte() { return header_.dataType() == wuHeader::BYTE; };   
    bool IsShort() { return header_.dataType() == wuHeader::SHORT; };   
    bool IsUshort() { return header_.dataType() == wuHeader::USHORT; };   
    bool IsInt() { return header_.dataType() == wuHeader::INT; };   
    bool IsFloat() { return header_.dataType() == wuHeader::FLOAT; };   
    bool IsDouble() { return header_.dataType() == wuHeader::DOUBLE; };   
    bool IsLongDouble() { return header_.dataType() == wuHeader::LONG_DOUBLE; };   
    wuHeader::type GetType( ) { return header_.dataType(); };
    void header( const wuHeader& val ); 
    const wuHeader& header( ) { return header_; };

    Array<unsigned char,3> GetByteArray( void );  
    Array<short,3> GetShortArray( void );  
    Array<unsigned short,3> GetUshortArray( void );  
    Array<int,3> GetIntArray( void );  
    Array<float,3> GetFloatArray( void );  
    Array<double,3> GetDoubleArray( void );  
    Array<long double,3> GetLongDoubleArray( void );  

    unsigned char* GetByteData(void); 
    short* GetShortData(void); 
    unsigned short* GetUshortData(void); 
    int* GetIntData(void); 
    float* GetFloatData(void); 
    double* GetDoubleData(void); 
    long double* GetLongDoubleData(void); 
    
protected:
    WUImage( const WUImage& ); // not allowed
    const WUImage& operator=( const WUImage& ); // not allowed

protected:
    wuHeader header_;
    Array<unsigned char,3> byteArray_;
    Array<short,3> shortArray_;
    Array<unsigned short,3> ushortArray_;
    Array<int,3> intArray_;
    Array<float,3> floatArray_;
    Array<double,3> doubleArray_;
    Array<long double,3> longDoubleArray_;

};

};

#endif // WU_IMAGE_H
