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

#ifndef _WU_HEADER_H
#define _WU_HEADER_H

//#include <blitz/config.h>
#include <blitz/array.h>
#include <iostream>

using namespace blitz;

namespace cosm {

#define swap4ByteInt(x) \
        ( ((x & 0xff000000) >> 24) \
        + ((x & 0x00ff0000) >>  8) \
        + ((x & 0x0000ff00) <<  8) \
        + ((x & 0x000000ff) << 24) )


const unsigned int COSM_WU_MAXDIM = 10000;
const unsigned int COSM_WU_MAX_STRINGS = 10;
const unsigned int COSM_WU_MAX_STRING_LEN = 80;

// Magic number in header to make sure header VALID 
// and byte swapped correctly
const int COSM_WU_TVAL = 1234;
// What TVAL looks like if byte-swapped
const int COSM_WU_SWAPPED_TVAL = -771489792;

// used to manipulate the dirty bit 
const unsigned int COSM_WU_STAT_FLAG_SET   = 0x00000001;
const unsigned int COSM_WU_STAT_FLAG_RESET = 0xFFFFFFFE;

// The header is a structure which contains all the elements found in an
// image file.  The first section of the structure is 1024 bytes.  This 
// corresponds to the MRC HEADER used by UCSF algorithms.  
class  wuHeader {


  public:

    //  Data types for 'mode' field
    enum type {
        BYTE   		= 0,
        SHORT  		= 1,
        FLOAT  		= 2,
	DOUBLE 		= 3,
        USHORT 		= 4,
        INT    		= 5,
	LONG_DOUBLE 	= 6
    };

    wuHeader();
    wuHeader( std::string filename );
    ~wuHeader();

    wuHeader( const wuHeader& header );
    const wuHeader& operator=( const wuHeader& header );

    void addHistory ( char* text );	
    void clearNotes ( );
    void printNotes ( int start, int end );
    void addComment ( char* text );
    void getComment ( );

    int read( FILE* pFile );
    int read( std::string filename );
    int write( FILE* pFILE );
    int write( std::ofstream* file );
    int write( std::string filename );

    unsigned int columns() { return nx; };
    void columns( unsigned int n ) { nx = n; };
    unsigned int rows() { return ny; };
    void rows( unsigned int n ) { ny = n; };
    unsigned int sections() { return nz; };
    void sections( unsigned int n ) { nz = n; };

    float minimum() { return amin; };
    void minimum( float m ) { amin = m; };
    float maximum() { return amax; };
    void maximum( float m ) { amax = m; };
    float mean() { return amean; };
    void mean( float m ) { amean = m; };
    
    void dataType( wuHeader::type type ) { mode = type; };
    wuHeader::type dataType() { return (wuHeader::type) mode; };

    wuHeader::type dataType(char tmp) { return wuHeader::BYTE; };
    wuHeader::type dataType(short tmp) { return wuHeader::SHORT; };
    wuHeader::type dataType(float tmp) { return wuHeader::FLOAT; };
    wuHeader::type dataType(unsigned short tmp) { return wuHeader::USHORT; };
    wuHeader::type dataType(int tmp) { return wuHeader::INT; };
    wuHeader::type dataType(double tmp) { return wuHeader::DOUBLE; };
    wuHeader::type dataType(long double tmp) { return wuHeader::LONG_DOUBLE; };

    std::string dataTypeStr() { 
        switch ( mode ) {
            case wuHeader::BYTE: return "BYTE";
            case wuHeader::SHORT: return "SHORT";
            case wuHeader::FLOAT: return "FLOAT";
            case wuHeader::USHORT: return "USHORT";
            case wuHeader::INT: return "INT";
            case wuHeader::DOUBLE: return "DOUBLE";
            case wuHeader::LONG_DOUBLE: return "LONGDOUBLE";
			default: return "UNKNOWN";
       }
	   return "UNKNOWN";
    };

    int dataSize() {
        switch ( mode ) {
            case wuHeader::BYTE: return nx*ny*nz*sizeof(char);
            case wuHeader::SHORT: return nx*ny*nz*sizeof(short);
            case wuHeader::FLOAT: return nx*ny*nz*sizeof(float);
            case wuHeader::USHORT: return nx*ny*nz*sizeof(unsigned short);
            case wuHeader::INT: return nx*ny*nz*sizeof(int);
            case wuHeader::DOUBLE: return nx*ny*nz*sizeof(double);
            case wuHeader::LONG_DOUBLE: return nx*ny*nz*sizeof(long double);
			default: return 0;
       }
	   return 0;
    };
    void swap( bool val ) { user15 = (val == true ? COSM_WU_SWAPPED_TVAL : COSM_WU_TVAL); };
    int swapping( void ) {
        return (user15 == COSM_WU_TVAL) ? 0 : (user15 == COSM_WU_SWAPPED_TVAL ? 1 : -1);
    };

  private:

    void swapHeader();

  private:

    unsigned int nx;		// number of columns 
    unsigned int ny;		// number of rows 
    unsigned int nz;		// number of sections 
    unsigned int mode;		// data type 
 
    unsigned int xstart; 	// number of first column in map 
    unsigned int ystart; 	// number of first row in map 
    unsigned int zstart; 	// number of first section in map

    unsigned int mx;		// number of intervals along X 
    unsigned int my;		// number of intervals along Y
    unsigned int mz;		// number of intervals along Z

    // cell dimensions (cell/mxyz = pixel spacing)
    float xlength;		
    float ylength;
    float zlength; 

    // cell angles (degrees)
    float alpha;		
    float beta;
    float gamma;

    // which axis corresponds to columns, rows, sections (X=1,Y=2,Z=3)
    unsigned int col_axis;	
    unsigned int row_axis;	
    unsigned int sect_axis;     

    float amin;			// minimum intensity value 
    float amax;			// maximum intensity value
    float amean;		// mean intensity value

    unsigned int ispg;		// space group number 
    unsigned int nsymbt;	// number of bytes used for symmetry operators
 
    int	user1_flags;	 	// used for dataset dirty bit  
    int user2;
    int user3;
    int user4;
    int user5;
    int user6;
    int user7;
    int user8;
    int user9;
    int user10;
    int user11;
    int user12;
    int user13;
    int user14;
    int user15;
    unsigned int user16_footer_size; // size of the footer in bytes

    short id_type;	// type of data set 
    short lens_type;	// lens type 

    short n1;		// unknown 
    short n2;		// unknown
    int data_value2;	// not a readable format 

    // two tilt sets (original 1-3 and current 3-6) 
    float tilt1;	
    float tilt2;	
    float tilt3;
    float tilt4;
    float tilt5;
    float tilt6;

    // multiple wave length information 
    int wave1;		
    int wave2;
    int wave3;

    // x,y,z origin of the image
    float xorigin;	
    float yorigin;
    float zorigin;

    // number of text strings being used
    int num_labels;		   
    // storage for first ten strings
    char label[COSM_WU_MAX_STRINGS][COSM_WU_MAX_STRING_LEN]; 
};


template<typename T>
void swapData( 
    T* data, 
    int size 
) {
     char* tmpData;
     char  tmp;
     int sizeT = sizeof(T)/2;
     for ( int i = 0; i < size; i++ ) {
        tmpData = (char*)&data[i];
        int last = sizeof(T)-1;
        for ( int j = 0; j < sizeT; j++ ) {
            tmp = tmpData[j];
            tmpData[j] = tmpData[last-j];
            tmpData[last-j] = tmp;
        }
    }
}

template<typename T, int N>
int wuDataRead(
    Array<T,N>& A, 
    std::string filename 
) {
    FILE *fp;
    if ( (fp=fopen(filename.data(),"rb")) == (FILE*)NULL ) {
	std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    int bytesRead = 0;
    wuHeader header;
    int swap = header.read(fp);
    T tmp = 0;
    if ( header.dataType( ) != header.dataType(tmp) ) {
	std::cout <<"Incorrect data type"<<std::endl;;
        fclose(fp);
	return -1;
    }

    // the order is oposite the Array order
    if ( N == 3 ) {
        A.resize(header.sections(), header.rows(), header.columns());
    } else if ( N == 2 && header.columns() == 1 ) {
        A.resize(header.sections(), header.rows());
    } else {
	std::cout <<"Incorrect dimensions"<<std::endl;;
        fclose(fp);
	return -1;
    }
    std::cout <<"["<<header.sections()<<","<<header.rows()<<","<<header.columns()<<"]"<<std::endl; 
    int data_size =
        header.columns() * header.rows() * header.sections() * sizeof(T);
    if ( (bytesRead = fread((void*)A.data(),1,data_size,fp)) != data_size ) {
        std::cout <<"WARNING: Incomplete read of image data, bytes read: "<<bytesRead<<std::endl;
        std::cout <<"         Data may be invalid (file="<<filename<<")"<<std::endl;
        fclose(fp);
	return -1;
    }
    // byte-swap image data if needed
    if ( swap == 1 ) {
        std::cout <<"Swapping data"<<std::endl;
        swapData((T*)A.data(), data_size/sizeof(T));
    } else {
        std::cout <<"No data swapping"<<std::endl;
    }
    fclose(fp);
    return 0;
}

template<typename T, int N>
int wuDataRead(
    Array<T,N>& A, 
    std::string filename,
    wuHeader& header  
) {
    FILE *fp;
    if ( (fp=fopen(filename.data(),"rb")) == (FILE*)NULL ) {
	std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    int bytesRead = 0;
    int swap = header.swapping();
    // the order is oposite the Array order
    if ( N == 3 ) {
        A.resize(header.sections(), header.rows(), header.columns());
    } else if ( N == 2 && header.columns() == 1 ) {
        A.resize(header.sections(), header.rows());
    } else {
	std::cout <<"Incorrect dimensions"<<std::endl;;
        fclose(fp);
	return -1;
    }
    std::cout <<"["<<header.sections()<<","<<header.rows()<<","<<header.columns()<<"]"<<std::endl; 
    int data_size =
        header.columns() * header.rows() * header.sections() * sizeof(T);
    if ( (bytesRead = fread((void*)A.data(),1,data_size,fp)) != data_size ) {
        std::cout <<"WARNING: Incomplete read of image data, bytes read: "<<bytesRead<<std::endl;
        std::cout <<"         Data may be invalid (file="<<filename<<")"<<std::endl;
        fclose(fp);
	return -1;
    }
    // byte-swap image data if needed
    if ( swap == 1 ) {
        std::cout <<"Swapping data"<<std::endl;
        swapData((T*)A.data(), data_size/sizeof(T));
    } else {
        std::cout <<"No data swapping"<<std::endl;
    }
    fclose(fp);
    return 0;
}

template<typename T, int N>
int wuDataWrite(
    Array<T,N>& A, 
    std::string filename,
    bool headerFlag = true
) {
    T tmp = 0;
    FILE *fp;
    if ( (fp=fopen(filename.data(),"wb")) == (FILE*)NULL ) {
	std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    wuHeader header;
    // the order is opposite the Array order
    header.columns(1);
    if ( N == 3 ) {
        header.columns(A.extent(2));
    }
    header.rows(A.extent(1));
    header.sections(A.extent(0));
    std::cout <<"columns: "<<header.columns()<<", rows: "<<header.rows() <<", sections: "<<header.sections() << std::endl; 
    header.minimum((float)(blitz::min)(A));
    header.maximum((float)(blitz::max)(A));
    header.mean((float)blitz::mean(A));
    std::cout <<"min: "<<header.minimum()<<", max: "<<header.maximum() <<", mean: "<<header.mean() << std::endl; 
    header.dataType(header.dataType(tmp));
    if ( headerFlag )
    {
        header.write(fp);
    }
    int data_size =
        header.columns() * header.rows() * header.sections() * sizeof(T);
    if ( fwrite((void*)A.data(),1,data_size,fp) != (size_t)data_size ) {
        std::cout <<"WARNING: Writing data for file "<<filename <<" incomplete"<<std::endl;
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

template<typename T, int N>
int wuDataRead(
    Array<complex<T>,N>& A, 
    std::string filename
) {
    FILE *fp;
    if ( (fp=fopen(filename.data(),"rb")) == (FILE*)NULL ) {
	std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    int bytesRead = 0;
    wuHeader header;
    int swap = header.read(fp);
    T tmp = 0;
    if ( header.dataType( ) != header.dataType(tmp) ) {
	std::cout <<"Incorrect data type"<<std::endl;;
	fclose(fp);
	return -1;
    }
    int typeSize = sizeof(complex<T>)/sizeof(T);
    // the order is oposite the Array order
    if ( N == 3 ) {
        A.resize(header.sections(), header.rows(), header.columns()/typeSize);
    } else if ( N == 2 && header.columns() == 1 ) {
        A.resize(header.sections(), header.rows()/typeSize);
    } else {
	std::cout <<"Incorrect dimensions"<<std::endl;;
	fclose(fp);
	return -1;
    }
    int data_size =
        header.columns() * header.rows() * header.sections() * sizeof(T);
    if ( (bytesRead = fread((void*)A.data(),1,data_size,fp)) != data_size ) {
        std::cout <<"WARNING: Incomplete read of image data, bytes read: "<<bytesRead<<std::endl;
        std::cout <<"         Data may be invalid (file="<<filename<<")"<<std::endl;
	fclose(fp);
	return -1;
    }
    // byte-swap image data if needed
    if ( swap == 1 ) {
        std::cout <<"Swapping data not implemented"<<std::endl;
	fclose(fp);
	return -1;
    } else {
        // std::cout <<"No data swapping"<<std::endl;
    }
    fclose(fp);
    return 0;
}

template<typename T, int N>
int wuDataWrite(
    Array<std::complex<T>,N>& A, 
    std::string filename,
    bool headerFlag = true 
) {
    T tmp = 0;
    FILE *fp;
    int writeBytes;
    if ( (fp=fopen(filename.data(),"wb")) == (FILE*)NULL ) {
	std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    wuHeader header;
    // the order is opposite the Array order
    int typeSize = sizeof(complex<T>)/sizeof(T);
    header.columns(1);
    header.sections(A.extent(0));
    if ( N == 3 ) {
        header.columns(A.extent(2)*typeSize);
        header.rows(A.extent(1));
    }
    else if ( N == 2 ) {
        header.rows(A.extent(1)*typeSize);
    } else {
	std::cout <<"Unsupported dimensions" <<std::endl;
    }
    header.minimum((float)(std::min)((blitz::min)(blitz::real(A)), (blitz::min)(blitz::imag(A))));
    header.maximum((float)(std::max)((blitz::max)(blitz::real(A)), (blitz::max)(blitz::imag(A))));
    header.mean((float)(blitz::mean(real(A))+blitz::mean(imag(A)))/2.0);
    header.dataType(header.dataType(tmp));
    if ( headerFlag )
    {
        header.write(fp);
    }
    int data_size =
        header.columns() * header.rows() * header.sections() * sizeof(T);
    if ( (writeBytes = fwrite((void*)A.data(),1,data_size,fp)) != data_size ) {
        std::cout <<"WARNING: Writing data for file "<<filename <<" incompletem writeBytes: "<<writeBytes<<std::endl;
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

}

#endif
