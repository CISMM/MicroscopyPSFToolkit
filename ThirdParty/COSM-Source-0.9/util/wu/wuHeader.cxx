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

#include "wuHeader.h"
#include <fstream>

using namespace cosm;

wuHeader::wuHeader(
) : nx(0), 
    ny(0), 
    nz(0), 
    mode(FLOAT),
    amin(0), 
    amax(0), 
    amean(0), 
    user15(COSM_WU_TVAL), 
    user16_footer_size(0), 
    num_labels(0) 
{
}

wuHeader::wuHeader(
    std::string filename
) {
    FILE *fp;
    if ( (fp=fopen(filename.data(),"rb")) == (FILE*)NULL ) {
        std::cout <<"ERROR: Cant open "<<filename<<" for read"<<std::endl;
        return;
    }
    this->read(fp);
    fclose(fp);
    return;
}

wuHeader::~wuHeader() {};

wuHeader::wuHeader( 
    const wuHeader& header 
) {
    nx = header.nx;
    ny = header.ny;
    nz = header.nz;
    mode = header.mode;
    amin = header.amin;
    amax = header.amax;
    amean = header.amean;
    user1_flags = swap4ByteInt(user1_flags);
    user2 = header.user2;
    user3 = header.user3;
    user4 = header.user4;
    user5 = header.user5;
    user6 = header.user6;
    user7 = header.user7;
    user8 = header.user8;
    user9 = header.user9;
    user10 = header.user10;
    user11 = header.user11;
    user12 = header.user12;
    user13 = header.user13;
    user14 = header.user14;
    user15 = header.user15;
    user16_footer_size = header.user16_footer_size;
    num_labels = header.num_labels;
}

const wuHeader& 
wuHeader::operator=( 
    const wuHeader& header 
) {
    nx = header.nx;
    ny = header.ny;
    nz = header.nz;
    mode = header.mode;
    amin = header.amin;
    amax = header.amax;
    amean = header.amean;
    user1_flags = swap4ByteInt(user1_flags);
    user2 = header.user2;
    user3 = header.user3;
    user4 = header.user4;
    user5 = header.user5;
    user6 = header.user6;
    user7 = header.user7;
    user8 = header.user8;
    user9 = header.user9;
    user10 = header.user10;
    user11 = header.user11;
    user12 = header.user12;
    user13 = header.user13;
    user14 = header.user14;
    user15 = header.user15;
    user16_footer_size = header.user16_footer_size;
    num_labels = header.num_labels;
    return *this;
}

void wuHeader::swapHeader() {
    nx = swap4ByteInt(nx);
    ny = swap4ByteInt(ny);
    nz = swap4ByteInt(nz);
    mode = swap4ByteInt(mode);
    int tmp;
    tmp = swap4ByteInt(*(int*)&amin);
    amin = *(float*)&tmp;
    tmp = swap4ByteInt(*(int*)&amax);
    amax = *(float*)&tmp;
    tmp = swap4ByteInt(*(int*)&amean);
    amean = *(float*)&tmp;
    user1_flags = swap4ByteInt(user1_flags);
    user2 = swap4ByteInt(user2);
    user3 = swap4ByteInt(user3);
    user4 = swap4ByteInt(user4);
    user5 = swap4ByteInt(user5);
    user6 = swap4ByteInt(user6);
    user7 = swap4ByteInt(user7);
    user8 = swap4ByteInt(user8);
    user9 = swap4ByteInt(user9);
    user10 = swap4ByteInt(user10);
    user11 = swap4ByteInt(user11);
    user12 = swap4ByteInt(user12);
    user13 = swap4ByteInt(user13);
    user14 = swap4ByteInt(user14);
    user15 = swap4ByteInt(user15);
    user16_footer_size = swap4ByteInt(user16_footer_size);
    num_labels = swap4ByteInt(num_labels);
}

int wuHeader::read(
    std::string filename
){
    FILE *fp;
    if ( (fp=fopen(filename.data(),"rb")) == (FILE*)NULL ) {
        std::cout <<"open file "<< filename << " failed"<<std::endl;
        return -1;
    }
    int ret = read(fp);
    fclose(fp);
    return ret;
}

int wuHeader::read(
    FILE* pFile
) {
    int bytesRead = 0;
    if ( (bytesRead = fread((void*)this,1,1024,pFile)) != 1024 ) {
        std::cout <<"ERROR: Incomplete read of header, bytes read "<<bytesRead<<std::endl;
        return -1;
    }
    if ( user15 != COSM_WU_TVAL ) {
        if( user15 == COSM_WU_SWAPPED_TVAL ) {
	    std::cout <<"swap header"<<std::endl;
            swapHeader();
	    return 1;
	}
        else {
	    return -1;
        }
    } 
    return 0;
}

int wuHeader::write(
    std::string filename
) {
    FILE *fp;
    if ( (fp=fopen(filename.data(),"wb")) == (FILE*)NULL ) {
        std::cout <<"ERROR: Can't open "<<filename<<" for write"<<std::endl;
        return -1;
    }
    this->write(fp);
    fclose(fp);
    return 0;
}
                                                                                
int wuHeader::write(
    FILE* pFile
) {
    if ( user15 != COSM_WU_TVAL ) {
        std::cout <<"ERROR: Invalid header"<<std::endl;
        return -1;
    }
                                                                                
    if ( nx < 0 || nx > COSM_WU_MAXDIM || 
        ny < 0 || ny > COSM_WU_MAXDIM || 
        nz < 0 || nz > COSM_WU_MAXDIM ) {
        std::cout <<"ERROR: Dimensions in header invalid"<<std::endl;
        std::cout <<"       Image dimensions must be in range [0.."<<COSM_WU_MAXDIM<<"]"<<std::endl;
        return -1;
    }
                                                                                
    if ( fwrite((void*)this,1,1024,pFile) != 1024 ) {
        std::cout <<"ERROR: Incomplete write of header"<<std::endl;
        return -1;
    }
    return 0;
}

int wuHeader::write(
    std::ofstream* file
) {
    if ( user15 != COSM_WU_TVAL ) {
        std::cout <<"ERROR: Invalid header"<<std::endl;
        return -1;
    }
                                                                                
    if ( nx < 0 || nx > COSM_WU_MAXDIM || 
        ny < 0 || ny > COSM_WU_MAXDIM || 
        nz < 0 || nz > COSM_WU_MAXDIM ) {
        std::cout <<"ERROR: Dimensions in header invalid"<<std::endl;
        std::cout <<"       Image dimensions must be in range [0.."<<COSM_WU_MAXDIM<<"]"<<std::endl;
        return -1;
    }
                                                                                
    if ( ! file->write((char*)this, 1024) ) {
        std::cout <<"ERROR: Incomplete write of header"<<std::endl;
        return -1;
    }
    return 0;
}
