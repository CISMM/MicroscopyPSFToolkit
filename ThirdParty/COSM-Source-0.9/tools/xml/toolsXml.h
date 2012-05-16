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

//Version 0.8.11::External-5
//Support whole XML mechanisim of Tools.

#ifndef TOOLS_XML_H
#define TOOLS_XML_H

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <blitz/array.h>
#include <map>

class TiXmlElement;

namespace cosm{

template<typename T>
class ToolsXml {

public:
    ToolsXml();
    ~ToolsXml();

    bool saveConvolve( const std::string& file );

    bool saveVariant( const std::string& file );

    void imgFile( const std::string& val ) { imgFile_ = val; };
    const std::string& imgFile( void ) const { return imgFile_; };

    void psfFile( const std::string& val ) { psfFile_ = val; };
    const std::string& psfFile( void ) const { return psfFile_; };

    void outFile( const std::string& val ) { outFile_ = val; };
    const std::string& outFile( void ) const { return outFile_; };

    void imgFile2( const std::string& val ) { imgFile2_ = val; };
    const std::string& imgFile2( void ) const { return imgFile2_; };

    void psfFile2( const std::string& val ) { psfFile2_ = val; };
    const std::string& psfFile2( void ) const { return psfFile2_; };

    void outFile2( const std::string& val ) { outFile2_ = val; };
    const std::string& outFile2( void ) const { return outFile2_; };

    void outImgSize( const std::string& val ) { outImgSize_ = val; };
    const std::string& outImgSize( void ) const { return outImgSize_; };

    void variantStrataValue( unsigned short val ) { variantStrataValue_ = val; };
    unsigned short variantStrataValue( void ) const { return variantStrataValue_; };

    void variantStrataStart( unsigned short val ) { variantStrataStart_ = val; };
    unsigned short variantStrataStart( void ) const { return variantStrataStart_; };

    void variantStrataSize( unsigned short val ) { variantStrataSize_ = val; };
    unsigned short variantStrataSize( void ) const { return variantStrataSize_; };

    void convolveCenteredPSF( bool centeredPsf ) { convolveCenteredPsf_ = centeredPsf; };
    bool convolveCenteredPSF( void ) const { return convolveCenteredPsf_; };

    void variantCenteredPSF( bool centeredPsf ) { variantCenteredPsf_ = centeredPsf; };
    bool variantCenteredPSF( void ) const { return variantCenteredPsf_; };

protected:
    std::string imgFile_;
    std::string psfFile_;
    std::string outFile_;
    std::string imgFile2_;
    std::string psfFile2_;
    std::string outFile2_;
    
    std::string outImgSize_;

    unsigned short variantStrataValue_;
    unsigned short variantStrataStart_;
    unsigned short variantStrataSize_;

    bool convolveCenteredPsf_;
    bool variantCenteredPsf_;

};
};

#include "xml/toolsXml.cxx"

#endif // TOOLS_XML_H

