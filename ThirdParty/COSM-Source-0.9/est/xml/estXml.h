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

#ifndef EST_XML_H
#define EST_XML_H

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include "est/estimate.h"
#include <map>

class TiXmlElement;

namespace cosm{

template<typename T>
class EstXml {
public:
    enum Algorithm {
        ALGORITHM_LLS  = 0x0,
        ALGORITHM_MAP  = 0x1,
        ALGORITHM_EM   = 0x2,
        ALGORITHM_JVC  = 0x3,
        ALGORITHM_EMSV = 0x4,
        ALGORITHM_EMOS = 0x5
    };

    enum PenaltyType {
        PENALTY_NONE = 0x0,
        PENALTY_INTENSITY = 0x1,
        PENALTY_ROUGHNESS = 0x2
    };

    enum OTFType {
        OTF_TYPE_MEMORY = 0x1,
        OTF_TYPE_DISK   = 0x2
    };

public:
    EstXml();
    ~EstXml();

    bool open( const std::string& file );
    bool save( const std::string& file );

    void algorithm( Algorithm val ) { algorithm_ = val; };
    void algorithm( const std::string& val ) { algorithm_ = stringToAlgorithmMap_[val]; };
    Algorithm algorithm( void ) const { return algorithm_; };
    const std::string& algorithmStr( void ) const { return algorithmToStringMap_[algorithm_]; };

    void imgFile( const std::string& val ) { imgFile_ = val; };
    const std::string& imgFile( void ) const { return imgFile_; };

    void psfFile( const std::string& val ) { psfFile_ = val; };
    const std::string& psfFile( void ) const { return psfFile_; };

    void estFile( const std::string& val ) { estFile_ = val; };
    const std::string& estFile( void ) const { return estFile_; };

    void version( const std::string& val ) { version_ = val; };
    const std::string& version( void ) const { return version_; };

    void centeredPSF( bool center ) { centeredPSF_ = center; };
    bool centeredPSF( void ) const { return centeredPSF_; };

    void linearParameter( double val ) { linearParameter_ = val; };
    double linearParameter( void ) const { return linearParameter_; };

    void iterationNumber( unsigned long val ) { iterationNumber_ = val; };
    unsigned long iterationNumber( void ) const { return iterationNumber_; };

    void iterationUpdate( unsigned long val ) { iterationUpdate_ = val; };
    unsigned long iterationUpdate( void ) const { return iterationUpdate_; };

    void iterationWrite( unsigned long val ) { iterationWrite_ = val; };
    unsigned long iterationWrite( void ) const { return iterationWrite_; };

    void invariantDoubleZ( bool val ) { invariantDoubleZ_ = val; };
    bool invariantDoubleZ( void ) const { return invariantDoubleZ_; };

    void invariantPenaltyType( PenaltyType val ) { invariantPenaltyType_ = val; };
    void invariantPenaltyType( const std::string& val ) { invariantPenaltyType_ = stringToPenaltyTypeMap_[val]; };
    PenaltyType invariantPenaltyType( void ) const { return invariantPenaltyType_; };
    const std::string& invariantPenaltyTypeStr( void ) const { return penaltyTypeToStringMap_[invariantPenaltyType_]; };

    void invariantPenaltyValue( double val ) { invariantPenaltyValue_ = val; };
    double invariantPenaltyValue( void ) const { return invariantPenaltyValue_; };

    void variantDoubleZ( bool val ) { variantDoubleZ_ = val; };
    bool variantDoubleZ( void ) const { return variantDoubleZ_; };

    void variantPenaltyType( PenaltyType val ) { variantPenaltyType_ = val; };
    void variantPenaltyType( const std::string& val ) { variantPenaltyType_ = stringToPenaltyTypeMap_[val]; };
    PenaltyType variantPenaltyType( void ) const { return variantPenaltyType_; };
    const std::string& variantPenaltyTypeStr( void ) const { return penaltyTypeToStringMap_[variantPenaltyType_]; };

    void variantPenaltyValue( double val ) { variantPenaltyValue_ = val; };
    double variantPenaltyValue( void ) const { return variantPenaltyValue_; };

    void variantOtfType( OTFType val ) { variantOtfType_ = val; };
    void variantOtfType( const std::string& val ) { variantOtfType_ = stringToOtfTypeMap_[val]; };
    OTFType variantOtfType( void ) const { return variantOtfType_; };
    const std::string& variantOtfTypeStr( void ) const { return otfTypeToStringMap_[variantOtfType_]; };

    void variantOtfName( const std::string& val ) { variantOtfName_ = val; };
    const std::string& variantOtfName( void ) const { return variantOtfName_; };

    void variantStrataValue( unsigned short val ) { variantStrataValue_ = val; };
    unsigned short variantStrataValue( void ) const { return variantStrataValue_; };

    void variantStrataStart( unsigned short val ) { variantStrataStart_ = val; };
    unsigned short variantStrataStart( void ) const { return variantStrataStart_; };

    void variantStrataSize( unsigned short val ) { variantStrataSize_ = val; };
    unsigned short variantStrataSize( void ) const { return variantStrataSize_; };

    void convergenceType( int val ) { convergenceType_ = val; };
    int convergenceType( void ) const { return convergenceType_; };

    void convergenceName( const std::string& val ) { convergenceName_ = val; };
    const std::string& convergenceName( void ) const { return convergenceName_; };

protected:
    // Not allowed
    EstXml( const EstXml& );
    // Not allowed
    const EstXml& operator=( const EstXml& );

    T extractValue(
        TiXmlElement* element,
        const std::string& attribute
    );

protected:
    std::string imgFile_;
    std::string psfFile_;
    std::string estFile_;
    std::string version_;    
    bool centeredPSF_;
    
    std::map<std::string, Algorithm> stringToAlgorithmMap_;
    std::map<std::string, PenaltyType> stringToPenaltyTypeMap_;
    std::map<std::string, OTFType> stringToOtfTypeMap_;
    std::map<std::string, ErrorType> stringToErrorTypeMap_;

    std::map< Algorithm, std::string> algorithmToStringMap_;
    std::map< PenaltyType, std::string> penaltyTypeToStringMap_;
    std::map< OTFType, std::string> otfTypeToStringMap_;
    std::map< ErrorType, std::string> errorTypeToStringMap_;

    Algorithm algorithm_;
    double linearParameter_;
    unsigned long iterationNumber_;
    unsigned long iterationUpdate_;
    unsigned long iterationWrite_;
    bool invariantDoubleZ_;
    PenaltyType invariantPenaltyType_;
    double invariantPenaltyValue_;
    bool variantDoubleZ_;
    PenaltyType variantPenaltyType_;
    double variantPenaltyValue_;
    unsigned short variantStrataValue_;
    unsigned short variantStrataStart_;
    unsigned short variantStrataSize_;
    OTFType variantOtfType_;
    std::string variantOtfName_;
    int convergenceType_;
    std::string convergenceName_;
    
};
};

#include "xml/estXml.cxx"

#endif // EST_XML_H

