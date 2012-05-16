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

#ifndef PSF_XML_H
#define PSF_XML_H

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include "psf/psfParameters.h"
#include "psf/psf.h"
#include <map>

class TiXmlElement;

namespace cosm{

template<typename T>
class PsfXml : public PsfParameters<T> {
public:
    enum Precision {
        PRECISION_SINGLE = 0x0,
        PRECISION_DOUBLE = 0x1,
        PRECISION_LONG_DOUBLE = 0x2
    };

public:
    PsfXml();
    ~PsfXml();

    bool open( const std::string& file );
    bool save( const std::string& file );

    void model( typename Psf<T>::Model val ) { model_ = val; };
    void model( const std::string& val ) { model_ = stringToModelMap_[val]; };
    typename Psf<T>::Model model( void ) const { return model_; };
    const std::string& modelStr( void ) { return modelToStringMap_[model_]; };

    void eval( typename Psf<T>::Eval val ) { eval_ = val; };
    void eval( const std::string& val ) { eval_ = stringToEvalMap_[val]; };
    typename Psf<T>::Eval eval( void ) const { return eval_; };
    const std::string& evalStr( void ) { return evalToStringMap_[eval_]; };

    void type( PsfType val ) { type_ = val; };
    void type( const std::string& val ) { type_ = stringToTypeMap_[val]; };
    PsfType type( void ) const { return type_; };
    const std::string& typeStr( void ) { return typeToStringMap_[type_]; };

    void precision( Precision val ) { precision_ = val; };
    void precision( const std::string& val ) { precision_ = stringToPrecisionMap_[val]; };
    Precision precision( void ) { return precision_; };
    const std::string& precisionStr( void ) { return precisionToStringMap_[precision_]; };

    void file( const std::string& val ) { file_ = val; };
    const std::string& file( void ) const { return file_; };

    void version( const std::string& val ) { version_ = val; };
    const std::string& version( void ) const { return version_; };

    void sum( bool sum ) { sum_ = sum; };
    bool sum( void ) const { return sum_; };
	
	void centered( bool center ) { centered_ = center; };
	bool centered( void ) const { return centered_; };
	
	void realAndImag( bool val ) {realAndImag_ = val; };
	bool realAndImag( void ) const { return realAndImag_; };

protected:
    T unitConversion( 
        T value, 
        const std::string& unit 
    ); 

    T extractValue(  
        TiXmlElement* element, 
        const std::string& attribute, 
        bool unitFlag = false 
    );

protected:
    // Not allowed
    PsfXml( const PsfXml& );
    // Not allowed
    const PsfXml& operator=( const PsfXml& );

protected:
    std::string file_;
    std::string version_;    
    typename Psf<T>::Model model_;
    typename Psf<T>::Eval eval_;
    PsfType type_;
    Precision precision_;
    bool sum_;
	bool centered_;
	bool realAndImag_;
    
    std::map<std::string, typename Psf<T>::Model> stringToModelMap_;
    std::map<typename Psf<T>::Model, std::string> modelToStringMap_;
    std::map<std::string, typename Psf<T>::Eval> stringToEvalMap_;
    std::map<typename Psf<T>::Eval, std::string> evalToStringMap_;
    std::map<std::string, PsfType> stringToTypeMap_;
    std::map< PsfType, std::string> typeToStringMap_;
    std::map<std::string, Precision> stringToPrecisionMap_;
    std::map< Precision, std::string> precisionToStringMap_;
};
};

#include "xml/psfXml.cxx"

#endif // PSF_XML_H
