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

#include "tinyxml/tinyxml.h"
#include <iostream>
#define STR(s) #s
#define XSTR(s) STR(s)
#define VERSION XSTR(COSM_VERSION)

using namespace cosm;

template<typename T>
PsfXml<T>::PsfXml() {
    stringToModelMap_["Gibson_Lanni1992"] = Psf<T>::MODEL_GIBSON_LANI;
    stringToModelMap_["Haeberle"] = Psf<T>::MODEL_HAEBERLE;
    modelToStringMap_[Psf<T>::MODEL_GIBSON_LANI] = "Gibson_Lanni1992"; 
    modelToStringMap_[Psf<T>::MODEL_HAEBERLE] = "Haeberle";

    stringToEvalMap_["Interpolation"] = Psf<T>::EVAL_INTERPOLATION;
    stringToEvalMap_["Exact"] = Psf<T>::EVAL_EXACT;
    evalToStringMap_[Psf<T>::EVAL_INTERPOLATION] = "Interpolation"; 
    evalToStringMap_[Psf<T>::EVAL_EXACT] = "Exact";

    stringToTypeMap_["Non-Confocal"] = cosm::OPTICAL_SECTIONING_WIDEFIELD;
    stringToTypeMap_["2-Photon"] = cosm::OPTICAL_SECTIONING_2_PHOTON;
    stringToTypeMap_["Circular-Confocal"] = CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE;
    stringToTypeMap_["Line-Confocal"] = cosm::CONFOCAL_ROTATING_DISK_LINE_APERTURE;
	stringToTypeMap_["DIC"] = cosm::DIC;
    typeToStringMap_[cosm::OPTICAL_SECTIONING_WIDEFIELD] = "Non-Confocal";
    typeToStringMap_[cosm::OPTICAL_SECTIONING_2_PHOTON] = "2-Photon";
    typeToStringMap_[CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE] = "Circular-Confocal";
    typeToStringMap_[cosm::CONFOCAL_ROTATING_DISK_LINE_APERTURE] = "Line-Confocal";
	typeToStringMap_[cosm::DIC] = "DIC";

    stringToPrecisionMap_["Single"] = PRECISION_SINGLE;
    stringToPrecisionMap_["Double"] = PRECISION_DOUBLE;
    stringToPrecisionMap_["LongDouble"] = PRECISION_LONG_DOUBLE;
    precisionToStringMap_[PRECISION_SINGLE] = "Single"; 
    precisionToStringMap_[PRECISION_DOUBLE] = "Double"; 
    precisionToStringMap_[PRECISION_LONG_DOUBLE] = "LongDouble"; 

};

template<typename T>
PsfXml<T>::~PsfXml() {};

template<typename T>
T PsfXml<T>::unitConversion( 
    T value, 
    const std::string& unit 
) {
    if      ( unit == "m"  ) value *= 1000; 
    else if ( unit == "cm" ) value *= 10; 
    else if ( unit == "mm" ) value *= 1.0; 
    else if ( unit == "um" ) value *= 1e-3; 
    else if ( unit == "nm" ) value *= 1e-6; 
    return value;
}

template<typename T>
T PsfXml<T>::extractValue( 
    TiXmlElement* element, 
    const std::string& attribute, 
    bool unitFlag 
) {
    T value = 0.0;
    std::stringstream elemStr;
    if ( !element) 
    {
        std::cout <<"Element NULL in XML file"<< std::endl;
        return -1;
    }
    elemStr << *element->Attribute(attribute);
    elemStr >> value;

    if ( unitFlag ) 
    {
        std::string unitStr = element->Attribute("unit");
        value = unitConversion( value, unitStr );
    }
    return value;
}

template<typename T>
bool PsfXml<T>::open( 
    const std::string& filename
) {
    TiXmlDocument doc( filename );
    bool loadOK = doc.LoadFile();
    if ( !loadOK ) 
    {
        if ( doc.Error() ) 
        {
            std::cout << "Error in " << doc.Value() <<": "<<doc.ErrorDesc() << std::endl;
        }
        std::cout << "XML file failed to load" << std::endl;
        return false;
    }
    TiXmlElement* psfElement = doc.FirstChildElement("PSF");
    if ( !psfElement )
    {
        std::cout <<"No PSF Element in XML file"<< std::endl;
        return false;
    }
    std::string model = psfElement->Attribute("model");
    model_ = stringToModelMap_[model];

    std::string eval = psfElement->Attribute("eval");
    eval_ = stringToEvalMap_[eval];

    std::string type = psfElement->Attribute("type");
    type_ = stringToTypeMap_[type];

    std::string precision = psfElement->Attribute("precision");
    precision_ = stringToPrecisionMap_[precision];

    file_ = psfElement->Attribute("file");
    version_ = psfElement->Attribute("version");

    std::stringstream dimension;
    dimension << psfElement->Attribute("dimension");
    dimension >> this->Nxy_ >> this->Nxy_ >> this->Nz_;

    std::stringstream spacing;
    std::string unit = psfElement->Attribute("unit");
    spacing <<  psfElement->Attribute("spacing");
    spacing >> this->deltaXY_ >> this->deltaXY_ >> this->deltaZ_;
    this->deltaXY_ = unitConversion( this->deltaXY_, unit );
    this->deltaZ_ = unitConversion( this->deltaZ_, unit );

    this->na_ = extractValue(psfElement->FirstChildElement("NumericalAperture"), "value");
    this->lambda_ = extractValue(psfElement->FirstChildElement("Wavelength"), "value", true);
    this->lm_ = extractValue(psfElement->FirstChildElement("LateralMagnification"), "value");
    this->absError_ = extractValue(psfElement->FirstChildElement("AbsoluteError"), "value");

    TiXmlElement* specimen =  psfElement->FirstChildElement("Specimen");;
    if ( !specimen ) 
    {
        std::cout <<"No Specimen Element in XML file"<< std::endl;
        return false; 
    }
    this->ts_ = extractValue(specimen->FirstChildElement("Thickness"), "value", true);
    this->ns_ = extractValue(specimen->FirstChildElement("RefractiveIndex"), "value");

    TiXmlElement* immersion = psfElement->FirstChildElement("ImmersionMedium");
    if ( !immersion ) 
    {
        std::cout <<"No Immersion Element in XML file"<< std::endl;
        return false;
    }
    this->tid_ = extractValue(immersion->FirstChildElement("Thickness"), "Design", true);
    this->nid_ = extractValue(immersion->FirstChildElement("RefractiveIndex"), "Design");
    this->tia_ = extractValue(immersion->FirstChildElement("Thickness"), "Actual", true);
    this->nia_ = extractValue(immersion->FirstChildElement("RefractiveIndex"), "Actual");

    TiXmlElement* coverslip = psfElement->FirstChildElement("CoverSlip"); 
    if ( !coverslip ) 
    {
        std::cout <<"No Coverslip Element in XML file"<< std::endl;
        return false;
    }
    this->tgd_ = extractValue(coverslip->FirstChildElement("Thickness"), "Design", true);
    this->ngd_ = extractValue(coverslip->FirstChildElement("RefractiveIndex"), "Design");
    this->tga_ = extractValue(coverslip->FirstChildElement("Thickness"), "Actual", true);
    this->nga_ = extractValue(coverslip->FirstChildElement("RefractiveIndex"), "Actual");

    this->tld_ = extractValue(psfElement->FirstChildElement("OTL"), "Design", true);
    this->tla_ = extractValue(psfElement->FirstChildElement("OTL"), "Actual", true);
	
    TiXmlElement* confocal = psfElement->FirstChildElement("Confocal"); 
    if ( confocal ) 
    {
        this->fsize_ = extractValue(confocal->FirstChildElement("Size"), "value");
		this->distance_ = extractValue(confocal->FirstChildElement("Distance"), "value");
		this->magY_ = extractValue(confocal->FirstChildElement("Magnification"), "value");
    }
    TiXmlElement* dic = psfElement->FirstChildElement("DIC"); 
    if ( dic ) 
    {
        this->shear_ = extractValue(dic->FirstChildElement("Shear"), "value");
		this->bias_ = extractValue(dic->FirstChildElement("Bias"), "value");
		this->amplitudeRatio_ = extractValue(dic->FirstChildElement("AmplitudeRatio"), "value");
    }
    return true;
}

template<typename T>
bool PsfXml<T>::save(
    const std::string& filename
){
    TiXmlDocument doc;
    TiXmlComment* comment;
    TiXmlElement* element;

    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );
                                    
    comment = new TiXmlComment();
	std::string v = " Cosm Version: ";
	v += VERSION;
    comment->SetValue(v);
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue(" Cosm Psf specification file ");
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     model = <\"Gibson_Lanni1992\" | \"Haeberle\" > " );
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     eval = <\"Interpolation\" | \"Exact\"> " );
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     type = <\"Non-confocal\" | \"Circular-confocal\" | \"Line-confocal\"> " );
    doc.LinkEndChild( comment );

    TiXmlElement* psfElement = new TiXmlElement( "PSF" );

    psfElement->SetAttribute("model", modelToStringMap_[model_]);
    psfElement->SetAttribute("eval", evalToStringMap_[eval_]);
    psfElement->SetAttribute("type", typeToStringMap_[type_]);
    psfElement->SetAttribute("precision", precisionToStringMap_[precision_]);
    psfElement->SetAttribute("file", file_);
    std::stringstream dimension;
    dimension << this->Nxy_ <<" "<< this->Nxy_ <<" "<< this->Nz_;
    psfElement->SetAttribute("dimension", dimension.str());
    std::stringstream spacing;
    spacing << this->deltaXY_ <<" "<< this->deltaXY_ <<" "<< this->deltaZ_;
    psfElement->SetAttribute("spacing", spacing.str());
    psfElement->SetAttribute("unit", "mm");
    psfElement->SetAttribute("version", version_);

    element = new TiXmlElement( "AbsoluteError" );
    element->SetDoubleAttribute("value", this->absError_);
    psfElement->LinkEndChild( element);

    element = new TiXmlElement( "LateralMagnification");
    element->SetDoubleAttribute("value", this->lm_);
    psfElement->LinkEndChild( element);

    element = new TiXmlElement( "NumericalAperture");
    element->SetDoubleAttribute("value", this->na_);
    psfElement->LinkEndChild( element);

    element = new TiXmlElement( "Wavelength");
    element->SetDoubleAttribute("value", this->lambda_);
    element->SetAttribute("unit", "mm");
    psfElement->LinkEndChild( element);

    TiXmlElement* coverSlip = new TiXmlElement( "CoverSlip");
    element = new TiXmlElement( "RefractiveIndex");
    element->SetDoubleAttribute("Design", this->ngd_);
    element->SetDoubleAttribute("Actual", this->nga_);
    coverSlip->LinkEndChild( element);
    element = new TiXmlElement( "Thickness");
    element->SetDoubleAttribute("Design", this->tgd_);
    element->SetDoubleAttribute("Actual", this->tga_);
    element->SetAttribute("unit", "mm");
    coverSlip->LinkEndChild( element);
    psfElement->LinkEndChild( coverSlip);

    TiXmlElement* immersionMedium = new TiXmlElement( "ImmersionMedium");
    element = new TiXmlElement( "RefractiveIndex");
    element->SetDoubleAttribute("Design", this->nid_);
    element->SetDoubleAttribute("Actual", this->nia_);
    immersionMedium->LinkEndChild( element);
    element = new TiXmlElement( "Thickness");
    element->SetDoubleAttribute("Design", this->tid_);
    element->SetDoubleAttribute("Actual", this->tia_);
    element->SetAttribute("unit", "mm");
    immersionMedium->LinkEndChild( element);
    psfElement->LinkEndChild( immersionMedium);

    TiXmlElement* specimen = new TiXmlElement( "Specimen");
    element = new TiXmlElement( "RefractiveIndex");
    element->SetDoubleAttribute("value", this->ns_);
    specimen->LinkEndChild( element);
    element = new TiXmlElement( "Thickness");
    element->SetDoubleAttribute("value", this->ts_);
    element->SetAttribute("unit", "mm");
    specimen->LinkEndChild( element);
    psfElement->LinkEndChild( specimen);

    element = new TiXmlElement( "OTL");
    element->SetDoubleAttribute("Design", this->tld_);
    element->SetDoubleAttribute("Actual", this->tla_);
    element->SetAttribute("unit", "mm");
    psfElement->LinkEndChild( element);
	
	if ( type_ == CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE || 
         type_ == CONFOCAL_ROTATING_DISK_LINE_APERTURE )
    {
        TiXmlElement* confocal = new TiXmlElement("Confocal");
		element = new TiXmlElement("Size");
        element->SetDoubleAttribute("value", this->fsize_);
		confocal->LinkEndChild( element );
		element = new TiXmlElement("Distance");
        element->SetDoubleAttribute("value", this->distance_);
		confocal->LinkEndChild( element );
		element = new TiXmlElement("Magnification");
        element->SetDoubleAttribute("value", this->magY_);
		confocal->LinkEndChild( element );
		psfElement->LinkEndChild( confocal);
    }
	if ( type_ == DIC )
    {
        TiXmlElement* dic = new TiXmlElement("DIC");
		element = new TiXmlElement("Shear");
        element->SetDoubleAttribute("value", this->shear_);
		dic->LinkEndChild( element );
		element = new TiXmlElement("Bias");
        element->SetDoubleAttribute("value", this->bias_);
		dic->LinkEndChild( element );
		element = new TiXmlElement("AmplitudeRatio");
        element->SetDoubleAttribute("value", this->amplitudeRatio_);
		dic->LinkEndChild( element );
		psfElement->LinkEndChild( dic );
    }

    doc.LinkEndChild( psfElement );

    return doc.SaveFile( filename );
}
