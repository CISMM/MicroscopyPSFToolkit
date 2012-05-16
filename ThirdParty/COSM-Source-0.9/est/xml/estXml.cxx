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
EstXml<T>::EstXml() {
    stringToAlgorithmMap_["LLS"]  = ALGORITHM_LLS;
    stringToAlgorithmMap_["MAP"]  = ALGORITHM_MAP;
    stringToAlgorithmMap_["EM"]   = ALGORITHM_EM;
    stringToAlgorithmMap_["JVC"]  = ALGORITHM_JVC;
    stringToAlgorithmMap_["EMSV"] = ALGORITHM_EMSV;
    stringToAlgorithmMap_["EMOS"] = ALGORITHM_EMOS;
    algorithmToStringMap_[ALGORITHM_LLS]  = "LLS";
    algorithmToStringMap_[ALGORITHM_MAP]  = "MAP";
    algorithmToStringMap_[ALGORITHM_EM]   =  "EM";
    algorithmToStringMap_[ALGORITHM_JVC]  = "JVC";
    algorithmToStringMap_[ALGORITHM_EMSV] = "EMSV";
    algorithmToStringMap_[ALGORITHM_EMOS] = "EMOS";

    stringToPenaltyTypeMap_["none"] = PENALTY_NONE;
    stringToPenaltyTypeMap_["intensity"] = PENALTY_INTENSITY;
    stringToPenaltyTypeMap_["roughness"] = PENALTY_ROUGHNESS;
    penaltyTypeToStringMap_[PENALTY_NONE] = "none";
    penaltyTypeToStringMap_[PENALTY_INTENSITY] = "intensity";
    penaltyTypeToStringMap_[PENALTY_ROUGHNESS] = "roughness";

    stringToOtfTypeMap_["memory"] = OTF_TYPE_MEMORY;
    stringToOtfTypeMap_["disk"] = OTF_TYPE_DISK;
    otfTypeToStringMap_[OTF_TYPE_MEMORY] = "memory";
    otfTypeToStringMap_[OTF_TYPE_DISK] = "disk";

    stringToErrorTypeMap_["maximum"] = MAXIMUM_ERROR;
    stringToErrorTypeMap_["mean"] = MEAN_ERROR;
    stringToErrorTypeMap_["meanSquare"] = MEAN_SQUARE_ERROR;
    stringToErrorTypeMap_["loglikelihood"] = LOG_LIKELYHOOD_ERROR;
    stringToErrorTypeMap_["idivergence"] = I_DIVERGENCE_ERROR;
    errorTypeToStringMap_[MAXIMUM_ERROR] = "maximum";
    errorTypeToStringMap_[MEAN_ERROR] = "mean";
    errorTypeToStringMap_[MEAN_SQUARE_ERROR] = "meanSquare";
    errorTypeToStringMap_[LOG_LIKELYHOOD_ERROR] = "loglikelihood";
    errorTypeToStringMap_[I_DIVERGENCE_ERROR] = "idivergence";
};

template<typename T>
EstXml<T>::~EstXml() {};

template<typename T>
T EstXml<T>::extractValue( 
    TiXmlElement* element, 
    const std::string& attribute
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
    return value;
}

template<typename T>
bool EstXml<T>::open( 
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
    TiXmlElement* estElement = doc.FirstChildElement("Estimation");
    if ( !estElement )
    {
        std::cout <<"No Estimation Element in XML file"<< std::endl;
        return false;
    }
    std::string algorithm = estElement->Attribute("algorithm");
    algorithm_ = stringToAlgorithmMap_[algorithm];

    imgFile_ = estElement->Attribute("imgFile");
    psfFile_ = estElement->Attribute("psfFile");
    estFile_ = estElement->Attribute("estFile");

    std::string centeredPSF = estElement->Attribute("centeredPSF");
    centeredPSF_ = centeredPSF == "1";

    version_ = estElement->Attribute("version");


    TiXmlElement* linearElement = estElement->FirstChildElement("Linear");
    if ( linearElement )
    {
        linearParameter_ = extractValue(linearElement->FirstChildElement("Parameter"), "value");
    }
    TiXmlElement* iterativeElement = estElement->FirstChildElement("Iterative");
    if ( iterativeElement )
    {
        iterationNumber_ = int(extractValue(iterativeElement, "iterations"));
        iterationUpdate_ = int(extractValue(iterativeElement, "update"));
        iterationWrite_ = int(extractValue(iterativeElement, "write"));

        TiXmlElement* invariantElement = iterativeElement->FirstChildElement("Invariant");
        if ( invariantElement )
        {
            invariantDoubleZ_ = extractValue(invariantElement->FirstChildElement("DoubleZ"), "value") > 0.001;
            TiXmlElement* invariantPenaltyElement = invariantElement->FirstChildElement("Penalty");
            if ( invariantPenaltyElement )
            {
                std::string type = invariantPenaltyElement->Attribute("type");
                invariantPenaltyType_ = stringToPenaltyTypeMap_[type];
            }
            invariantPenaltyValue_ = extractValue(invariantElement->FirstChildElement("Penalty"), "value");
        }

        TiXmlElement* variantElement = iterativeElement->FirstChildElement("Variant");
        if ( variantElement )
        {
            variantDoubleZ_ = extractValue(invariantElement->FirstChildElement("DoubleZ"), "value") > 0.001;
            TiXmlElement* variantPenaltyElement = variantElement->FirstChildElement("Penalty");
            if ( variantPenaltyElement )
            {
                std::string type = variantPenaltyElement->Attribute("type");
                variantPenaltyType_ = stringToPenaltyTypeMap_[type];
            }
            variantPenaltyValue_ = extractValue(variantElement->FirstChildElement("Penalty"), "value");

            variantStrataValue_ = extractValue(variantElement->FirstChildElement("Strata"), "value");
            variantStrataStart_ = extractValue(variantElement->FirstChildElement("Strata"), "start");
            variantStrataSize_ = extractValue(variantElement->FirstChildElement("Strata"), "size");
            TiXmlElement* variantOtfElement = variantElement->FirstChildElement("OTF");
            if ( variantOtfElement )
            {
                std::string type = variantOtfElement->Attribute("type");
                variantOtfType_ = stringToOtfTypeMap_[type];
                variantOtfName_ = variantOtfElement->Attribute("name");
            }
        }

        TiXmlElement* convergenceElement = iterativeElement->FirstChildElement("Convergence");
        if ( convergenceElement )
        {
			convergenceType_ = 0;
		    TiXmlElement* convergenceTypeElement = convergenceElement->FirstChildElement("ConvergenceType");
			while( convergenceTypeElement ) 
			{
			    std::string type = convergenceTypeElement->Attribute("value");
				convergenceType_ |= stringToErrorTypeMap_[type];
				convergenceTypeElement = convergenceTypeElement->NextSiblingElement("ConvergenceType");
			}
            convergenceName_ = convergenceElement->Attribute("name");
        }
    }

    return true;
}

template<typename T>
bool EstXml<T>::save(
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
    comment->SetValue(" Cosm Estimation specification file ");
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     algorithm    = < \"LLS\" | \"MAP\" | \"EM\" | \"JVC\" | \"EMSV\" | \"EMOS\" > " );
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     penalty type = < \"none\" | \"intensity\" | \"roughness\" > " );
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     OTF type     = < \"memory\" | \"disk\" > " );
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     convergence type = <\"maximum\" | \"mean\" | \"meanSquare\" | \"loglikelihood\" | \"idivergence\"> " );
    doc.LinkEndChild( comment );

    TiXmlElement* estElement = new TiXmlElement( "Estimation" );

    estElement->SetAttribute("algorithm", algorithmToStringMap_[algorithm_]);
    estElement->SetAttribute("imgFile", imgFile_);
    estElement->SetAttribute("psfFile", psfFile_);
    estElement->SetAttribute("estFile", estFile_);

    estElement->SetAttribute("centeredPSF", (centeredPSF_ ? "1" : "0"));
    estElement->SetAttribute("version", version_);

    TiXmlElement* linearElement = new TiXmlElement( "Linear" );
    {
        element = new TiXmlElement( "Parameter" );
        element->SetDoubleAttribute("value", linearParameter_);
        linearElement->LinkEndChild( element);
    }
    estElement->LinkEndChild( linearElement);

    TiXmlElement* iterativeElement = new TiXmlElement( "Iterative" );
    {
        iterativeElement->SetAttribute("iterations", iterationNumber_);
        iterativeElement->SetAttribute("update", iterationUpdate_);
        iterativeElement->SetAttribute("write", iterationWrite_);

        TiXmlElement* invariantElement = new TiXmlElement( "Invariant" );
        {
            TiXmlElement* doubleZElement = new TiXmlElement( "DoubleZ" );
            {
                doubleZElement->SetAttribute("value", (invariantDoubleZ_ ? "1" : "0"));
            }
            invariantElement->LinkEndChild( doubleZElement);
            TiXmlElement* penaltyElement = new TiXmlElement( "Penalty" );
            {
                penaltyElement->SetAttribute("type", penaltyTypeToStringMap_[invariantPenaltyType_]);
                penaltyElement->SetDoubleAttribute("value", invariantPenaltyValue_);
            }
            invariantElement->LinkEndChild( penaltyElement);
        }
        iterativeElement->LinkEndChild( invariantElement);

        TiXmlElement* variantElement = new TiXmlElement( "Variant" );
        {
            TiXmlElement* doubleZElement = new TiXmlElement( "DoubleZ" );
            {
                doubleZElement->SetAttribute("value", (variantDoubleZ_ ? "1" : "0"));
            }
            variantElement->LinkEndChild( doubleZElement);
            TiXmlElement* penaltyElement = new TiXmlElement( "Penalty" );
            {
                penaltyElement->SetAttribute("type", penaltyTypeToStringMap_[variantPenaltyType_]);
                penaltyElement->SetDoubleAttribute("value", variantPenaltyValue_);
            }
            variantElement->LinkEndChild( penaltyElement);
            TiXmlElement* strataElement = new TiXmlElement( "Strata" );
            {
                strataElement->SetAttribute("value", variantStrataValue_);
                strataElement->SetAttribute("start", variantStrataStart_);
                strataElement->SetAttribute("size", variantStrataSize_);
            }
            variantElement->LinkEndChild( strataElement);
            TiXmlElement* otfElement = new TiXmlElement( "OTF" );
            {
                otfElement->SetAttribute("type", otfTypeToStringMap_[variantOtfType_]);
                otfElement->SetAttribute("name", variantOtfName_);
            }
        }
        iterativeElement->LinkEndChild( variantElement);

        TiXmlElement* convergenceElement = new TiXmlElement( "Convergence" );
        {
            convergenceElement->SetAttribute("name", convergenceName_);
		    int i = 0;
			TiXmlElement* convergenceTypeElement = 0;
		    for ( i = 0; i < 5; i++ )
			{
				int errorType = (0x1 << i);
			    if ( (convergenceType_ & errorType) != 0 )
				{
					convergenceTypeElement = new TiXmlElement( "ConvergenceType" );
					convergenceTypeElement->SetAttribute("value", errorTypeToStringMap_[ErrorType(errorType)]);
					convergenceElement->LinkEndChild( convergenceTypeElement );
				}
			}
        }
        iterativeElement->LinkEndChild( convergenceElement);
    }
    estElement->LinkEndChild( iterativeElement);
    doc.LinkEndChild( estElement );
    return doc.SaveFile( filename );
}
