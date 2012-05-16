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

#include "tinyxml/tinyxml.h"
#include <iostream>
#define STR(s) #s
#define XSTR(s) STR(s)
#define VERSION XSTR(COSM_VERSION)

using namespace cosm;

template<typename T>
ToolsXml<T>::ToolsXml() {};

template<typename T>
ToolsXml<T>::~ToolsXml() {};

template<typename T>
bool ToolsXml<T>::saveConvolve(
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
    comment->SetValue(" Cosm tools of Convolve specification file ");
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue("     Output Image Size of Convolve    = < \"Normal Size\" | \"Size of 1st Image\" | \"Size of 2nd Image\" > " );
    doc.LinkEndChild( comment );


    TiXmlElement* toolsElement = new TiXmlElement( "Tools" );

    TiXmlElement* ConvolveElement = new TiXmlElement( "Convolve" );
    ConvolveElement->SetAttribute("imgFile", imgFile_);
    ConvolveElement->SetAttribute("psfFile", psfFile_);
    ConvolveElement->SetAttribute("outFile", outFile_);
    {
        element = new TiXmlElement( "OutputImageSize" );
        element->SetAttribute("value", outImgSize_);
        ConvolveElement->LinkEndChild( element);

        element = new TiXmlElement( "CenteredPsf" );
        element->SetAttribute("value", convolveCenteredPsf_);
        ConvolveElement->LinkEndChild( element);
    }
    toolsElement->LinkEndChild( ConvolveElement);

    doc.LinkEndChild( toolsElement );

    return doc.SaveFile( filename );
}

template<typename T>
bool ToolsXml<T>::saveVariant(
    const std::string& filename
){

    TiXmlDocument doc;
    TiXmlComment* comment;
    TiXmlElement* element;

    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

	comment = new TiXmlComment();
	std::string v = "Cosm Version: ";
	v += VERSION;
    comment->SetValue(v);
    doc.LinkEndChild( comment );

    comment = new TiXmlComment();
    comment->SetValue(" Cosm tools of Variant specification file ");
    doc.LinkEndChild( comment );


    TiXmlElement* toolsElement = new TiXmlElement( "Tools" );

    TiXmlElement* variantElement = new TiXmlElement( "Variant" );
    variantElement->SetAttribute("imgFile", imgFile2_);
    variantElement->SetAttribute("psfFile", psfFile2_);
    variantElement->SetAttribute("outFile", outFile2_);
    {
		TiXmlElement* strataElement = new TiXmlElement( "Strata" );
		{
			strataElement->SetAttribute("value", variantStrataValue_);
			strataElement->SetAttribute("start", variantStrataStart_);
			strataElement->SetAttribute("size", variantStrataSize_);
		}
		variantElement->LinkEndChild( strataElement);

        element = new TiXmlElement( "CenteredPsf" );
        element->SetAttribute("value", variantCenteredPsf_);
        variantElement->LinkEndChild( element);
    }
    toolsElement->LinkEndChild( variantElement);

    doc.LinkEndChild( toolsElement );

    return doc.SaveFile( filename );
}
