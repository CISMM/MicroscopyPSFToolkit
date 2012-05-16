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

#include "toolMgr.h"
#include "wu/wuHeader.h"
#include <tclap/CmdLine.h>
#include <string>
#include <iostream>

#define XSTR(s) STR(s)
#define STR(s) #s
#define VERSION XSTR(COSM_VERSION)

using namespace TCLAP;
using namespace cosm;

int main( int argc, char* argv[] ) {

    // Define command line object 
    CmdLine cmdLine("COSM Tools",' ', VERSION);
 
    // Define the argument options
    ValueArg<std::string> dataArg("i","input", "Input Image", false, "img.wu", "image");
    ValueArg<std::string> refArg("g","reference", "Reference Image", false, "ref.wu", "reference");
    ValueArg<std::string> outputArg("o","output", "Output Image", false, "out.wu", "output");
    ValueArg<int> xDimArg("x","xdim", "X Dimension", false, 64, "xdim");
    ValueArg<int> yDimArg("y","ydim", "Y Dimension", false, 64, "ydim");
    ValueArg<int> zDimArg("z","zdim", "Z Dimension", false, 64, "zdim");
    ValueArg<int> xShiftArg("j","xshift", "X Shift", false, 31, "xshift");
    ValueArg<int> yShiftArg("k","yshift", "Y Shift", false, 31, "yshift");
    ValueArg<int> zShiftArg("l","zshift", "Z Shift", false, 31, "zshift");
    ValueArg<int> xCenterArg("m","xcenter", "X Center", false, 31, "xcenter");
    ValueArg<int> yCenterArg("n","ycenter", "Y Center", false, 31, "ycenter");
    ValueArg<int> zCenterArg("q","zcenter", "Z Center", false, 31, "zcenter");
    ValueArg<int> xRadiusArg("r","xradius", "X Radius", false, 8, "xradius");
    ValueArg<int> yRadiusArg("s","yradius", "Y Radius", false, 8, "yradius");
    ValueArg<int> zRadiusArg("w","zradius", "Z Radius", false, 8, "zradius");
    ValueArg<double> valueArg("v","value", "Fill Value", false, 0, "fill value");
    ValueArg<int> typeArg("t","type", "Data Type", false, 0, "data type");
    SwitchArg endianArg("e","endian", "Big Endian", false);

    SwitchArg boxArg("b","box", "Box", false);
    SwitchArg pointArg("p","point", "Point", false);
    SwitchArg ellipsoidArg("c","ellipsoid", "Ellipsoid", false);

    SwitchArg infoArg("I","info", "Image Information", false);
    SwitchArg addArg("A","add", "Add Image Header", false);
    SwitchArg removeArg("R","remove", "Remove Image Header", false);
    SwitchArg transformArg("T","tranform", "Transform Image", false);
    ValueArg<bool> resampleArg("N","resample", "Resample Image", false, false, "up/down");
    SwitchArg convertArg("U","convert", "Convert Image Type", false);
    ValueArg<bool> shiftArg("S","shift", "Shift Image", false, false, "normal/circular");
    ValueArg<int> scaleArg("X","scale", "Scale Image", false, 2, "0-sum, 1-max, 2-value");
    SwitchArg compareArg("M","convolve", "Compare Images", false);
    SwitchArg convolveArg("C","convolve", "Convolve Images", false);
    SwitchArg objectArg("O","object", "Create Object Image", false);

    std::vector<Arg*> objArg(3);
    objArg[0] = &ellipsoidArg;
    objArg[1] = &boxArg;
    objArg[2] = &pointArg;
    cmdLine.xorAdd(objArg);

    std::vector<Arg*> toolArg(8);
    toolArg[0] = &infoArg;
    toolArg[1] = &addArg;
    toolArg[2] = &removeArg;
    toolArg[3] = &transformArg;
    toolArg[4] = &resampleArg;
    toolArg[5] = &convertArg;
    toolArg[6] = &convolveArg;
    toolArg[7] = &objectArg;
    cmdLine.xorAdd(toolArg);

    // add arguments to command line options
    cmdLine.add(dataArg);
    cmdLine.add(refArg);
    cmdLine.add(outputArg);
    cmdLine.add(xDimArg);
    cmdLine.add(yDimArg);
    cmdLine.add(zDimArg);
    cmdLine.add(xShiftArg);
    cmdLine.add(yShiftArg);
    cmdLine.add(zShiftArg);
    cmdLine.add(xCenterArg);
    cmdLine.add(yCenterArg);
    cmdLine.add(zCenterArg);
    cmdLine.add(xRadiusArg);
    cmdLine.add(yRadiusArg);
    cmdLine.add(zRadiusArg);
    cmdLine.add(valueArg);
    cmdLine.add(typeArg);
    cmdLine.add(endianArg);

    // parse command line
    cmdLine.parse(argc, argv);

    int xDim = xDimArg.getValue();
    int yDim = yDimArg.getValue();
    int zDim = zDimArg.getValue();

    int xShift = xShiftArg.getValue();
    int yShift = yShiftArg.getValue();
    int zShift = zShiftArg.getValue();

    int xCenter = xCenterArg.getValue();
    int yCenter = yCenterArg.getValue();
    int zCenter = zCenterArg.getValue();

    int xRadius = xRadiusArg.getValue();
    int yRadius = yRadiusArg.getValue();
    int zRadius = zRadiusArg.getValue();

    double value = valueArg.getValue();
    wuHeader::type type = (wuHeader::type)typeArg.getValue();
    bool endian = endianArg.getValue();

    if ( infoArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        wuHeader header;
        bool res = toolMgr.info(header);
        if ( res )
        {
            std::cout <<"Image Information:"<< std::endl;
            std::cout <<"  X: "<<header.columns()<<", Y: "<<header.rows()<<", Z: "<< header.sections() << std::endl;
            std::cout <<"  Maximum: "<<header.maximum()<<", Minimum: "<<header.minimum() <<", Mean: "<< header.mean() << std::endl;
            std::cout <<"  Data Type: "<<header.dataTypeStr()<<", Data Size: "<< header.dataSize() << std::endl;
            return 0;
        }
    }
    if ( addArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.addHeader(xDim, yDim, zDim, type, endian) )
        {
            return 0;
        }
    }
    if ( removeArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.removeHeader() ) 
        {
            return 0;
        }
    }
    if ( resampleArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.resample(xDim, yDim, zDim, resampleArg.getValue()) ) 
        {
            return 0;
        }
    }
    if ( transformArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.transformation(xDim, yDim, zDim, xShift, yShift, zShift, value) ) 
        {
            return 0;
        }
    }
    if ( convertArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.convert(type) ) 
        {
            return 0;
        }
    }
    if ( shiftArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.shift(xShift, yShift, zShift, shiftArg.getValue()) )
        {
            return 0;
        }
    }
    if ( scaleArg.getValue() )
    {
        bool sumFlag = scaleArg.getValue() == 0;
        bool maxFlag = scaleArg.getValue() == 1;
        double val = value != 0.0 ? value : 1.0;
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.scale( sumFlag, maxFlag, val) )
        {
            return 0;
        }
    }
    if ( compareArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), refArg.getValue(), outputArg.getValue() );
        double maximum, mean, meanSquare;
        if ( toolMgr.compare(maximum, mean, meanSquare) )
        {
            return 0;
        }
        std::cout <<"Maximum: "<< maximum <<", Mean: "<< mean <<", Mean square: "<< meanSquare << std::endl;
    }
    if ( convolveArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), refArg.getValue(), outputArg.getValue() );
        if ( toolMgr.convolve() )
        {
            return 0;
        }
    }
    if ( objectArg.getValue() && !(ellipsoidArg.getValue() || boxArg.getValue() || pointArg.getValue()) )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.create(xDim, yDim, zDim, type, value) )
        {
            return 0;
        }
    }
    if ( objectArg.getValue() && ellipsoidArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.ellipsoid(xCenter, yCenter, zCenter, xRadius, yRadius, zRadius, value) )
        {
            return 0;
        }
    }
    if ( objectArg.getValue() && boxArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.box(xCenter, yCenter, zCenter, xRadius, yRadius, zRadius, value) )
        {
            return 0;
        }
    }
    if ( objectArg.getValue() && pointArg.getValue() )
    {
        ToolMgr toolMgr( dataArg.getValue(), outputArg.getValue() );
        if ( toolMgr.point(xCenter, yCenter, zCenter, value) )
        {
            return 0;
        }
    }
    return -1;
}
