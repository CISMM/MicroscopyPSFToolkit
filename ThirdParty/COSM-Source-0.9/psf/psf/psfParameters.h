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
                                                                                
#ifndef _PSF_PARAMETERS_H
#define _PSF_PARAMETERS_H

namespace cosm {

template <typename T>
class PsfParameters {

public:
    PsfParameters() : 
         Nxy_(0),          // size of the image in x and y
         Nz_(0),           // size of the image in z
         deltaXY_(0),      // pixel size in z and y of the image (mm)
         deltaZ_(0),       // pixel size in z of the image (mm)
         ts_(0),           // specimen thickness
         tid_(0),          // immersion thickness design
         tia_(0),          // immersion thickness actual
         tgd_(0),          // coverglass thickness design
         tga_(0),          // coverglass thickness actual
         ns_(0),           // specimen refractive index
         nid_(0),          // immersion refractive index design
         nia_(0),          // immersion refractive index actual
         ngd_(0),          // coverglass refractive index design
         nga_(0),          // coverglass refractive index actual
         tld_(0),          // tube length design
         tla_(0),          // tube length actual
         lm_(0),           // lateral magnification
         na_(0),           // numerical aperture (NA)
         lambda_(0),
         absError_(0),
         fsize_(0),
         distance_(0),
         magY_(0),
         shear_(0),
         bias_(0),
         amplitudeRatio_(0),
         rotation_(0) {};

    virtual ~PsfParameters() {};
   
    void Nxy( int val ) { Nxy_ = val; };
    int Nxy( void ) const { return Nxy_; };
    void Nz( int val ) { Nz_ = val; };
    int Nz( void ) const { return Nz_; };
    void deltaXY( T val ) { deltaXY_ = val; };
    T deltaXY( void ) const { return deltaXY_; };
    void deltaZ( T val ) { deltaZ_ = val; };
    T deltaZ( void ) const { return deltaZ_; };
    void ts( T val ) { ts_ = val; };
    T ts( void ) const { return ts_; };
    void tid( T val ) { tid_ = val; };
    T tid( void ) const { return tid_; };
    void tia( T val ) { tia_ = val; };
    T tia( void ) const { return tia_; };
    void tgd( T val ) { tgd_ = val; };
    T tgd( void ) const { return tgd_; };
    void tga( T val ) { tga_ = val; };
    T tga( void ) const { return tga_; };
    void ns( T val ) { ns_ = val; };
    T ns( void ) const { return ns_; };
    void nid( T val ) { nid_ = val; };
    T nid( void ) const { return nid_; };
    void nia( T val ) { nia_ = val; };
    T nia( void ) const { return nia_; };
    void ngd( T val ) { ngd_ = val; };
    T ngd( void ) const { return ngd_; };
    void nga( T val ) { nga_ = val; };
    T nga( void ) const { return nga_; };
    void tld( T val ) { tld_ = val; };
    T tld( void ) const { return tld_; };
    void tla( T val ) { tla_ = val; };
    T tla( void ) const { return tla_; };
    void lm( T val ) { lm_ = val; };
    T lm( void ) const { return lm_; };
    void na( T val ) { na_ = val; };
    T na( void ) const { return na_; };
    void lambda( T val ) { lambda_ = val; };
    T lambda( void ) const { return lambda_; };
    void absError( T val ) { absError_ = val; };
    T absError( void ) const { return absError_; };
	void fsize( T val ) { fsize_ = val; };
	T fsize( void ) const { return fsize_; };
	void distance( T val ) { distance_ = val; };
	T distance( void ) const { return distance_; };
	void magY( T val ) { magY_ = val; };
	T magY( void ) const { return magY_; };
	void shear( T val ) { shear_ = val; };
	T shear( void ) const { return shear_; };
	void bias( T val ) { bias_ = val; };
	T bias( void ) const { return bias_; };
	void amplitudeRatio( T val ) { amplitudeRatio_ = val; };
	T amplitudeRatio( void ) const { return amplitudeRatio_; };
	void rotation( T val ) { rotation_ = val; };
	T rotation( void ) const { return rotation_; }

protected:
    // not allowed
    PsfParameters(const PsfParameters&);
    // not allowed
    const PsfParameters& operator=(const PsfParameters&);

protected:
    int Nxy_;        // size of the image in x and y
    int Nz_;         // size of the image in z
    T deltaXY_;      // pixel size in z and y of the image (mm)
    T deltaZ_;       // pixel size in z of the image (mm)
    T ts_;           // specimen thickness
    T tid_;          // immersion thickness design
    T tia_;          // immersion thickness actual
    T tgd_;          // coverglass thickness design
    T tga_;          // coverglass thickness actual
    T ns_;           // specimen refractive index
    T nid_;          // immersion refractive index design
    T nia_;          // immersion refractive index actual
    T ngd_;          // coverglass refractive index design
    T nga_;          // coverglass refractive index actual
    T tld_;          // tube length design
    T tla_;          // tube length actual
    T lm_;           // lateral magnification
    T na_;           // numerical aperture (NA)
    T lambda_;
    T absError_;
	T fsize_;
	T distance_;
	T magY_;
	T shear_;
	T bias_;
	T amplitudeRatio_;
	T rotation_;
};

}

#endif // _PSF_PARAMETERS_H;
