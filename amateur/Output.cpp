/***************************************************************************
 *   Copyright (C) 2010 by Pascal Stephan Philipp Steger                   *
 *   psteger@phys.ethz.ch                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, ma  02111-1307, USA.             *
 ***************************************************************************/

/**
 * \file Output.cpp
 * \brief all output routines
 * provides routines for accessing and writing AHF output as well as
 * calculated Halo outputs.
 */

#include <iostream>
#include <vector>
#include <fstream>
#include "Output.h"


Output::Output( void ){
  std::cerr << "Trying to create output without a simulation" << std::endl;
}

Output::Output(Simulation sim){
  sim_ = sim;
}

Output::~Output( void ){
  ar_.clear();
  hr_.clear();
}

void
Output::getAR( void ) {
  std::cout << " in getAR..." << std::endl;
  unsigned npart, nvpart;
  float Xc, Yc, Zc, VXc, VYc, VZc;
  float Mvir, Rvir, Vmax, Rmax, sigV;
  float lambda, Jx, Jy, Jz;
  float a, Eax, Eay, Eaz;
  float b, Ebx, Eby, Ebz;
  float c, Ecx, Ecy, Ecz;
  float ovdens, Redge;
  int nbins;
  float Ekin, Epot;
  float mbp_offset, com_offset;
  float r2, lambdaE;

  std::fstream ahfhalos((Folder("output") + Snaps(sim_.getIsnap()) 
			 + "/ahf_out_halos2").c_str(), std::ios::in);
  for(int i=0; i<sim_.getNhalo(); ++i){
    npart = nvpart = 0;
    Xc=Yc=Zc=VXc=VYc=VZc=Mvir=Rvir=Vmax=Rmax=sigV=lambda=Jx=Jy=Jz=0.0;
    a=Eax=Eay=Eaz=b=Ebx=Eby=Ebz=c=Ecx=Ecy=Ecz=ovdens=Redge=0.0;
    nbins=0;
    Ekin=Epot=mbp_offset=com_offset=r2=lambdaE=0.0;
    ahfhalos >> npart >> nvpart>> Xc    >> Yc   >> Zc
	     >> VXc   >> VYc   >> VZc
	     >> Mvir  >> Rvir
	     >> Vmax  >> Rmax
	     >> sigV
	     >> lambda
	     >> Jx    >> Jy    >> Jz
	     >> a     >> Eax   >> Eay   >> Eaz
	     >> b     >> Ebx   >> Eby   >> Ebz
	     >> c     >> Ecx   >> Ecy   >> Ecz
	     >> ovdens
	     >> Redge
	     >> nbins
	     >> Ekin   >> Epot 
	     >> mbp_offset >> com_offset 
	     >> r2 
	     >> lambdaE;
    AResults arn = AResults(npart,nvpart,Xc,Yc,Zc,VXc,VYc,VZc,
			    Mvir,Rvir,Vmax,Rmax,sigV,
			    lambda,Jx,Jy,Jz,
			    a,Eax,Eay,Eaz,b,Ebx,Eby,Ebz,c,Ecx,Ecy,Ecz,
			    ovdens,Redge,nbins,Ekin,Epot,
			    mbp_offset,com_offset,
			    r2,lambdaE);
    ar_.push_back(arn);
  }
  ahfhalos.close();
}

void
Output::writeAR( void ){
  std::cout << " in writeAR()..." << std::endl;
  std::string filename = Folder("output")+Snaps(sim_.getIsnap())
    +"/Aprop.txt";
  std::ofstream file(filename.c_str());
  for(int i=0; i<sim_.getNhalo(); ++i){
    file << ar_.at(i).toString() << std::endl;
  }
  file.close();
}

void
Output::fillHRglobal( int Ntot, real Mtot, real Rtot, real Mvir, real Rvir,
		      Vector3 xcm, Vector3 vcm, Vector3 xmb, Vector3 vmb ){
  tmp_.setNtot(Ntot);
  tmp_.setMtot(Mtot);
  tmp_.setRtot(Rtot);
  tmp_.setMvir(Mvir);
  tmp_.setRvir(Rvir);
  tmp_.setXCM(xcm);
  tmp_.setVCM(vcm);
  tmp_.setXMB(xmb);
  tmp_.setVMB(vmb);
}

void 
Output::fillHRtype( int type, int npart, real mpart, real lambda, 
		    real sigma, Vector3 J, Eigensystem3 E ){
  tmp_.setNpart(npart,type);
  tmp_.setMpart(mpart,type);
  tmp_.setJ(J,type);
  tmp_.setE(E,type);
  tmp_.setLambda(lambda,type);
  tmp_.setSigma(sigma,type);
}

void
Output::appendHR( void ){
  hr_.push_back(tmp_);
}

void
Output::writeHR( void ) {
  int nhalo = sim_.getNhalo();
  std::string filename = Folder("output") + Snaps(sim_.getIsnap()) 
    + "/Hprop.txt";
  std::ofstream file(filename.c_str());

  file << sim_.getBoxlength() << " " << sim_.getZsnap() 
       << " " << nhalo << std::endl;
  for(int i=0;i<nhalo;++i){
    file << i << " " << hr_.at(i).toString() << std::endl;
  }
  file.close();
}

void
Output::writeHRstripped( void ) {
  int nhalo = sim_.getNhalo();
  std::string filename = Folder("output") + Snaps(sim_.getIsnap()) 
    + "/HpropStripped.txt";
  std::ofstream file(filename.c_str());

  for(int i=0;i<nhalo;++i){
    file << i << " " << hr_.at(i).toString() << std::endl;
  }
  file.close();
}
