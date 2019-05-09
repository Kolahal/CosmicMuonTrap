#include "mmtMagField.hh"
#include "g4hdf5.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
using namespace std;

mmtMagneticField::mmtMagneticField(): G4MagneticField(), fMessenger(nullptr), fBz(1.0*tesla)
{
	//open the magnetic field map file here
	ifstream Mag_in;
	Mag_in.open("../MagFieldMap.txt");
	/*
	double xi = 0.0;
	double yj = 0.0;
	double zk = 0.0;
	double Bx = 0.0;
	double By = 0.0;
	double Bz = 0.0;
	
	while(!Mag_in.eof())
	{
		Mag_in>>xi>>yj>>zk>>Bx>>By>>Bz;
	}
	*/
}

mmtMagneticField::~mmtMagneticField()
{
}

void mmtMagneticField::GetFieldValue(const G4double* point,double *bField) const
{
	bField[0] = 0.;
	bField[1] = 0.;
	bField[2] = fBz;
}
