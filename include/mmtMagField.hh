#ifndef mmtMagField_HH
#define mmtMagField_HH

#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4MagneticField.hh"

class G4GenericMessenger;

class mmtMagneticField : public G4MagneticField
{
	public:
		mmtMagneticField();
		virtual ~mmtMagneticField();

		virtual void GetFieldValue(const G4double point[4],double* bField ) const;

		void SetField(G4double val) { fBy = val; }
		G4double GetField() const { return fBy; }

	private:
		//void DefineCommands();

		G4GenericMessenger* fMessenger;
		
		G4double fBx;
		G4double fBy;
		G4double fBz;
};
#endif
