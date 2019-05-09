#ifndef mmtSteppingAction_h
#define mmtSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class mmtDetectorConstruction;
class mmtEventAction;

/// Stepping action class.
//In UserSteppingAction() there are collected the energy deposit and track 
//lengths of charged particles in Absober and Gap layers and
//updated in mmtEventAction.

class mmtSteppingAction : public G4UserSteppingAction
{
	public:
		mmtSteppingAction(const mmtDetectorConstruction* detectorConstruction, mmtEventAction* eventAction);
		virtual ~mmtSteppingAction();

		virtual void UserSteppingAction(const G4Step* step);

	private:
		const mmtDetectorConstruction* fDetConstruction;
		mmtEventAction*  fEventAction;
};

#endif
