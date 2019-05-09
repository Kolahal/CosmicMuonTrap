#ifndef mmtEventAction_h
#define mmtEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class mmtEventAction : public G4UserEventAction
{
	public:
	mmtEventAction();
	virtual ~mmtEventAction();
	
	virtual void BeginOfEventAction(const G4Event* event);
	virtual void EndOfEventAction(const G4Event* event);
	
	void AddAbs(G4double de, G4double dl);
	void AddGap(G4double de, G4double dl);
	
	private:
	G4double  fEnergyAbs;
	G4double  fEnergyGap;
	G4double  fTrackLAbs;
	G4double  fTrackLGap;
};

inline void mmtEventAction::AddAbs(G4double de, G4double dl) {
	fEnergyAbs += de;
	fTrackLAbs += dl;
}

inline void mmtEventAction::AddGap(G4double de, G4double dl) {
	fEnergyGap += de;
	fTrackLGap += dl;
}

#endif
