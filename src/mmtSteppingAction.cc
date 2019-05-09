#include "mmtSteppingAction.hh"
#include "mmtEventAction.hh"
#include "mmtDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

mmtSteppingAction::mmtSteppingAction(const mmtDetectorConstruction* detectorConstruction,mmtEventAction* eventAction): G4UserSteppingAction(), fDetConstruction(detectorConstruction), fEventAction(eventAction)
{
}

mmtSteppingAction::~mmtSteppingAction()
{
}

void mmtSteppingAction::UserSteppingAction(const G4Step* step)
{
	// Collect energy and track length step by step
	
	// get volume of the current step
	auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	
	// energy deposit
	auto edep = step->GetTotalEnergyDeposit();
	
	// step length
	G4double stepLength = 0.;
	if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. )
	{
		stepLength = step->GetStepLength();
	}
	
	if ( volume == fDetConstruction->GetMMTAirPV() )
	{
		fEventAction->AddAbs(edep,stepLength);
	}
	
	if ( volume == fDetConstruction->GetPPS1PV() )
	{
		fEventAction->AddGap(edep,stepLength);
	}
	
	if ( volume == fDetConstruction->GetPPS1PV() )
	{
		fEventAction->AddGap(edep,stepLength);
	}
}
