#include "mmtActionInitialization.hh"
#include "mmtPrimaryGeneratorAction.hh"
#include "mmtRunAction.hh"
#include "mmtEventAction.hh"
#include "mmtSteppingAction.hh"
#include "mmtDetectorConstruction.hh"

mmtActionInitialization::mmtActionInitialization(mmtDetectorConstruction* detConstruction): G4VUserActionInitialization(),fDetConstruction(detConstruction)
{
}
mmtActionInitialization::~mmtActionInitialization()
{
}
void mmtActionInitialization::BuildForMaster() const
{
	SetUserAction(new mmtRunAction);
}

void mmtActionInitialization::Build() const
{
	SetUserAction(new mmtPrimaryGeneratorAction);
	SetUserAction(new mmtRunAction);
	auto eventAction = new mmtEventAction;
	SetUserAction(eventAction);
	SetUserAction(new mmtSteppingAction(fDetConstruction,eventAction));
}
