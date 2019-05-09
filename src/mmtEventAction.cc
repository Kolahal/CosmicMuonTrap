#include "mmtEventAction.hh"
#include "mmtRunAction.hh"
//#include "mmtAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

mmtEventAction::mmtEventAction(): G4UserEventAction(), fEnergyAbs(0.), fEnergyGap(0.), fTrackLAbs(0.), fTrackLGap(0.)
{
}

mmtEventAction::~mmtEventAction()
{
}

void mmtEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	fEnergyAbs = 0.;
	fEnergyGap = 0.;
	fTrackLAbs = 0.;
	fTrackLGap = 0.;
}

void mmtEventAction::EndOfEventAction(const G4Event* event)
{
	auto eventID = event->GetEventID();
	auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	
	if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
	{
		G4cout << "---> End of event: " << eventID << G4endl;      
		
		G4cout
		<< "   Absorber: total energy: " << std::setw(7)
                << G4BestUnit(fEnergyAbs,"Energy")
		<< "       total track length: " << std::setw(7)
                << G4BestUnit(fTrackLAbs,"Length")
       		<< G4endl                        
       		<< "        Gap: total energy: " << std::setw(7)
                << G4BestUnit(fEnergyGap,"Energy")
       		<< "       total track length: " << std::setw(7)
                << G4BestUnit(fTrackLGap,"Length")
       		<< G4endl;                       
	}
}
