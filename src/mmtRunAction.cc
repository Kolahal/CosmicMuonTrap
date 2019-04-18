#include "mmtRunAction.hh"
//#include "mmtAnalysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

mmtRunAction::mmtRunAction(): G4UserRunAction()
{
	// set printing event number per each event
	G4RunManager::GetRunManager()->SetPrintProgress(1);
}
mmtRunAction::~mmtRunAction()
{
}
/*
void mmtRunAction::BeginOfRunAction(const G4Run* //run//)
{
}
void mmtRunAction::EndOfRunAction(const G4Run* //run//)
{
}
*/
