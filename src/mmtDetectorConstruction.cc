#include "mmtDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4GlobalMagFieldMessenger* mmtDetectorConstruction::fMagFieldMessenger = nullptr;

mmtDetectorConstruction::mmtDetectorConstruction():G4VUserDetectorConstruction(),fCheckOverlaps(true)
{
}

mmtDetectorConstruction::~mmtDetectorConstruction()
{
}

G4VPhysicalVolume* mmtDetectorConstruction::Construct()
{
	DefineMaterials();
	return DefineVolumes();
}
void mmtDetectorConstruction::DefineMaterials()
{
	// Lead material defined using NIST Manager
	auto nistManager = G4NistManager::Instance();
	
	//G4Material* vac= new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
	G4Material* air=nistManager->FindOrBuildMaterial("G4_AIR");
	G4Material* pps=nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	G4cout<<"Kolahal: MMTDetectorConstruction::DefineMaterials() -> "<< *(G4Material::GetMaterialTable()) << G4endl;
}
G4VPhysicalVolume* mmtDetectorConstruction::DefineVolumes()
{
        // Geometry parameters
	G4int nofLayers = 1;
        G4double airThickness = 100.0*cm;
        G4double ppsThickness = 1.0*cm;
        G4double calorSizeXY  = 100.0*cm;
	
        auto layerThickness = ppsThickness + airThickness + ppsThickness;
        auto calorThickness = nofLayers * layerThickness;
        auto worldSizeXY = 1.2 * calorSizeXY;
        auto worldSizeZ  = 1.2 * calorThickness;
	
        // Get materials
        //auto defaultMaterial = G4Material::GetMaterial("Galactic");
	auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
	auto MagTrapMaterial = G4Material::GetMaterial("G4_AIR");
	auto ScintlrMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
	
	 
	G4RotationMatrix* Rot = new G4RotationMatrix;
	Rot->rotateX(M_PI/2.*rad);
	
	//World:
	auto worldS = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2);
	auto worldLV= new G4LogicalVolume(worldS, defaultMaterial, "World");
	auto worldPV= new G4PVPlacement(0, G4ThreeVector(), worldLV, "World", 0, false, 0, fCheckOverlaps);
	
	// Calorimeter
	auto calorimeterS= new G4Box("Calorimeter",calorSizeXY/2, calorSizeXY/2, calorThickness/2);
	auto calorLV	 = new G4LogicalVolume(calorimeterS, defaultMaterial,"Calorimeter");
	auto calorPV	 = new G4PVPlacement(Rot, G4ThreeVector(), calorLV, "Calorimeter", worldLV, false, 0, fCheckOverlaps);
	
	// Layer
	auto layerS = new G4Box("Layer", calorSizeXY/2, calorSizeXY/2, layerThickness/2);
	auto layerLV= new G4LogicalVolume(layerS, defaultMaterial, "Layer");
	auto replica= new G4PVReplica("Layer", layerLV, calorLV, kZAxis, nofLayers, layerThickness);
	
	// Air
	auto mmtAirS 	= new G4Box("mmtAir", calorSizeXY/2, calorSizeXY/2, airThickness/2);
	auto mmtAirLV	= new G4LogicalVolume(mmtAirS, MagTrapMaterial, "mmtAir");
	mmtAirPV	= new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), mmtAirLV, "mmtAir", layerLV, false, 0, fCheckOverlaps);
	
	// Scintillator - up
	auto ppsPlate1S	= new G4Box("pps_up", calorSizeXY/2, calorSizeXY/2, ppsThickness/2);
	auto ppsPlate1LV= new G4LogicalVolume(ppsPlate1S, ScintlrMaterial, "pps_up");
	ppsPlate1PV	= new G4PVPlacement(0, G4ThreeVector(0., 0., (airThickness/2+ppsThickness/2)), ppsPlate1LV, "pps_up", layerLV, false, 0, fCheckOverlaps);
	
	// Scintillator - dn
	auto ppsPlate2S = new G4Box("pps_dn", calorSizeXY/2, calorSizeXY/2, ppsThickness/2);
	auto ppsPlate2LV= new G4LogicalVolume(ppsPlate2S, ScintlrMaterial, "pps_dn");
	ppsPlate2PV     = new G4PVPlacement(0, G4ThreeVector(0., 0.,-(airThickness/2+ppsThickness/2)), ppsPlate2LV, "pps_dn", layerLV, false, 0, fCheckOverlaps);
	
  	G4cout
    	<< G4endl
    	<< "------------------------------------------------------------" << G4endl
    	<< "---> The Magnetic mirror trap is " << nofLayers << " layers of: [ "
	<< ppsThickness/cm << "mm of " << ScintlrMaterial->GetName()
	<< " + "
    	<< airThickness/cm << "mm of " << MagTrapMaterial->GetName()
    	<< " + "
    	<< ppsThickness/cm << "mm of " << ScintlrMaterial->GetName() << " ] " << G4endl
    	<< "------------------------------------------------------------" << G4endl;
	
	// Visualization attributes
	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
	auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	simpleBoxVisAtt->SetVisibility(true);
	calorLV->SetVisAttributes(simpleBoxVisAtt);
	
	return worldPV;
}
void mmtDetectorConstruction::ConstructSDandField()
{
	// Create global magnetic field messenger.
	G4ThreeVector fieldValue;
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);
	
	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}
