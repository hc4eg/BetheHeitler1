#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "MagneticField.hh"
#include "G4NistManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "WireChamberSD.hh"
#include "PaddleSD.hh"
#include "MonitorSD.hh"

//shapes
class G4Box;
class G4Tubs;
class G4SubtractionSolid;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4ExtrudedSolid;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

//struct Section;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:

	G4bool fVerbose;
	void SetMagField(G4LogicalVolume * , G4double);
     
	G4VPhysicalVolume* Construct();

	void UpdateGeometry();
    	void DumpGeometricalTree(G4VPhysicalVolume*, G4int);
                 
	inline G4VPhysicalVolume* GetPhysWorld() {return fPhysWorld;};

	G4NistManager* fNistMan; //used for getting materials from Geant4 database                 
	// a few methods to get variables that might be of use to other classes
	inline G4double GetTargetDistance() {return fTargetDistance;};
	inline G4double GetTargetThickness() {return fTargetThick;};
	inline G4String GetTargetMaterialName() {return fTargetMaterial->GetName();};
	inline G4double GetAirTargetThick() {return fAirTargetThick;};
	inline G4double GetVDCSpacing() {return fVDCSpacing;};
	inline G4double GetTargetContainerThick() {return fTargetContainerThick;};
	inline G4double GetMagneticField() {return fBcenter;};

	// methods to set target parameters
	inline void SetVerbose(G4bool val) { fVerbose = val;};
	inline void SetTargetThickness(G4double val) { fTargetThick = val;};
	inline void SetB0(G4double field) { fBcenter = field; if(fMagField) fMagField->SetFieldB0(field); };
	inline void SetTargetMaterial(G4String mat) { fTargetMaterial = G4Material::GetMaterial(mat);};
	inline void SetTargetMaterial( G4Material * mat) { fTargetMaterial = mat; }
	void Print();

  private:
     
	G4Material * fBackgroundMaterial;

	G4ThreeVector * fWorldSize;
	G4ThreeVector * fWorldCenter;

	// Magnetic Field
	G4double fBcenter;

	// Materials
	G4Material * Helium;
	G4Material * Vacuum;
	G4Material * Aluminum;
	G4Material * Copper;
	G4Material * Poly;
	G4Material * Mylar;
	G4Material * Uranium;
	G4Material * Air;
	G4Material * Sci;
	G4Material * Iron;
	G4Material * Ethane_Argon;
	G4Material * Carbon;

	//coil
	G4double fCoilInnerDiameter;
	G4double fCoilOuterDiameter;
	G4double fCoilHeight;
	G4Tubs * fSolidCoil;
	G4LogicalVolume * fLogicCoil;
	G4VPhysicalVolume * fPhysCoil;
	G4Material * fCoilMaterial;

	//target
	G4double fTargetDistance;
	G4double fTargetThick;
	G4double fTargetWidth;
	G4double fTargetHeight;
	G4double fTargetContainerThick;

	G4LogicalVolume * fLogicTargetContainer;
	G4LogicalVolume * fLogicTarget;
	G4LogicalVolume * fAirContainerLogic;
	G4Material * fTargetMaterial;
	G4Material * fTargetToBagMaterial;
	G4VPhysicalVolume * fPhysTarget;

	//helium bag
	G4Material * fBagMaterial;
	G4double fAirTargetThick;
	G4double fMagnetToBagDistance;
	G4double fBagThick;
	G4LogicalVolume * fLogicBagWall;

	//world
	G4Box*             fSolidWorld;    //pointer to the solid World 
	G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
	G4VPhysicalVolume* fPhysWorld;    //pointer to the physical World

	///magnet
	G4ThreeVector * fMagnetSize;
	G4Box * fSolidMagnet;
	G4LogicalVolume * fLogicMagnet;
	G4VPhysicalVolume * fPhysMagnet;

	//pole
	G4double fPoleOuterDiameter;
	G4double fPoleHeight;
	G4double fPoleBevelHeight;
	G4double fPoleBevelledDiameter;
	G4Tubs * fSolidPoleBase;
	G4LogicalVolume * fLogicPoleBase;
	G4VPhysicalVolume * fPhysPoleBase;
	G4Cons * fSolidPoleBevel;
	G4LogicalVolume * fLogicPoleBevel;
	G4VPhysicalVolume * fPhysPoleBevel;
	G4Material * fPoleMaterial;
	
	G4double fMagnetX, fMagnetY, fMagnetZ;
	G4double fYokeInnerX, fYokeInnerZ;
	//yoke
	
	G4double fYokeSideX;
	G4double fYokeSideY;
	G4double fYokeSideZ;
	G4double fYokeTopX;
	G4double fYokeTopY;
	G4double fYokeTopZ;
	G4double fYokeBeamHoleDiameter;

	G4SubtractionSolid * fYokeWithHole;	
	G4Tubs * fSolidYokeBeamHole;
	G4LogicalVolume * fLogicYokeBeamHole;
	G4VPhysicalVolume * fPhysYokeBeamHole;
	G4Box * fSolidYokeUpstreamSide;
	G4LogicalVolume * fLogicYokeUpstreamSide;
	G4VPhysicalVolume * fPhysYokeUpstreamSide;
	G4Box * fSolidYokeDownstreamSide;
	G4LogicalVolume * fLogicYokeDownstreamSide;
	G4VPhysicalVolume * fPhysYokeDownstreamSide;
	G4Box * fSolidYokeTop;
	G4LogicalVolume * fLogicYokeTop;
	G4VPhysicalVolume * fPhysYokeTop;
	G4Box * fSolidYokeBottom;
	G4LogicalVolume * fLogicYokeBottom;
	G4VPhysicalVolume * fPhysYokeBottom;
	G4Material * fYokeMaterial;
     
	//wire chambers
	G4LogicalVolume * fLogicChamberHolder;
	G4LogicalVolume * fLogicVDCHolder;
	G4LogicalVolume * fLogicVDCsheet;
	G4LogicalVolume * fLogicVDCgas1;
	G4LogicalVolume * fLogicVDCgasSD[2];
	G4LogicalVolume * fLogicVDCBag;

	//Sensitive Detectors
	WireChamberSD * fchamberSD;
	PaddleSD * fPaddleSD;
	MonitorSD * fMonitorSD;

	G4VPhysicalVolume * fPhysVDCGasSD[2];
	G4VPhysicalVolume * fPhysVDCpackage[2];

	G4Material * fChamberGasMaterial;
	G4Material * fChamberSheetMaterial;
	//G4Material * fWireMaterial;

	// dimensions
	G4double fVDCSizeX, fVDCSizeY, fVDCSizeZ ;
	G4double fVDCSideWidth, fVDCTopWidth;
	G4double fVDCSpacing;
	G4double fVDCSheetThick;
	G4double fVDCLayerThick;
	G4double fVDCAlLayerThick;
	G4double fVDCChamberHolderX;
	G4double fVDCChamberHolderY;
	G4double fVDCChamberHolderZ;

	G4double fVDCDistance1; // dist from magnet center to center of first VDC

	// Hodoscope
	G4double fPaddleHeight, fPaddleWidth, fPaddleThick;
	G4int fNumPaddles;
	G4double fHodoscopeWidth;
	G4double fPaddleWrapThick;
	G4double fHodoscopeDistance;
	G4double fHodoscopeOffset;

	G4LogicalVolume * fLogicHodoscopeHolder;
	G4VPhysicalVolume * fHodoscope[2];

	//misc
	DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
	MagneticField * fMagField;//pointer to the magnetic field

      
  private: //functions
    
	G4bool fCheckOverlaps;  

	void DefineMaterials();
	void ConstructYoke(G4ThreeVector center);
	void ConstructPole(G4ThreeVector center, G4bool bevelUp, G4int copy);
	void ConstructCoil(G4ThreeVector center, G4int copy);
	void ConstructBag(G4ThreeVector center);
	void ConstructChambers();
	void ConstructTarget();
	void ConstructHodoscope();

	G4VPhysicalVolume* ConstructAll();

 	//mag field stepper
	G4MagIntegratorStepper* fLocalStepper;
};

#endif

