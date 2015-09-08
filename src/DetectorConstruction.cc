//
//Default alignment: beam direction is x, detectors are to +y and z is up/down
//alter alignment by changing magnet rotation matrix
//axes: red: x, green: y, and blue: z

//******WARNING: use a pointer for G4RotationMatrix!  If you use &variable then you'll crash G4Tubs
//Notes:
	//rotateX/Y/Z means rotated ABOUT x/y/z axis
	//(extended) example F04 has non-uniform magnetic field

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4MagErrorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4EquationOfMotion.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Transform3D.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "WireChamberSD.hh"
#include "G4SDManager.hh"

//C++ includes
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

DetectorConstruction::DetectorConstruction()
:fMagField(0)
{
static const G4double inch = 2.54*cm;

//options:
fVerbose = false;

// Define materials first so we can use them in setting defaults
DefineMaterials();
  
//world
//sizes are based on the RayTrace diagram (z is conservative guess) *5!
fWorldSize = new G4ThreeVector(150.*cm*5, 200.*cm*5, 100.*cm*5);
fWorldCenter = new G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
fBackgroundMaterial = Helium; //
//fBackgroundMaterial = Vacuum; //
//
// Magnetic field
fBcenter = -0.4536*tesla;

//pole
fPoleOuterDiameter = 15.75*inch;
fPoleHeight = 4.75*inch; //total pole height
fPoleBevelHeight = 0.5*inch;
fPoleBevelledDiameter = fPoleOuterDiameter - 0.5*inch;
fPoleMaterial = Iron;

//coil
fCoilInnerDiameter = 15.75*inch;
fCoilOuterDiameter = 18.5*inch;
fCoilHeight = 3.625*inch; 
fCoilMaterial = Copper;
  
//yoke
fMagnetX = 30.75*inch;
fMagnetY = 22.75*inch;
fMagnetZ = 17.50*inch;
fYokeInnerX = 22.75*inch;
fYokeInnerZ = 11.50*inch;

fYokeSideX = (fMagnetX - fYokeInnerX)/2.;
fYokeSideY = fMagnetY;
fYokeSideZ = fYokeInnerZ;
fYokeTopX = fMagnetX;
fYokeTopY = fMagnetY;
fYokeTopZ = (fMagnetZ - fYokeInnerZ)/2.;
fYokeBeamHoleDiameter = 7.*cm;
fYokeMaterial = Iron;

// Bag
fBagThick = 0.1*mm;
fBagMaterial = Poly;
fMagnetToBagDistance = 1.0*mm;

//target
fTargetDistance = 42.26*cm;
fTargetWidth = 4.0*cm;
fTargetHeight = 4.0*cm;
fTargetToBagMaterial = Air;

// Target material and thickness can be changed through messenger
fTargetThick = 0.025*mm; // Thickness of Uranium target
fTargetMaterial = Uranium;  // Target material


//VDC Chambers
fVDCSizeX = 96*inch;
fVDCSizeY = 5.5*inch;
fVDCSizeZ = 36*inch;
fVDCSpacing = 15.5*inch;
fVDCSideWidth = 7.75*inch;
fVDCTopWidth = 7.5*inch;
fVDCSheetThick = 0.0005*inch;
fVDCLayerThick = 0.5*inch;
fVDCAlLayerThick = 0.75*inch;
fChamberSheetMaterial = Mylar;
fChamberGasMaterial = Ethane_Argon;

fVDCDistance1 = 85.0*cm; // distance of first VDC from magnet center

// Hodoscope
fPaddleHeight = 35.*cm;
fPaddleWidth = 5.*cm;
fPaddleThick = 1.*cm;
fNumPaddles = 29;
fHodoscopeWidth = 118.*cm;
fPaddleWrapThick = 1.0*mm;
fHodoscopeDistance = 185.*cm; // distance of Hodoscope centre from magnet centre
fHodoscopeOffset = 0.*cm; // offset of Hodoscope forward (+ = beam direction)


  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}



DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}



G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructAll();
}

void DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  

G4int ncomponents, natoms;
G4double fractionmass;
fNistMan = G4NistManager::Instance();
//
// define Elements
//
//from CRC handbook (online) http://www.hbcpnetbase.com/

G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.007975*g/mole);
//G4Element* He = new G4Element("Helium",symbol = "He", z = 2., a = 4.002602*g/mole);
G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.0106*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.006855*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 15.9994*g/mole);
//G4Element* Al = new G4Element("Aluminum",symbol="Al", z= 13., a= 26.9815386*g/mole);
//G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.0085*g/mole);
//G4Element* Fe = new G4Element("Iron",symbol = "Fe", z = 26., a = 55.845*g/mole);
//G4Element* Cu = new G4Element("Copper",symbol="Cu", z = 29., a = 63.546*g/mole);
//G4Element* U  = new G4Element("Uranium",symbol="U",z = 92., a = 238.02891*g/mole);
G4Element* Ar = new G4Element("Argon",symbol="Ar",z= 18.,a= 39.948*g/mole);



Sci = new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

Mylar = new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
Mylar->AddElement(C, natoms=10);
Mylar->AddElement(H, natoms= 8);
Mylar->AddElement(O, natoms= 4);

Poly = 
new G4Material("Polyethylene", density= 0.93*g/cm3, ncomponents=2);
Poly->AddElement(C, natoms=1);
Poly->AddElement(H, natoms=2);

Air = new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass = 0.7);
Air->AddElement(O, fractionmass = 0.3);

//50/50 mix partial pressures, assuming ideal gas this means they have 50/50 molar ratio
// moles of argon (Ar) = moles of Ethane (C2H6)
G4double ethane_density = 1.3562*mg/cm3; // 0 deg C (273 K) 1 atm
G4double ethane_molar_mass = 30.069*g/mole;
G4double argon_density = 1.784*mg/cm3; // 0 deg C (273 K) 1 atm
G4double argon_molar_mass = 39.95*g/mole;
G4double molar_mass = (ethane_molar_mass + argon_molar_mass)/2.; //ave of mix
// Convert to 25 deg C at 1 atm (298 K)
G4double temp = 273. + 25.;
G4double gas_density = 0.5*(argon_density*273./temp/argon_molar_mass
			+ ethane_density*273./temp/ethane_molar_mass)*molar_mass;
if(fVerbose) G4cout << "Ethane/Argon gas density = " << gas_density/(mg/cm3) << " mg/cm3" << G4endl;

Ethane_Argon = new G4Material("Ethane/Argon", gas_density, ncomponents=3);
Ethane_Argon->AddElement(Ar, natoms=1);
Ethane_Argon->AddElement(C, natoms=2);
Ethane_Argon->AddElement(H, natoms=6);


// Use G4-NIST materials data base
//

Helium = fNistMan->FindOrBuildMaterial("G4_He");
Aluminum = fNistMan->FindOrBuildMaterial("G4_Al");
Copper = fNistMan->FindOrBuildMaterial("G4_Cu");
Carbon = fNistMan->FindOrBuildMaterial("G4_C");
Uranium = fNistMan->FindOrBuildMaterial("G4_U");
Iron = fNistMan->FindOrBuildMaterial("G4_Fe");
Vacuum = fNistMan->FindOrBuildMaterial("G4_Galactic");


// print table
//
//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}


void DetectorConstruction::ConstructCoil(G4ThreeVector center, G4int copy)
//tubes of copper filled with water wrapped around magnetic pole
//I simplify by treating them as simply solid copper 
{

G4String name = "coil";
name += copy; //sprintf(name,"Coil%d",copy); //sprintf: printf except stores in a string, returns an integer equal to the # chars

fSolidCoil = new G4Tubs(	name,//name
				fCoilInnerDiameter/2., //min radius
				fCoilOuterDiameter/2., //max radius
				fCoilHeight/2., // z/2
				0*degree, //min phi
				360*degree); //size in phi

fLogicCoil = new G4LogicalVolume(	fSolidCoil,//solid volume
					fCoilMaterial,//material
					name);

fPhysCoil = new G4PVPlacement(	0,//rotation
				center,//origin
				fLogicCoil,//logical volume
				name,//name
				fLogicMagnet, //mother volume
				false,//bool operation?
				copy);//copy number

//visible attributes

  G4VisAttributes* coilVisAtt = new G4VisAttributes(G4Colour::Red());
  coilVisAtt->SetVisibility(true);
  //coilVisAtt->SetForceSolid(true); //must be solid
  coilVisAtt->SetForceWireframe(true); //must be wireframe
  fLogicCoil->SetVisAttributes(coilVisAtt);

if (fVerbose) 
	{
	printf("Just constructed Coil%d at x = %E cm, y = %E cm, z = %E cm.\n", copy, center.getX(), center.getY(), center.getZ()); 
	}
}

void DetectorConstruction::ConstructPole(G4ThreeVector center, G4bool bevelUp, G4int copy)
//center is center of BOTH bevel and base
{
//first, the unbevelled base

G4double zShift;
//*if bevel up: bottom edge is at -(baseHeight + bevelHeight)/2
	//therefore center of base is at: baseHeight/2 - (baseHeight + bevelHeight)/2 = bevelHeight/2. relative to center
if (bevelUp) zShift = -fPoleBevelHeight/2.;
//if bevel down then we swap sign: 
else zShift = fPoleBevelHeight/2.;
G4ThreeVector poleBaseCenter(center.getX(),center.getY(),center.getZ() + zShift);
G4String name = "poleBase";
fSolidPoleBase = 0; fLogicPoleBase = 0; fPhysPoleBase = 0;
name += copy;

	fSolidPoleBase = new G4Tubs(	name,
				0., //min radius
				fPoleOuterDiameter/2., //max radius
				(fPoleHeight - fPoleBevelHeight)/2.,//z/2
				0.*degree,//start phi
				360*degree);//size in phi
	fLogicPoleBase = new G4LogicalVolume(	fSolidPoleBase,
						fPoleMaterial, 
						name);
	fPhysPoleBase = new G4PVPlacement(	0,
						poleBaseCenter,
						fLogicPoleBase,
						name,
						fLogicMagnet,
						false,
						copy);
//now the bevelled top
//*if bevel UP: top edge is at (baseHeight + bevelHeight)/2
	//therefore center of bevel is at: (baseHeight + bevelHeight)/2 - bevelHeight/2 = fPoleHeight/2. - fPoleBevelHeight/2. relative to center
if (bevelUp) zShift = (fPoleHeight - fPoleBevelHeight)/2.;
//if bevel down then we swap sign:
else zShift = -(fPoleHeight - fPoleBevelHeight)/2.;
G4ThreeVector poleBevelCenter(center.getX(),center.getY(),center.getZ() + zShift);
name = "poleBevel";
name += copy;//sprintf(name,"PoleBevel%d",copy); 
//cerr << name << endl;
fSolidPoleBevel = 0; fLogicPoleBevel = 0; fPhysPoleBevel = 0;
	fSolidPoleBevel = new G4Cons(	name,
					0.,//r min1
					fPoleOuterDiameter/2.,//r max1
					0.,//r min2
					fPoleBevelledDiameter/2.,//r max2
					fPoleBevelHeight/2.,//z/2
					0.,//start phi
					360.*degree);//end phi
	fLogicPoleBevel = new G4LogicalVolume(	fSolidPoleBevel,
						fPoleMaterial,
						name);
	G4RotationMatrix * rotMatrix = new G4RotationMatrix();
	if (!bevelUp) rotMatrix->rotateX(180.*degree);
	fPhysPoleBevel = new G4PVPlacement(	rotMatrix, //flip upside down
						poleBevelCenter,
						fLogicPoleBevel,
						name,
						fLogicMagnet,
						false,
						copy);
//colouring (for both):

	G4VisAttributes* poleVisAtt = new G4VisAttributes(G4Colour::Gray());
	poleVisAtt->SetVisibility(true);
	//poleVisAtt->SetForceSolid(true); //must be solid
	poleVisAtt->SetForceWireframe(true); //must be wireframe
	fLogicPoleBevel->SetVisAttributes(poleVisAtt);	
	fLogicPoleBase->SetVisAttributes(poleVisAtt);	
				
	
}

void DetectorConstruction::ConstructYoke(G4ThreeVector center)
{
	//G4String name; //used to name volumes since apparently "" doesn't work right

	//bottom
	G4ThreeVector yokeBottomCenter = center;
	yokeBottomCenter.setZ(yokeBottomCenter.getZ() - fYokeSideZ/2. - fYokeTopZ/2.);

	fSolidYokeBottom = 0; fLogicYokeBottom = 0; fPhysYokeBottom = 0;
	fSolidYokeBottom = new G4Box(	"yokeBase",
					fYokeTopX/2.,
					fYokeTopY/2.,
					fYokeTopZ/2.);
	fLogicYokeBottom = new G4LogicalVolume( fSolidYokeBottom,
						fYokeMaterial, 
						"yokeBase");
	fPhysYokeBottom = new G4PVPlacement(	0,
						yokeBottomCenter,
						fLogicYokeBottom,
						"yokeBase",
						fLogicMagnet,
						false,
						0);
	//top
	G4ThreeVector yokeTopCenter = center;
	yokeTopCenter.setZ(yokeTopCenter.getZ() + fYokeSideZ/2. + fYokeTopZ/2.);

	// REP should be able to use just a copy of bottom
	/*
	fSolidYokeTop = new G4Box(	"yokeTop",
					fYokeTopX/2.,
					fYokeTopY/2.,
					fYokeTopZ/2.);
	fLogicYokeTop = new G4LogicalVolume(	fSolidYokeTop,
						fYokeMaterial, 
						"yokeTop");
	*/
	fPhysYokeTop = 0;
	fPhysYokeTop = new G4PVPlacement(	0,
						yokeTopCenter,
						fLogicYokeBottom,
						"yokeTop",
						fLogicMagnet,
						false,
						1);

	//upstream side
	G4ThreeVector yokeUpstreamSideCenter = center;
	yokeUpstreamSideCenter.setX(yokeUpstreamSideCenter.getX() - fYokeTopX/2. + fYokeSideX/2.);

	fSolidYokeUpstreamSide = 0;
	fSolidYokeUpstreamSide = new G4Box(	"yokeUpstreamSide",
						fYokeSideX/2.,
						fYokeSideY/2.,
						fYokeSideZ/2.);
	//hole (in upstream side)
	G4RotationMatrix * beamHoleRot = new G4RotationMatrix();
	beamHoleRot->rotateY(90.*degree);


	fSolidYokeBeamHole = 0;
	fSolidYokeBeamHole = new G4Tubs(	"yokeBeamHole",
						0., //min radius
						fYokeBeamHoleDiameter/2., //max radius
						fYokeSideX/2.*1.0001,//z/2
						0.*degree,
						360.*degree);

	fYokeWithHole = 0;
	fYokeWithHole = new G4SubtractionSolid(	"yokeWithBeamHole",//name
						fSolidYokeUpstreamSide,//G4VSolid pointer (base to cut from)
						fSolidYokeBeamHole,//G4VSolid pointer (hole)
						beamHoleRot, //rotation
						G4ThreeVector(0,0,0)); //in center of 1st solid
	//----------done hole--------------
	fLogicYokeUpstreamSide = 0;
	fLogicYokeUpstreamSide = new G4LogicalVolume(	fYokeWithHole, //solid vol
							fYokeMaterial,
							"yokeUpstreamSide"); //name
	fPhysYokeUpstreamSide = 0;
	fPhysYokeUpstreamSide = new G4PVPlacement(	0,
							yokeUpstreamSideCenter,
							fLogicYokeUpstreamSide,
							"yokeUpstreamSide",
							fLogicMagnet,
							false,
							0);
	//downstream side
	G4ThreeVector yokeDownstreamSideCenter = center;
	yokeDownstreamSideCenter.setX(yokeDownstreamSideCenter.getX() + fYokeTopX/2. - fYokeSideX/2.);
	// Downstream should just be a copy of upstream
	fPhysYokeDownstreamSide = 0;
	fPhysYokeDownstreamSide = new G4PVPlacement(	0,
							yokeDownstreamSideCenter,
							fLogicYokeUpstreamSide,
							"yokeDownstreamSide",
							fLogicMagnet,
							false,
							1);
	
	//colour it all up nice and fancy-like
	G4VisAttributes* yokeVisAtt = new G4VisAttributes(G4Colour::White());
	yokeVisAtt->SetVisibility(true);
	yokeVisAtt->SetForceWireframe(true); //must be wireframe
	//yokeVisAtt->SetForceSolid(true); //must be solid
	fLogicYokeUpstreamSide->SetVisAttributes(yokeVisAtt);	
	fLogicYokeBottom->SetVisAttributes(yokeVisAtt);	
	
}

G4VPhysicalVolume* DetectorConstruction::ConstructAll()
{

// Clean old geometry, if any
	G4cout << "Cleaning up ..." << G4endl;
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


     
fSolidWorld = 0; fLogicWorld = 0; fPhysWorld = 0;
// World

fSolidWorld = new G4Box(	"World",//its name
		                fWorldSize->getX()/2.,
				fWorldSize->getY()/2.,
				fWorldSize->getZ()/2.);	
                         
fLogicWorld = new G4LogicalVolume(fSolidWorld,		//its solid
                                   fBackgroundMaterial,	//its material
                                   "World");		//its name
                                   
fPhysWorld = new G4PVPlacement(		0,			//no rotation allowed
  					*fWorldCenter,		//center
                                 	fLogicWorld,		//its logical volume				 
                                 	"World",		//its name
                                 	0,			//its mother volume (None)
                                 	false,			//no boolean operation
                                 	0);			//copy number
// Visualization attributes
fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);

if (fVerbose) 
	{
	printf("Just constructed World at x = %E cm, y = %E cm, z = %E cm.\n",
	fWorldCenter->getX(), fWorldCenter->getY(), fWorldCenter->getZ()); 
	}

//magnet container: must be at least big enough to fit everything
fMagnetSize = new G4ThreeVector(fYokeTopX , fYokeTopY, fYokeSideZ + 2.*fYokeTopZ);
/*+ 2.*(fPipeVerticalLength + fPipeStartX  - fYokeTopX)*/
//add a bit to ensure everything fits (should be at least ~.1% to prevent f.p. errors)
// REP I don't this is necessary anymore - haven't followed this elsewhere.
//fMagnetSize->setX(fMagnetSize->getX()*1.01);
//fMagnetSize->setY(fMagnetSize->getY()*1.01);
//fMagnetSize->setZ(fMagnetSize->getZ()*1.01);


fSolidMagnet = 0; fLogicMagnet = 0;
//now build the magnet:
G4ThreeVector magnetCenter = *fWorldCenter; //magnet is at center
fSolidMagnet = new G4Box(	"Magnet",//its name
		                fMagnetSize->getX()/2.,
				fMagnetSize->getY()/2.,
				fMagnetSize->getZ()/2.);	
                         
fLogicMagnet = new G4LogicalVolume(fSolidMagnet,		//its solid
                                   fBackgroundMaterial,	//its material
                                   "Magnet");		//its name

//"yoke"
ConstructYoke(magnetCenter);

//poles

	//upper pole:
	G4ThreeVector poleCenter = magnetCenter;
	G4bool bevelup = false;
	poleCenter.setZ(poleCenter.getZ() + (fYokeSideZ - fPoleHeight)/2.);
	ConstructPole(poleCenter, bevelup, 0); //second is bevel up

	//lower pole:
	poleCenter = magnetCenter;
	bevelup = true;
	poleCenter.setZ(poleCenter.getZ() - (fYokeSideZ - fPoleHeight)/2.);
	ConstructPole(poleCenter, bevelup, 1); //second is bevel up

//coils

	//upper coil:
	G4ThreeVector coilCenter = magnetCenter;
	coilCenter.setZ(coilCenter.getZ() + (fYokeSideZ - fCoilHeight)/2. ); 
	ConstructCoil(coilCenter, 0);
                                          
	//lower coil:
	coilCenter = magnetCenter;
	coilCenter.setZ(coilCenter.getZ() - (fYokeSideZ - fCoilHeight)/2. ); 
	ConstructCoil(coilCenter, 1);

//now add the magnet:
fPhysMagnet = 0;
fPhysMagnet = new G4PVPlacement(
					//fMagnetRotation,
					0,			// no rotation
  					magnetCenter,		//center
                                 	fLogicMagnet,		//its logical volume				 
                                 	"Magnet",		//its name
                                 	fLogicWorld,		//its mother  volume
                                 	false,			//no boolean operation
                                 	0);	

//****colour magnet logical volume for debug purposes: 
G4VisAttributes* magVisAtt = new G4VisAttributes(G4Colour(1,0,0,0.9));//last number is opacity
magVisAtt->SetVisibility(true);
magVisAtt->SetForceSolid(true);
//fLogicMagnet->SetVisAttributes (magVisAtt);
fLogicMagnet->SetVisAttributes(G4VisAttributes::Invisible);

//chambers:
// Constructs a single (2 VDC) package
  ConstructChambers();

// Now place two copies of the chamber package into the world
  G4RotationMatrix ChamberRot;
  ChamberRot.rotateX(45.*deg); //to angle the chamber at 45 deg
  G4double PositionY = fVDCDistance1 + sqrt(2.)*(fVDCChamberHolderY/2. - fVDCSizeY/2. - fBagThick/2.);

  // Placement of the holder volume within the world volume
  // Beam left (copy number 0)
fPhysVDCpackage[0] = 0;
fPhysVDCpackage[0] =  new G4PVPlacement( G4Transform3D(ChamberRot, G4ThreeVector(0.,PositionY,0.)),
		    fLogicChamberHolder,
		    "Drift Chamber Package 1",
		    fLogicWorld, 
		    false,
		    0); 
		    //fCheckOverlaps);
  // Beam right (copy number 1)
  ChamberRot.set(0.,0.,0.); // back to identity
  ChamberRot.rotateZ(180.*deg);
  ChamberRot.rotateX(-45.*deg);
fPhysVDCpackage[1] = 0;
fPhysVDCpackage[1] =  new G4PVPlacement( G4Transform3D(ChamberRot, G4ThreeVector(0.,-PositionY,0.)),
		    fLogicChamberHolder,
		    "Drift Chamber Package 2",
		    fLogicWorld, 
		    false,
		    1); 

// Construct the Hodoscope.
	ConstructHodoscope();

// place two copies
// Beam left
  G4RotationMatrix HodoscopeRot(0.,0.,0.);
  HodoscopeRot.rotateX(-45.*deg);
  HodoscopeRot.rotateZ(180.*deg);
  fHodoscope[0] =  new G4PVPlacement(
	G4Transform3D(HodoscopeRot, G4ThreeVector(fHodoscopeOffset,fHodoscopeDistance,0.)),
	    fLogicHodoscopeHolder,
	    "Hodoscope 1",
	    fLogicWorld, 
	    false,
	    0); 

// Beam right
  HodoscopeRot.rotateZ(-180.*deg);
  fHodoscope[1] =  new G4PVPlacement(
	G4Transform3D(HodoscopeRot, G4ThreeVector(fHodoscopeOffset,-fHodoscopeDistance,0.)),
	    fLogicHodoscopeHolder,
	    "Hodoscope 2",
	    fLogicWorld, 
	    false,
	    1); 


//target
  ConstructTarget();

fPhysTarget = 0;
fPhysTarget =  new G4PVPlacement(0,
	    G4ThreeVector(-fTargetDistance, 0., 0.),
	    fLogicTargetContainer,
	    "Target",
	    fLogicWorld,
	    false,
	    0,
		    fCheckOverlaps);

//let's see if this works: whole world has magnetic field
SetMagField(fLogicWorld, fBcenter);

  //
  //always return the physical World
  //
  return fPhysWorld;
}

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4LogicalVolume * logicalMag, G4double Bcenter)
{
  if(fMagField) delete fMagField;		//delete the existing magn field

//  if(!getenv("MagneticFieldTable"))
//	{
//	G4Exception("SetMagField","Fatal", FatalException,
//	"Please setenv MagneticFieldTable to point to the magnetic field table.\n");
//	}
  //fMagField =  new G4UniformMagField(G4ThreeVector(0.,0.,magStrength));
//  G4String fieldTable = getenv("MagneticFieldTable");
//  fMagField = new MagField(fieldTable,0);  

//  MagneticField uses only a parameterization of the magnetic field.
  fMagField = new MagneticField();  
  fMagField->SetFieldB0(Bcenter);
  if(fVerbose) G4cout << "Magnetic Field B0 = "<< fMagField->GetFieldB0()/tesla << " Tesla." << G4endl;

  G4Mag_UsualEqRhs* fLocalEquation = new G4Mag_UsualEqRhs(fMagField);

  //steppers: http://geant4.web.cern.ch/geant4/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/Detector/electroMagneticField.html
/*
   fLocalStepper = new G4SimpleHeum( fEquation, nvar );         
            //  3rd  order, a good alternative to ClassicalRK
   fLocalStepper = new G4SimpleRunge( fEquation, nvar ); 
            //  2nd  order, for less smooth fields
   fLocalStepper = new G4CashKarpRKF45( fEquation );     
            // 4/5th order for very smooth fields 
*/

  G4int nvar = 8; //12 for spin tracking, 8 if no spin tracking
  fLocalStepper = new G4ClassicalRK4( fLocalEquation, nvar );
	//my field should be fairly smooth since it uses
	//the parameterization which is valid at all real points

  G4double fMinStep = 0.1*mm;
  
  G4ChordFinder* fLocalChordFinder = new G4ChordFinder(fMagField,
                                        fMinStep, fLocalStepper);
  
  G4double minEps= 1.0e-5;  //   Minimum & value for smallest steps
  G4double maxEps= 1.0e-4;  //   Maximum & value for largest steps

  G4FieldManager* fieldMgr = new G4FieldManager();
  fieldMgr->SetMinimumEpsilonStep(minEps);
  fieldMgr->SetMaximumEpsilonStep(maxEps);
  fieldMgr->SetDetectorField(fMagField);
  fieldMgr->SetChordFinder(fLocalChordFinder);

  G4double intersectionPrecision = 0.001*mm, delOneStep = 0.001*mm;
  fieldMgr->SetDeltaIntersection(intersectionPrecision);
  fieldMgr->SetDeltaOneStep(delOneStep);
  
  G4bool allLocal = true;
  logicalMag->SetFieldManager(fieldMgr,allLocal);

}

void DetectorConstruction::ConstructTarget()
{

// ================================================================================= 
// Target Construction
// Target package containes Air, target, polyethelene bag window
// Needs to be recalculated in case Target thickness has been changed.
	fAirTargetThick = fTargetDistance - fMagnetX/2.
		- fTargetThick/2. - fBagThick - fMagnetToBagDistance;
  if(fVerbose) G4cout << "AirTargetThick = " << fAirTargetThick/cm << " cm " << G4endl;

  fTargetContainerThick = fTargetThick + 2.*fAirTargetThick + 2.*fBagThick;
  G4Box* TargetContainerBox = new G4Box("Targ Cont. box", fTargetContainerThick/2., fTargetWidth, fTargetWidth);
  fLogicTargetContainer = 0;
  fLogicTargetContainer = new G4LogicalVolume(TargetContainerBox, fTargetToBagMaterial, "TargetContainer");
  fLogicTargetContainer->SetVisAttributes(G4VisAttributes::Invisible);

  G4Box* targetSolid = new G4Box("Target", fTargetThick/2., fTargetWidth/2., fTargetHeight/2.);

  fLogicTarget = 0;
  fLogicTarget = new G4LogicalVolume(targetSolid, fTargetMaterial, "Target");
  G4ThreeVector position(0., 0., 0.);
  new G4PVPlacement(0, 
			position,
			fLogicTarget,
			"Target", 
			fLogicTargetContainer, 
			false, 
			0); 
  G4Box* AirSolid = new G4Box("Target", fAirTargetThick/2., fTargetWidth/2., fTargetHeight/2.);
  fAirContainerLogic = 0;
  fAirContainerLogic = new G4LogicalVolume(AirSolid, fTargetToBagMaterial, "AirContainer");
  // we want to make this sensitive to recored e- e+ parameters just after target.
  // Declare Sensitive Detector
  G4String MonitorSDname = "Monitor";
  G4String MonitorHitsCollectionName = "MonitorHitCollection";
  if(fMonitorSD) delete fMonitorSD;
  fMonitorSD = new MonitorSD(MonitorSDname, MonitorHitsCollectionName);
  G4SDManager::GetSDMpointer()->AddNewDetector(fMonitorSD);
  fAirContainerLogic->SetSensitiveDetector(fMonitorSD);
  fAirContainerLogic->SetVisAttributes(G4VisAttributes::Invisible);

  position.setX(fTargetThick/2. + fAirTargetThick/2.);
  new G4PVPlacement(0, 
			position,
			fAirContainerLogic,
			"Target", 
			fLogicTargetContainer, 
			false, 
			0); 
  // Polyethelene windows
  G4Box* windowBox = new G4Box("Window Box", fBagThick/2., fTargetWidth, fTargetWidth);
  fLogicBagWall = 0;
  fLogicBagWall = new G4LogicalVolume(windowBox, fBagMaterial, "Bag Wall");
  position.setX(fTargetThick/2. + fAirTargetThick + fBagThick/2.);
  new G4PVPlacement(0, 
			position,
			fLogicBagWall,
			"Target", 
			fLogicTargetContainer, 
			false, 
			0); 

	

	G4VisAttributes* targetAtt = new G4VisAttributes(G4Colour::Yellow());
	targetAtt->SetVisibility(true);
	targetAtt->SetForceWireframe(true); //must be wireframe
	fLogicTarget->SetVisAttributes(targetAtt);	
	G4VisAttributes* bagAtt = new G4VisAttributes(G4Colour::Blue());
	bagAtt->SetVisibility(true);
	bagAtt->SetForceWireframe(true); //must be wireframe
	fLogicBagWall->SetVisAttributes(bagAtt);	

}

void DetectorConstruction::ConstructChambers()
{
//all dimensions from JP Leckey's PhD thesis
/*Layers:	
	Al frame (0.75",  0.5" away from foil)
	Gas (0.5")
	HV Foil (0.5" away from wire, 0.0005" thick)
	Wire (U)
	Spacer (0.5")
	HV Foil (0.5" away from wire, 0.0005" thick)
	Wire (V)
	HV Foil (0.5" away from wire, 0.0005" thick)
	Gas frame (0.5", for letting gas in/out)
	Al frame (0.75", 0.5" away from foil)
Notes:	-Gas frame has same dimensions as Al frame except for width and # holes (which I'll ignore for the sim)
*/
  // =================================================================================
  // Chamber Construction

  // Chamber package Holder for the two VDCs

  fVDCChamberHolderX = fVDCSizeX;
  fVDCChamberHolderY = fVDCSpacing + fVDCSizeY + 2.*fBagThick;
  fVDCChamberHolderZ = fVDCSpacing + fVDCSizeZ;

  G4Box* chamberHolderBox = 
    new G4Box("Chamber Holder", fVDCChamberHolderX/2., fVDCChamberHolderY/2., fVDCChamberHolderZ/2.);

  fLogicChamberHolder = 0;
  fLogicChamberHolder = 
    new G4LogicalVolume(chamberHolderBox,
 			fBackgroundMaterial, 
			"2 VDC Chamber Holder");
  // make package holder invisible
  fLogicChamberHolder->SetVisAttributes(G4VisAttributes::Invisible);

  // A single VDC
  G4Box* VDCBox = 
    new G4Box("VDC Box", fVDCSizeX/2., fVDCSizeY/2. + fBagThick, fVDCSizeZ/2.);
  fLogicVDCHolder = 0;
  fLogicVDCHolder = new G4LogicalVolume(VDCBox,
 			Air, 
			"VDC Holder");
  // =================================================================================
  //aluminum frame for a single VDC:

  G4double	al_layer_X = fVDCSideWidth;
  G4double	al_layer_Y = fVDCSizeY;
  G4double	al_layer_Z = fVDCSizeZ - 2.*fVDCTopWidth;

  //left:
  G4Box* alFrameLeft =    new G4Box("Al Frame Left",  //its name
        al_layer_X/2., al_layer_Y/2., al_layer_Z/2.); //its size
  
  G4LogicalVolume* logicAlFrameLeft =                         
    new G4LogicalVolume(alFrameLeft,            //its solid
                        Aluminum,
                        "Al Frame Left/Right");         //its name    
  G4ThreeVector position = G4ThreeVector(-fVDCSizeX/2. + fVDCSideWidth/2.,0.,0.);
  new G4PVPlacement(0, 
			position,
			logicAlFrameLeft,
			"AlFrameLeft", 
			fLogicVDCHolder, 
			false, 
			0); 
			//fCheckOverlaps);

  //right:
  position.setX( position.getX() + fVDCSizeX - fVDCSideWidth);
  new G4PVPlacement(	0, 
			position,
			logicAlFrameLeft,
			"Al Frame Right", 
			fLogicVDCHolder, 
			false, 
			1); 
			//fCheckOverlaps);

  //top:
	al_layer_X = fVDCSizeX;
	al_layer_Y = fVDCSizeY;
	al_layer_Z = fVDCTopWidth;
  G4Box* alFrameTop =   new G4Box("Al Frame Top",                    //its name
		        al_layer_X/2., al_layer_Y/2., al_layer_Z/2.); //its size
  
  G4LogicalVolume* logicAlFrameTop = new G4LogicalVolume(
			alFrameTop,            //its solid
                        Aluminum,
                        "Al Frame Top/Bottom");         //its name    

  position.setX(0.);
  position.setY(0.);
  position.setZ(fVDCSizeZ/2. - fVDCTopWidth/2.);
  new G4PVPlacement(	0, 
			position,
			logicAlFrameTop,
			"Al Frame Top", 
			fLogicVDCHolder, 
			false, 
			0); 
			//fCheckOverlaps);

  //bottom:
  position.setZ( position.getZ() - fVDCSizeZ + fVDCTopWidth);
  new G4PVPlacement(	0, 
			position,
			logicAlFrameTop,
			"Al Frame Bottom", 
			fLogicVDCHolder, 
			false, 
			1); 
			//fCheckOverlaps);

	G4VisAttributes* VDCFrameVisAtt = new G4VisAttributes(G4Colour::Yellow());
	VDCFrameVisAtt->SetVisibility(true);
	//VDCFrameVisAtt->SetForceWireframe(true); //must be wireframe
	VDCFrameVisAtt->SetForceSolid(true); //must be Solid
	logicAlFrameTop->SetVisAttributes(VDCFrameVisAtt);	
	logicAlFrameLeft->SetVisAttributes(VDCFrameVisAtt);	
  // ===================================================================================
  
  // Declare Sensitive Detector
  G4String wireChamberSDname = "WireChamber";
  G4String HitsCollectionName = "HitCollection";
  if(fchamberSD) delete fchamberSD;
  fchamberSD = new WireChamberSD(wireChamberSDname, HitsCollectionName);
  G4SDManager::GetSDMpointer()->AddNewDetector( fchamberSD );

  // The wire chamber layers = 
  // Individual Wire Segments
  G4double sheetX, sheetY, sheetZ;
  G4double gas1X, gas1Y, gas1Z;
  G4double gas2X, gas2Y, gas2Z;
  sheetX = fVDCSizeX - 2.*fVDCSideWidth; sheetY = fVDCSheetThick; sheetZ = fVDCSizeZ - 2.*fVDCTopWidth;
  gas1X = sheetX; gas1Y = fVDCLayerThick - fVDCSheetThick; gas1Z = sheetZ;
  gas2X = sheetX; gas2Y = 2.*fVDCLayerThick - fVDCSheetThick; gas2Z = sheetZ;

  G4Box * sheetBox = new G4Box("Sheet Box", sheetX/2., sheetY/2., sheetZ/2.);
  G4Box * gas1Box = new G4Box ("Gas layer 1 box", gas1X/2., gas1Y/2., gas1Z/2.); 
  G4Box * gas2Box = new G4Box ("Gas layer 2 box", gas2X/2., gas2Y/2., gas2Z/2.); 

   fLogicVDCsheet = 0;
   fLogicVDCgas1 = 0;
   fLogicVDCgasSD[0] = 0;
   fLogicVDCgasSD[1] = 0;

   fLogicVDCsheet = new G4LogicalVolume( sheetBox, fChamberSheetMaterial, "Sheet");
   fLogicVDCgas1 = new G4LogicalVolume( gas1Box, fChamberGasMaterial, "Gas 1");
   fLogicVDCgasSD[0] = new G4LogicalVolume( gas2Box, fChamberGasMaterial, "Gas SD 0");
   fLogicVDCgasSD[1] = new G4LogicalVolume( gas2Box, fChamberGasMaterial, "Gas SD 1");
   fLogicVDCgasSD[0]->SetSensitiveDetector(fchamberSD); 
   fLogicVDCgasSD[1]->SetSensitiveDetector(fchamberSD); 

	G4VisAttributes* VDCLayerVisAtt = new G4VisAttributes(G4Colour::Grey());
	VDCLayerVisAtt->SetVisibility(true);
	VDCLayerVisAtt->SetForceWireframe(true); //must be wireframe
	G4VisAttributes* VDCGasVisAtt = new G4VisAttributes(G4Colour::Blue());
	VDCGasVisAtt->SetVisibility(true);
	VDCGasVisAtt->SetForceWireframe(true); //must be wireframe
	fLogicVDCsheet->SetVisAttributes(VDCLayerVisAtt);	
	fLogicVDCgas1->SetVisAttributes(VDCGasVisAtt);	
	fLogicVDCgasSD[0]->SetVisAttributes(VDCGasVisAtt);
	fLogicVDCgasSD[1]->SetVisAttributes(VDCGasVisAtt);
	//fLogicVDCgas1->SetVisAttributes(G4VisAttributes::Invisible);
	//fLogicVDCgasSD[0]->SetVisAttributes(G4VisAttributes::Invisible);
	//fLogicVDCgasSD[1]->SetVisAttributes(G4VisAttributes::Invisible);
 
  // now layer them up
  position.setX(0.);
  position.setY(fVDCSizeY/2. - fVDCAlLayerThick - sheetY/2.);
  position.setZ(0.);

  // sheet 0
  //G4VPhysicalVolume * PhysSheet0 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 0",
						fLogicVDCHolder, false, 0);
  // gas 1 layer 0
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * Physgas0 =
  new G4PVPlacement(0, position, fLogicVDCgas1, "VDC Gas 0",
						fLogicVDCHolder, false, 0);
  // sheet 1
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * PhysSheet1 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 1",
						fLogicVDCHolder, false, 1);
  // gas 1 layer 1
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * Physgas1 =
  new G4PVPlacement(0, position, fLogicVDCgas1, "VDC Gas 1",
						fLogicVDCHolder, false, 1);
  // sheet 2
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * PhysSheet2 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 2",
						fLogicVDCHolder, false, 2);
  // gas 2 layer 0 SD
  position.setY(position.getY() - (sheetY/2. + gas2Y/2.) );
  //		   fPhysVDCGasSD[0] =
  new G4PVPlacement(0, position, fLogicVDCgasSD[0], "VDC Gas SD 0",
						fLogicVDCHolder, false, 0);
  // sheet 3
  position.setY(position.getY() - (sheetY/2. + gas2Y/2.) );
  //G4VPhysicalVolume * PhysSheet3 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 3",
						fLogicVDCHolder, false, 3);
  // gas 2 layer 1 SD
  position.setY(position.getY() - (sheetY/2. + gas2Y/2.) );
  //		   fPhysVDCGasSD[1] =
  new G4PVPlacement(0, position, fLogicVDCgasSD[1], "VDC Gas SD 1",
						fLogicVDCHolder, false, 1);
  // sheet 4
  position.setY(position.getY() - (sheetY/2. + gas2Y/2.) );
  //G4VPhysicalVolume * PhysSheet4 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 4",
						fLogicVDCHolder, false, 4);
  // gas 1 layer 2
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * Physgas2 =
  new G4PVPlacement(0, position, fLogicVDCgas1, "VDC Gas 2",
						fLogicVDCHolder, false, 2);
  // sheet 5
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * PhysSheet5 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 5",
						fLogicVDCHolder, false, 5);
  // gas 1 layer 3
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * Physgas3 =
  new G4PVPlacement(0, position, fLogicVDCgas1, "VDC Gas 3",
						fLogicVDCHolder, false, 3);
  // sheet 6
  position.setY(position.getY() - (sheetY/2. + gas1Y/2.) );
  //G4VPhysicalVolume * PhysSheet6 =
  new G4PVPlacement(0, position, fLogicVDCsheet, "VDC Sheet 6",
						fLogicVDCHolder, false, 6);
  // Now add helium bag windows outside the Al Frame
  sheetX = fVDCSizeX; sheetY = fBagThick; sheetZ = fVDCSizeZ;
  G4Box * BagBox = new G4Box("Bag Box", sheetX/2., sheetY/2., sheetZ/2.);
   fLogicVDCBag = 0;
   fLogicVDCBag = new G4LogicalVolume( BagBox, fBagMaterial, "VDC Bag");
  // two copies - one on each side
  position.setX(0.);
  position.setY(fVDCSizeY/2. + fBagThick/2.);
  position.setZ(0.);
  new G4PVPlacement(0, position, fLogicVDCBag, "VDC Bag 0",
						fLogicVDCHolder, false, 0);
  position.setY(-(fVDCSizeY/2. + fBagThick/2.));
  new G4PVPlacement(0, position, fLogicVDCBag, "VDC Bag 1",
						fLogicVDCHolder, false, 1);
	fLogicVDCBag->SetVisAttributes(VDCLayerVisAtt);	
  // ======================================================
  // One VDC is now constructed.
  // Place two copies in the Chamber holder

  position.setX(0.);
  position.setY(+fVDCChamberHolderY/2.-fVDCSizeY/2.-fBagThick/2.);
  position.setZ(-fVDCChamberHolderZ/2.+fVDCSizeZ/2.+fBagThick/2.);

  new G4PVPlacement(	0,  // no rotation
			position,
			fLogicVDCHolder,
			"Rear VDC", 
			fLogicChamberHolder, 
			false, 
			1);  // Copy number
			//fCheckOverlaps);


  position.setY(-fVDCChamberHolderY/2.+fVDCSizeY/2.+fBagThick/2.);
  position.setZ(+fVDCChamberHolderZ/2.-fVDCSizeZ/2.-fBagThick/2.);

  new G4PVPlacement(	0,  // no rotation
			position,
			fLogicVDCHolder,
			"Front VDC", 
			fLogicChamberHolder, 
			false, 
			0);  // Copy number
			//fCheckOverlaps);
// we are done ConstructChambers()
}


void
DetectorConstruction::ConstructHodoscope()
{
  // Logical volume to hold one hodoscope
  G4double PaddleX = fPaddleWidth + fPaddleWrapThick*2.;
  G4double PaddleY = fPaddleThick + fPaddleWrapThick*2.;
  G4double PaddleZ = fPaddleHeight;
  G4double PaddleSpacing = (fHodoscopeWidth - PaddleX )/(fNumPaddles - 1);

  G4double HodHolderX = fHodoscopeWidth;
  G4double HodHolderY = 2.*PaddleY;
  G4double HodHolderZ = PaddleZ;

  G4Box* hodoscopeHolderBox = 
    new G4Box("Hodoscope Holder box", HodHolderX/2., HodHolderY/2., HodHolderZ/2.);
  fLogicHodoscopeHolder = 0;
  fLogicHodoscopeHolder = 
    new G4LogicalVolume(hodoscopeHolderBox, Air, "Hodoscope Holder");
  // make holder invisible
  fLogicHodoscopeHolder->SetVisAttributes(G4VisAttributes::Invisible);

  // create a single paddle Russian doll style
  // first the wrapping
  G4Box* paddleBox = new G4Box("Paddle Box", PaddleX/2., PaddleY/2., PaddleZ/2.);
  G4LogicalVolume * paddleLogic = new G4LogicalVolume(paddleBox, Poly, "Paddle Wrap");
  // then the scintillator
  G4Box* scintBox = new G4Box("Paddle Scint", fPaddleWidth/2., fPaddleThick/2., fPaddleHeight/2.);
  G4LogicalVolume * scintLogic = new G4LogicalVolume(scintBox, Sci, "Paddle Scintillator");
  G4ThreeVector position(0.,0.,0.);
   new G4PVPlacement(0, position, scintLogic, "PaddleScint", paddleLogic, false, 0);  

  // Declare Sensitive Detector
  G4String PaddleSDname = "Paddle";
  G4String PaddleHitsCollectionName = "PaddleHitCollection";
  if(fPaddleSD) delete fPaddleSD;
  fPaddleSD = new PaddleSD(PaddleSDname, PaddleHitsCollectionName);
  G4SDManager::GetSDMpointer()->AddNewDetector(fPaddleSD);
  scintLogic->SetSensitiveDetector(fPaddleSD);
	// visual attributes
	G4VisAttributes* paddleWrap_att = new G4VisAttributes(G4Colour::Blue());
	paddleWrap_att->SetVisibility(true);
	paddleWrap_att->SetForceSolid(true); //must be Solid
	G4VisAttributes* paddleScint_att = new G4VisAttributes(G4Colour::Yellow());
	paddleScint_att->SetVisibility(true);
	paddleScint_att->SetForceSolid(true); //must be Solid
	paddleLogic->SetVisAttributes(paddleWrap_att);
	scintLogic->SetVisAttributes(paddleScint_att);

  // now put the paddle into the holder
  position.setX(-HodHolderX/2. - PaddleSpacing/2.);
  position.setY( PaddleY/2.);
  position.setZ( 0.); 
  char names[20];
  for(G4int i = 0; i < fNumPaddles; i++)
		{
		position.setX( position.getX() + PaddleSpacing);
		if(i%2 == 0) position.setY( PaddleY/2.);
		else position.setY( -PaddleY/2.);
		sprintf(names,"Paddle %d",i);
		G4String name(names);
		new G4PVPlacement(0, position, paddleLogic, name, fLogicHodoscopeHolder, false, i);  
		}

}





#if 0
  G4LogicalVolume* logicWireSeg; 
  G4VisAttributes* wireSegAttrib = new G4VisAttributes(G4Colour(1.,1.,0));
  G4VisAttributes* wireSegAttrib2 = new G4VisAttributes(G4Colour(1.,0.,0)); //different colour so it's easier to see that the two are different

  // Declare Sensitive Detector
  G4String wireChamberSDname = "WireChamberSDname";
  WireChamberSD* chamberSD = new WireChamberSD(wireChamberSDname,
                                            "TrackerHitsCollection");
  G4double angle = 90.*deg - 26.45*deg; //Leckey's thesis
  G4double xLoc, zLoc;
  G4int index = 1;
  xLoc = wireSpacing; //temp

  //wire layers: I combine layers 1 (u) and 2 (v) into one big layer with twice the width
  //note: in xz plane
/*  G4RotationMatrix *rot = new G4RotationMatrix;
  rot->rotateX(M_PI/2.*rad);



  G4Box * box = new G4Box	("box",
				chamb_innerX/2.,
				chamb_innerZ/2.,
				1.*inch
				);


  G4Box * wireChamb = new G4Box("wire", 
				wireSpacing/2.,
				sqrt(chamb_innerX/2.*chamb_innerX/2. + chamb_innerZ/2.*chamb_innerZ/2.), //make sure it is always big enough
				1.*inch);



	G4RotationMatrix * rotMatrix1 = new G4RotationMatrix;
	rotMatrix1->rotateZ(angle);

	G4RotationMatrix * rotMatrix2 = new G4RotationMatrix;
	rotMatrix2->rotateZ(angle);
	rotMatrix2->rotateZ(90.*deg);

	G4IntersectionSolid *solidWireSeg1; //one intersection
	G4IntersectionSolid *solidWireSeg2; //both intersections
	
  G4double shift = chamb_innerZ/2./tan(90.*deg - angle) - wireSpacing/2.; //this is to make sure the entire plane is covered: from trig
  G4double xLoc2; //for second dimension
  G4double xLength = chamb_innerZ/tan(90.*deg - angle); //length along x direction of first intersection solid
  xLoc = -chamb_innerX/2. - shift; 
  while (xLoc < chamb_innerX/2. + shift + wireSpacing/4.)
	{	
	solidWireSeg1 = new G4IntersectionSolid("Wire  " + index, box, wireChamb, rotMatrix1, G4ThreeVector(xLoc,0.,0.));
	xLoc2 = xLoc - xLength/2.; 
	while (xLoc2 < xLoc + xLength/2.)//chamb_innerX/2. + shift + wireSpacing/4.)
		{
		solidWireSeg2 = new G4IntersectionSolid("Wire Box " + index, solidWireSeg1, wireChamb, rotMatrix2, G4ThreeVector(xLoc,0.,0.));



		logicWireSeg =                         
		    new G4LogicalVolume(solidWireSeg2,            //its solid
                        fChamberMaterial,             //its material 
                        "Wire Box "+ index);         //its name

		// Declare Sensitive Detector
		G4SDManager::GetSDMpointer()->AddNewDetector( chamberSD );
		logicWireSeg->SetSensitiveDetector(chamberSD);

	
		new G4PVPlacement(rot, G4ThreeVector(0,0,0), logicWireSeg,"Wire Box " + index,
				fLogicChamberHolder, false, index, fCheckOverlaps);
		xLoc2 += wireSpacing;
		index++;
		logicWireSeg->SetVisAttributes(wireSegAttrib);
		//wireSegAttrib->SetForceSolid(true); //must be solid
		}
	xLoc += wireSpacing;
	}
*/
 //faster
 //I opt to rotate the entire coordinate system
  G4PVPlacement * pPVPlacement;
  angle = 90.*deg - 26.45*deg;
  G4RotationMatrix *rot = new G4RotationMatrix;
  rot->rotateX(M_PI/2.*rad);



  G4Box * box = new G4Box	("Wire Layer",
				chamb_innerX/2.,
				chamb_innerZ/2.,
				1.*inch //each wire layer is 0.5" long, so the half length of the two is 1"
				);


  G4Box * wireChamb = new G4Box("Wire Box", 
				wireSpacing/2.,
				wireSpacing/2.,
				1.*inch);



	G4RotationMatrix * rotMatrix = new G4RotationMatrix;
	rotMatrix->rotateZ(angle);

	G4IntersectionSolid *solidWireSeg; //one intersection

  G4double shiftX = chamb_innerX + wireSpacing; //this is to make sure the entire plane is covered: worst case scenario
  G4double shiftZ = chamb_innerZ + wireSpacing; //this is to make sure the entire plane is covered: worst case scenario

  G4double xLocP, zLocP; //rotated coordinates, P for prime
  xLoc = -chamb_innerX/2. - shiftX;
  ofstream outFile; //debug
  outFile.open ("WirePositions.txt");
  G4ThreeVector threeVec;//debug
  while (xLoc < chamb_innerX/2. + shiftX)
	{	
	zLoc = -chamb_innerZ/2. - shiftZ;
	while (zLoc < chamb_innerZ/2. + shiftZ)
		{
		xLocP = xLoc*cos(90.*deg - angle) - zLoc*sin(90.*deg - angle);
		zLocP = xLoc*sin(90.*deg - angle) + zLoc*cos(90.*deg - angle);

		//check if the boolean solid would actually exist (i.e. volume > 0):
		if (abs(zLocP) > chamb_innerZ/2. + wireSpacing/2.*sqrt(2.)*1.01) //z is too large/small
			{
			zLoc += wireSpacing;
			continue;
			}
		else if (abs(xLocP) > chamb_innerX/2. + wireSpacing/2.*sqrt(2.)*1.01) //x is too large/small, x' depends on z so we can't just break this loop
			{
			zLoc += wireSpacing;
			continue;
			}
		threeVec.setX(xLocP);
		threeVec.setY(zLocP);
		threeVec.setZ(0.);
		outFile << index << "\t" << threeVec.getX() << "\t" << threeVec.getY() << "\t" << threeVec.getZ() << "\n";
	


		solidWireSeg = new G4IntersectionSolid("Wire Box " + index, box, wireChamb, G4Transform3D(rotMatrix->inverse(), G4ThreeVector(xLocP,zLocP,0.))); //g4transform3d is so that it will store the data instead of just using the pointer (which I may change)

		//if (solidWireSeg->GetCubicVolume() < wireSpacing*wireSpacing*1.0*inch/10000.) continue; //sanity check to make sure volumes actually overlap
		//if (solidWireSeg == 0) continue;



		logicWireSeg =                         
		    new G4LogicalVolume(solidWireSeg,            //its solid
                        fChamberMaterial,             //its material
                        "Wire Box "+ index);         //its name

		// Declare Sensitive Detector
		G4SDManager::GetSDMpointer()->AddNewDetector( chamberSD );
		logicWireSeg->SetSensitiveDetector(chamberSD);

	
		pPVPlacement = new G4PVPlacement(rot, G4ThreeVector(0,0,0), logicWireSeg,"Wire Box " + index,
				fLogicChamberHolder, false, index, fCheckOverlaps);
		//pPVPlacement->SetTranslation(G4ThreeVector(xLocP,zLocP,0.));
		zLoc += wireSpacing;
		index++;
		logicWireSeg->SetVisAttributes(wireSegAttrib);
		//wireSegAttrib->SetForceSolid(true); //must be solid
		}
	xLoc += wireSpacing;
	}
  outFile.close();//debug

}
#endif


#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructAll());
}

void DetectorConstruction::Print()
	{
	G4cout << "----------------------------------" << G4endl;
	G4cout << "  Changeable Detector Parameters  " << G4endl;
	G4cout << "  Target Thickness =  " << fTargetThick/mm << " mm" << G4endl;
	G4cout << "  Target Material  =  " << fTargetMaterial->GetName() << G4endl;
	G4cout << "  Magnetic Field   =  " << fBcenter/tesla << " Tesla" << G4endl;
	G4cout << "----------------------------------" << G4endl;
	}

//Prints the information about each volume: Taken from example A01
void DetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth)
{
	//spacing to make things look nice
  	for(G4int isp=0;isp<depth;isp++) { G4cout << "  "; }
  
	G4cout << aVolume->GetName() << "[";
	
	//if there is a number, print it.
	G4int copyNo = aVolume->GetCopyNo();
	if(copyNo >= 0) {
		G4cout << copyNo;}
	else if( aVolume->IsParameterised() ){
		G4cout << "P";}
	else if( aVolume->IsReplicated() ){
		G4cout << "R";}

  	G4cout << "] " << aVolume->GetLogicalVolume()->GetName() << " ";
        G4cout << aVolume->GetLogicalVolume()->GetNoDaughters() << " ";
        G4cout << aVolume->GetLogicalVolume()->GetMaterial()->GetName();
	
	//if a volume is a sensitive detector
  	if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
  	{
    		G4cout << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()->GetFullPathName();
  	}
  
	G4cout << G4endl;

	//recursivly call this function for all volumes
  	for(G4int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++) {
  		DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),depth+1); }
		
} //End of DumpeGeometricalTree()

