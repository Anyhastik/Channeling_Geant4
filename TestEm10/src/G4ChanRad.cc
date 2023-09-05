//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// History:
//2022.12 - created for TR simulation by ***
//

#include "G4ChanRad.hh"

// #include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4EmProcessSubType.hh"
//#include "G4Integrator.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4RotationMatrix.hh"
// #include "G4SandiaTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Timer.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
//#include "G4VSolid.hh"
#include "G4PhysicsModelCatalog.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>*/

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// #include <boost/math/tools/assert.hpp>
//#include <boost/math/special_functions/next.hpp> // for float_distance





////////////////////////////////////////////
//#define BOOST_NO_EXCEPTIONS
//////////////////////////////////////////////////

// #include <boost/multiprecision/cpp_dec_float.hpp>

// using boost::multiprecision::number;
// using boost::multiprecision::cpp_dec_float;

// // Re-compute using 5 extra decimal digits precision (22) than double (17).
// #define MP_DIGITS10 unsigned (std::numeric_limits<double>::max_digits10 + 5)

// typedef cpp_dec_float<MP_DIGITS10> mp_backend;
// typedef number<mp_backend> mp_type;

///////////////////////////////////////////////////////////////////////////////////////////////////////////



#include "phys_const.hh"
#include "crystal_functions.hh"
#include "globals.hh"


using namespace std;


struct radiation
{
    vector<ld1> w;
    vector<ld1> rad;
};
////////////////////////////////////////////////////////////////////////////
// Constructor, destructor
G4ChanRad::G4ChanRad(G4LogicalVolume* anEnvelope,
                                   G4Material* foilMat, G4Material* gasMat,
                                   G4double a, G4double b, G4int n, G4double EGamma,
                                   const G4String& processName,
                                   G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
  , fGammaCutInKineticEnergy(nullptr)
  , fAngleDistrTable(nullptr)
  , fEnergyDistrTable(nullptr)
  , fAngleForEnergyTable(nullptr)
  , fPlatePhotoAbsCof(nullptr)
  , fGasPhotoAbsCof(nullptr)
  , fGammaTkinCut(0.0)
{
  verboseLevel = 1;
  secID = G4PhysicsModelCatalog::GetModelID("model_XTRenergyLoss"); //TO CHANGE
  SetProcessSubType(fTransitionRadiation);							//TO CHANGE

  //понаписанное мной
	//int nbeam = 5;
	// num = 5000 * nbeam;
	// std::unique_ptr<d[]> x_path(new d[(G4int)num]);
    // std::unique_ptr<d[]> time_moving(new d[(G4int)num]);
    // std::unique_ptr<d[]> x_velocity(new d[(G4int)num]);

    // std::unique_ptr<d[]> z_path(new d[(G4int)num]);
    // std::unique_ptr<d[]> z_velocity(new d[(G4int)num]);
	// poi = G4int (num);
	// d x_path[poi];
	// d time_moving[poi];
	// d x_velocity[poi];
	// d z_path[poi];
	// d z_velocity[poi];
	//defenition of the crystal
	x_path.resize(poi);
	time_moving.resize(poi);
	x_velocity.resize(poi);
	z_path.resize(poi);
	z_velocity.resize(poi);
	
	G4String crystaltype = "si"; 
    G4int kk = 1, ll = 1, nn = 1;
    G4int sigen = 1;
	// ld1 dp;
	// array<ld1, 2*nmax+1> cpot;
    tuple<ld1, array<ld1, 2*nmax+1>> tup = GPHYS_CONST::cpotss(sigen, kk, ll, nn, crystaltype);
    tie(dp, cpot) = tup;
	crystaltickness = 20; // 10-20mkm
    //ld1 t0 = crystaltickness * pow(10, -6) / c0;
    //d n = 5000 * nbeam;
    // G4int n1 = 5000;
    // ld1 h = (t0 - a) / n1;
	//функция,которая тут вызвается и выдает траектории
	cout<<"tetta="<<tetta<<endl;
	GetTrajectories(x_path, time_moving, x_velocity, z_path, z_velocity,
					crystaltype, sigen, num, kk, ll, nn, cpot, crystaltickness);
		
	
	//периоды
	// min_TT;
    // array<ld1, nbeam> TT;
   // tuple<ld1, array<ld1, nbeam>> tup1 = 
	GetPeriods(sigen, kk, ll, nn, crystaltype);
	//tie(min_TT, TT) = tup1;
	
	
	
  //
  fPtrGamma    = nullptr;
  fMinEnergyTR = fMaxEnergyTR = fMaxThetaTR = fGamma = fEnergy = 0.0;
  fVarAngle = fLambda = fTotalDist = fPlateThick = fGasThick = 0.0;
  fAlphaPlate = 100.;
  fAlphaGas = 40.;
  fMinEnergyTR = 0.01 *  hbar * 4 * pi * pow(gama, 2) / min_TT;
  fMaxEnergyTR =  hbar * 4 * pi * pow(gama, 2) / min_TT;
  G4cout << "fMaxEnergyTR = " << fMaxEnergyTR << G4endl;



  fTheMinEnergyTR = CLHEP::keV * 1.; //  1.; // 
  fTheMaxEnergyTR = CLHEP::keV * 100.; // 40.; //

  fTheMinAngle    = 1.e-8;  //
  fTheMaxAngle    = 4.e-4;
  tables_count = 0;
  fEgamma = EGamma;
  fTotBin = 2;  //  number of bins in linear scale
  fBinTR  =  100; //   number of bins in TR vectors 

  // min/max angle2 in log-vectors

  fMinThetaTR = 3.0e-9; 
  fMaxThetaTR = 1.0e-4;

  
  // Proton energy vector initialization
  fProtonEnergyVector =												//когда-нибудь поменять на электроны!!!
    new G4PhysicsLogVector(fMinProtonTkin, fMaxProtonTkin, fTotBin);
// linear wektor we need above
  fXTREnergyVector =
    new G4PhysicsLogVector(fTheMinEnergyTR, fTheMaxEnergyTR, fBinTR);

  fEnvelope = anEnvelope;

  fPlateNumber = n;
  if(verboseLevel > 0)
    G4cout << "### G4ChanRad: the number of TR radiator plates = "
           << fPlateNumber << G4endl;
  if(fPlateNumber == 0)
  {
    G4Exception("G4ChanRad::G4ChanRad()", "VXTRELoss01",
                FatalException, "No plates in X-ray TR radiator");
  }
  // default is XTR dEdx, not flux after radiator
  // fExitFlux      = false; 							//for Post step doiing
  // default angle distribution according numerical G4integration
  fFastAngle     = false; // no angle according sum of delta-functions by default
  // fAngleRadDistr = false;
  // fCompton       = false;

  fLambda = DBL_MAX;

  // Mean thicknesses of plates and gas gaps
  fPlateThick = a;
  fGasThick   = b;
  fTotalDist  = fPlateNumber * (fPlateThick + fGasThick);
  if(verboseLevel > 0)
    G4cout << "total radiator thickness = " << fTotalDist / cm << " cm"
           << G4endl;

  // // index of plate material
  // fMatIndex1 = foilMat->GetIndex();
  // if(verboseLevel > 0)
    // G4cout << "plate material = " << foilMat->GetName() << G4endl;

  // // index of gas material
  // fMatIndex2 = gasMat->GetIndex();
  // if(verboseLevel > 0)
    // G4cout << "gas material = " << gasMat->GetName() << G4endl;

  // // plasma energy squared for plate material
  // fSigma1 = fPlasmaCof * foilMat->GetElectronDensity();
  // if(verboseLevel > 0)
    // G4cout << "plate plasma energy = " << std::sqrt(fSigma1) / eV << " eV"
           // << G4endl;

  // // plasma energy squared for gas material
  // fSigma2 = fPlasmaCof * gasMat->GetElectronDensity();
  // if(verboseLevel > 0)
    // G4cout << "gas plasma energy = " << std::sqrt(fSigma2) / eV << " eV"
           // << G4endl;

  // // Compute cofs for preparation of linear photo absorption
  // ComputePlatePhotoAbsCof();
  // ComputeGasPhotoAbsCof();

  pParticleChange = &fParticleChange;
}



// ⣿⣿⣿⣿⣿⢙⠛⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠛⢋⠉⠘⣿⣿
// ⣿⣿⣿⣿⡏⠄⠄⠄⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⡿⠟⠉⠄⠄⠄⠄⠁⢸⣿
// ⣿⣿⣿⣿⠄⠄⠄⠄⠄⠄⠈⠙⢿⣿⣿⣿⡿⠋⠄⠄⠄⠄⠄⢀⠄⠄⢹⣿
// ⣿⣿⣿⣿⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⢀⠄⠄⢸⣿
// ⣿⣿⣿⠿⢤⣤⡔⠒⣶⣶⣶⣶⣶⡖⠒⠒⠒⠒⣶⣖⠒⣶⣦⣤⣤⣄⣸⣿
// ⣿⣷⡀⠄⠈⣿⣿⡄⢻⣿⣿⣿⣿⣿⡄⠄⠄⠄⢹⣿⡆⢹⣿⣿⣿⣿⣿⣿
// ⣿⣿⣷⠄⠄⠸⣿⣷⡈⣿⣿⡿⠹⣿⣿⡀⠄⠄⠄⢻⣿⡄⢻⣿⣿⣿⣿⣿
// ⣿⣿⣿⣷⣾⣷⣿⣿⣿⣿⡿⠁⠄⢻⣿⣿⣿⣿⣿⣿⣿⣷⣼⣿⣿⣿⣿⣿
// ⣿⣿⣿⣿⠄⣿⣿⣿⣿⡿⠁⠄⣤⠈⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⣿
// ⣿⣿⣿⠟⠄⠄⠉⠉⠉⠄⠄⡿⠛⢻⠈⠙⠟⠛⠛⠻⠿⠿⢿⠿⠿⠁⢀⣿
// ⣿⣿⣿⣿⡝⠐⠄⠄⠄⠄⠄⣀⠄⠄⠄⠄⠩⠐⡀⠄⠄⠈⠉⢿⣿⣿⣿⣿
// ⣿⣿⣿⣿⣿⡖⠄⣀⠄⠄⠛⠻⣿⣿⠟⠛⠄⠄⠄⠄⠄⣀⣴⣿⣿⣿⣿⣿
// ⣿⣿⣿⣿⣿⠇⣠⠊⠈⠐⠤⢤⣀⣀⣠⠤⣀⣀⣀⠛⠉⠄⠉⠻⣿⣿⣿⣿
// ⣿⣿⣿⡿⢋⡰⠁⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠈⢿⣿⣿
// ⣿⣿⡋⠐⢰⠁⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠘⣿⣿

///////////////////////////////////////////////////////////////////////////
G4ChanRad::~G4ChanRad()
{
  delete fProtonEnergyVector;
  delete fXTREnergyVector;
  if(fEnergyDistrTable)
  {
    fEnergyDistrTable->clearAndDestroy();
    delete fEnergyDistrTable;
  } 
  // if(fAngleRadDistr)
  // {
    // fAngleDistrTable->clearAndDestroy();
    // delete fAngleDistrTable;
  // }
  // if(fAngleForEnergyTable)
  // {
    // fAngleForEnergyTable->clearAndDestroy();
    // delete fAngleForEnergyTable;
  // }
}


void G4ChanRad::ProcessDescription(std::ostream& out) const
{
  out << "Base class for 'fast' parameterisation model describing X-ray "
         "Channeling\n"
         "radiation.\n";
}

///////////////////////////////////////////////////////////////////////////////
// Returns condition for application of the model depending on particle type
G4bool G4ChanRad::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.0);
}

/////////////////////////////////////////////////////////////////////////////////
// Calculate step size for XTR process inside raaditor
G4double G4ChanRad::GetMeanFreePath(const G4Track& aTrack, G4double,
                                           G4ForceCondition* condition)
{
  G4int iTkin, iPlace;
  G4double lambda, sigma, kinEnergy, mass, gamma;
  G4double charge, chargeSq, massRatio, TkinScaled;
  G4double E1, E2, W, W1, W2;

  *condition = NotForced;

  if(aTrack.GetVolume()->GetLogicalVolume() != fEnvelope)
    lambda = DBL_MAX;
  else
  {
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    kinEnergy                          = aParticle->GetKineticEnergy();
    mass  = aParticle->GetDefinition()->GetPDGMass();
    gamma = 1.0 + kinEnergy / mass; //fEgamma
    if(verboseLevel > 1)
    {
      G4cout << " gamma = " << gamma << ";   fGamma = " << fGamma << G4endl;
    }

    if(std::fabs(gamma - fGamma) < 0.05 * gamma)
      lambda = fLambda;
    else
    {
		G4cout << " we are here " << G4endl;
		
	  charge     = aParticle->GetDefinition()->GetPDGCharge();
      chargeSq   = charge * charge;
      massRatio  = proton_mass_c2 / mass;
      TkinScaled = kinEnergy * massRatio;

      for(iTkin = 0; iTkin < fTotBin; ++iTkin)
      {
       		G4cout << " we are here 1 " << G4endl;
	   if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))
          break;
      }
      iPlace = iTkin - 1;

		G4cout << " we are here 2 " << G4endl;

      if(iTkin == 0)
        lambda = DBL_MAX;  // Tkin is too small, neglect of TR photon generation
	
      else  // general case: Tkin between two vectors of the material
      {
		  
		 G4cout << " we are here 3 " << G4endl;
		 
		 
        if(iTkin == fTotBin)
        {
			
		   G4cout << " we are here 4" << G4endl;
		   
          sigma = (*(*fEnergyDistrTable)(iPlace))(0) * chargeSq;
        }
        else
        {
		
		  G4cout << " we are here 5 " << G4endl;
			
          E1    = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1);
          E2    = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
          W     = 1.0 / (E2 - E1);
          W1    = (E2 - TkinScaled) * W;
          W2    = (TkinScaled - E1) * W;
          sigma = ((*(*fEnergyDistrTable)(iPlace))(0) * W1 +
                   (*(*fEnergyDistrTable)(iPlace + 1))(0) * W2) *
                  chargeSq;
        }
        if(sigma < DBL_MIN)
          lambda = DBL_MAX;
        else
          lambda = 1. / sigma;
        fLambda = lambda;
        fGamma  = gamma;
        // if(verboseLevel > 0)
        // {
          G4cout << " lambda = " << lambda / mm << " mm" << G4endl;
        // }
      }
    }
  }
  return lambda;
}

//////////////////////////////////////////////////////////////////////////
// G4interface for build table from physics list
void G4ChanRad::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  if(pd.GetPDGCharge() == 0.)
  {
    G4Exception("G4ChanRad::BuildPhysicsTable", "Notification",
                JustWarning, "XTR initialisation for neutral particle ?!");
  }
  BuildEnergyTable();

  // if(fAngleRadDistr)
  // {
    // if(verboseLevel > 0)
    // {
      // G4cout
        // << "Build angle for energy distribution according the current radiator"
        // << G4endl;
    // }
    // BuildAngleForEnergyBank();
  // }
}

//my thing

ld1 G4ChanRad::integrator_trap(ld1 y2, ld1 y1, ld1 w2, ld1 w1)
{
	return (y2+y1) * 0.5 * fabs(w1-w2);
}


//end of my thing

//////////////////////////////////////////////////////////////////////////
// Build G4integral energy distribution of XTR photons
void G4ChanRad::BuildEnergyTable()
{
  G4int iTkin, iTR, iPlace;
  G4double radiatorCof = 1.0;  // for tuning of XTR yield
  G4double energySum   = 0.0;

  fEnergyDistrTable = new G4PhysicsTable(fTotBin);
  // if(fAngleRadDistr)
    // fAngleDistrTable = new G4PhysicsTable(fTotBin);

  fGammaTkinCut = 0.0;

  // setting of min/max TR energies
  // if(fGammaTkinCut > fTheMinEnergyTR)
    // fMinEnergyTR = fGammaTkinCut;
  // else
    // fMinEnergyTR = fTheMinEnergyTR;

  // if(fGammaTkinCut > fTheMaxEnergyTR)
    // fMaxEnergyTR = 2.0 * fGammaTkinCut;
  // else
    // fMaxEnergyTR = fTheMaxEnergyTR;

  // G4Integrator<G4ChanRad, G4double (G4ChanRad::*)(G4double)>
    // G4integral;

  G4cout.precision(4);
  G4Timer timer;
  timer.Start();

  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "Lorentz Factor"
           << "\t"
           << "TR photon number" << G4endl;
    G4cout << G4endl;
  }
  //downstairs
   
    tables_count += 1;
	if (tables_count == 1)
	{
		x_path_copy = x_path;
		time_moving_copy = time_moving;
		x_velocity_copy= x_velocity;
		z_path_copy = z_path;
		z_velocity_copy = z_velocity;
	}
	if (tables_count == 2)
	{
		x_path = x_path_copy;
		time_moving = time_moving_copy;
		x_velocity = x_velocity_copy;
		z_path = z_path_copy;
		z_velocity = z_velocity_copy;
	}
	ld1 beta_z = z_velocity.at(0)/ (c0 * pow(10, 10));
	G4int file_counter = 1;
	vector<ld1> times;
    times.clear();
    vector<ld1> velocities;
    velocities.clear();
	
	ld1 t1[nbeam];
	
	for(G4int i = 0; i < nbeam; i++)
        t1[i] = crystaltickness * pow(10, -6) / (beta_z * c0);
	
	G4int Nt[nbeam];
    for(G4int i = 0; i < nbeam; i++)
        Nt[i] = G4int(t1[i] / TT[i]); 
	
	
	 G4cout<< "here0" << G4endl;
	 
	ld1 max_w = 4 * pi * pow(gama, 2) / min_TT;
    ld1 rad_counter_w = 0.01*max_w;
    radiation emission;
    emission.w.clear();
    emission.rad.clear();
    ofstream fout89;
    fout89.open("radiation100.txt");
    char ch1, ch2;

	ld1 interr[nbeam];
	//downstairs
	print(TT);
	for (int j = 0; j < 100; j++)
    {
        emission.w.push_back(rad_counter_w * hbar);
        rad_counter_w += 0.01 * max_w;
    }
	for(G4int i = 0; i < nbeam; i++)
    { G4cout<< "here1" << G4endl;
	    while (*time_moving.begin() <= TT[i])
        {
            times.push_back(*time_moving.begin());
			time_moving.erase(time_moving.begin());
            velocities.push_back(*x_velocity.begin());
			x_velocity.erase(x_velocity.begin());
            file_counter++;
        }
        
        times.push_back(*time_moving.begin());
		time_moving.erase(time_moving.begin());
        velocities.push_back(*x_velocity.begin());
		x_velocity.erase(x_velocity.begin());

        for (G4int j = 0; j < 5; j++) //число гармоник
        {
            interr[j] = trapezoid(times, velocities, j + 1, TT[i]);
			
        }
		rad_counter_w = 0.01 * max_w;
        for(int l=0; l<100; l++)
        {
            emission.rad.push_back(dWdE(Nt[i], rad_counter_w, TT[i], interr)/(nbeam * rad_counter_w * hbar * 10));
            //cout<<"back="<<emission.rad.back()<<endl;
            fout89 <<  rad_counter_w*hbar << " " << emission.rad.back() << "\n";
            rad_counter_w += 0.01 * max_w;
        }
        // while (rad_counter_w <= max_w)
        // {
			// //G4cout<< "here star" << G4endl;
            // emission.w.push_back(rad_counter_w * hbar);
			// ld1 en_en = dWdE(Nt[i], rad_counter_w, TT[i], interr);
            // emission.rad.push_back(en_en/(nbeam * rad_counter_w * hbar * 10));
            // fout89 << rad_counter_w * hbar << " " << /*fixed << setprecision(6) <<*/ emission.rad.back() << "\n";
            // rad_counter_w += 0.01 * max_w;
			// //G4cout<< "here end" << G4endl;
        // }
        fout89 << "-" << " " << "-" << "\n";
        times.clear();
        velocities.clear();
        while (file_counter <= (n_steps * (i+1) ) && i != nbeam-1)
        {
            //ifs >> time >> velocity;
			time_moving.erase(time_moving.begin());
			x_velocity.erase(x_velocity.begin());
            file_counter++;
        }
        //ifs >> time >> velocity;
		time_moving.erase(time_moving.begin());
		x_velocity.erase(x_velocity.begin());
		
        file_counter++;
        //rad_counter_w = 0.01 *max_w;
        //emission.w.clear();
    }
	G4cout<< "here2" << G4endl;
	fout89.close();
	//выше расчитывается для каждой частицы излучение
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // for(iTkin = 0; iTkin < fTotBin; ++iTkin)
  // {}	  // Lorentz factor loop
   // energyVector->PutValue(fBinTR - 1, energySum);//why -1?
   
   
   //reverse(emission.rad.begin(), emission.rad.end());
   std::ofstream total("total_N_integrated_for_Geant4.txt");
   
  for(iTkin = 0; iTkin < fTotBin;++iTkin)//why do we need this loop? 0,1...
  {
    G4PhysicsLinearVector* energyVector =
      new G4PhysicsLinearVector(fMinEnergyTR, fMaxEnergyTR, fBinTR);   //vector for energy distribution

    fGamma = fEgamma;
	
    energySum = 0.0;
	energyVector->PutValue(fBinTR - 1, energySum);
	total << energyVector->GetLowEdgeEnergy(fBinTR - 1) << " " << energySum/nbeam << "\n";
	for(iTR = fBinTR - 2; iTR >= 0; iTR--)
	{
		//if (energyVector->GetLowEdgeEnergy(iTR) !=0)
		//{
		//Number = SpectralXTRdEdx(energyVector.at(iTR),&Nt, &t1, emission.rad, min_TT);
		ld1 spectra1 = SpectralXTRdEdx(iTR, energyVector->GetLowEdgeEnergy(iTR), emission.rad, min_TT);
		ld1 spectra2 = SpectralXTRdEdx(iTR+1, energyVector->GetLowEdgeEnergy(iTR+1), emission.rad, min_TT);
		energySum += integrator_trap(
				spectra1,
				spectra2,
				energyVector->GetLowEdgeEnergy(iTR),
				energyVector->GetLowEdgeEnergy(iTR+1));			
		//ld1 y2, ld1 y1, ld1 w2, ld1 w1
		energyVector->PutValue(iTR, energySum);//fTotalDist
		total << energyVector->GetLowEdgeEnergy(iTR) << " " << energySum << "\n";
		//}
		//else energyVector->PutValue(iTR, energySum/fTotalDist);
	}
	// cout<<"energy vec:"<<endl;
	// for (G4int i=0; i<fBinTR; i++)
		// // cout<<energyVector[i];
	// cout<<energyVector->GetLowEdgeEnergy(i)<< " ";
	// cout<<"energy vec end."<<endl;
	//iTkin = 0;
	//total << energyVector->GetLowEdgeEnergy(iTR) << " " << energySum << "\n";
    iPlace = iTkin;
    fEnergyDistrTable->insertAt(iPlace, energyVector);
	
	//fEnergyDistrTable->insertAt(iPlace+1, energyVector);

    if(verboseLevel > 0)
    {
      G4cout << fGamma << "\t" << energySum << " disttable0 = "<< (*(*fEnergyDistrTable)(iPlace))(0) << G4endl;
    }
  }
  total.close();
  cout<<"energy vec:"<<endl;
	for (G4int i=0; i<fBinTR; i++)
		// cout<<energyVector[i];
	cout<< "hw = " << (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(i) << "   " << (*(*fEnergyDistrTable)(iPlace))(i)<< " ";
	cout<<"energy vec end."<<endl;
  timer.Stop();
  G4cout.precision(6);
  if(verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << "total time for build X-ray TR energy = "
           << timer.GetUserElapsed() << " s" << G4endl;
  }
 fGamma = 0.;
  return; //возвращает полное число фотонов
}

//////////////////////////////////////////////////////////////////////////////
// The main function which is responsible for the treatment of a particle
// passage through G4Envelope with discrete generation of G4Gamma
G4VParticleChange* G4ChanRad::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step& aStep)
{
  G4int iTkin;
  G4double energyTR, theta, theta2, phi, dirX, dirY, dirZ;

  fParticleChange.Initialize(aTrack);

  if(verboseLevel > 1)
  {
    G4cout << "Start of G4ChanRad::PostStepDoIt " << G4endl;
    G4cout << "name of current material =  "
           << aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()
           << G4endl;
  }
  if(aTrack.GetVolume()->GetLogicalVolume() != fEnvelope)
  {
    if(verboseLevel > 0)
    {
      G4cout << "Go out from G4ChanRad::PostStepDoIt: wrong volume "
             << G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  else
  {
    G4StepPoint* pPostStepPoint        = aStep.GetPostStepPoint();
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    // Now we are ready to Generate one TR photon
    G4double kinEnergy = aParticle->GetKineticEnergy();
    G4double mass      = aParticle->GetDefinition()->GetPDGMass();
    G4double gamma     = fEgamma; //1.0 + kinEnergy / mass;

    if(verboseLevel > 1)
    {
      G4cout << "gamma = " << gamma << G4endl;
    }
    G4double massRatio           = proton_mass_c2 / mass;
    G4double TkinScaled          = kinEnergy * massRatio;
    G4ThreeVector position       = pPostStepPoint->GetPosition();
    G4ParticleMomentum direction = aParticle->GetMomentumDirection();
    G4double startTime           = pPostStepPoint->GetGlobalTime();

    for(iTkin = 0; iTkin < fTotBin; ++iTkin)
    {
      if(TkinScaled < fProtonEnergyVector->GetLowEdgeEnergy(iTkin))
        break;
    }

    if(iTkin == 0)  // Tkin is too small, neglect of TR photon generation
    {
      if(verboseLevel > 0)
      {
        G4cout << "Go out from G4ChanRad::PostStepDoIt:iTkin = " << iTkin
               << G4endl;
      }
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    else  // general case: Tkin between two vectors of the material
    {
      fParticleChange.SetNumberOfSecondaries(1);

      energyTR = GetXTRrandomEnergy(TkinScaled, iTkin);

      if(verboseLevel > 1)
      {
        G4cout << "energyTR = " << energyTR / keV << " keV" << G4endl;
      }
      if(fAngleRadDistr)
      {
        // theta2 = GetRandomAngle(energyTR, iTkin);
        // if(theta2 > 0.)
          // theta = std::sqrt(theta2);
        // else
          // theta = 0.;
      }
      else
        theta = std::fabs(G4RandGauss::shoot(0.0, pi / gamma));

      if(theta >= 0.1)
        theta = 0.1;

      phi = twopi * G4UniformRand();

      dirX = std::sin(theta) * std::cos(phi);
      dirY = std::sin(theta) * std::sin(phi);
      dirZ = std::cos(theta);

      G4ThreeVector directionTR(dirX, dirY, dirZ);
      directionTR.rotateUz(direction);
      directionTR.unit();

      G4DynamicParticle* aPhotonTR =
        new G4DynamicParticle(G4Gamma::Gamma(), directionTR, energyTR);

      // A XTR photon is set on the particle track inside the radiator
      // and is moved to the G4Envelope surface for standard X-ray TR models
      // only. The case of fExitFlux=true

      // if(fExitFlux)
      // {
        // const G4RotationMatrix* rotM =
          // pPostStepPoG4int->GetTouchable()->GetRotation();
        // G4ThreeVector transl = pPostStepPoG4int->GetTouchable()->GetTranslation();
        // G4AffineTransform transform = G4AffineTransform(rotM, transl);
        // transform.Invert();
        // G4ThreeVector localP = transform.TransformPoG4int(position);
        // G4ThreeVector localV = transform.TransformAxis(directionTR);

        // G4double distance =
          // fEnvelope->GetSolid()->DistanceToOut(localP, localV);
        // if(verboseLevel > 1)
        // {
          // G4cout << "distance to exit = " << distance / mm << " mm" << G4endl;
        // }
        // position += distance * directionTR;
        // startTime += distance / c_light;
      // }
      G4Track* aSecondaryTrack = new G4Track(aPhotonTR, startTime, position);
      aSecondaryTrack->SetTouchableHandle(
        aStep.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID(aTrack.GetTrackID());

      fParticleChange.AddSecondary(aSecondaryTrack);
      fParticleChange.ProposeEnergy(kinEnergy);
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

/////////////////////////////////////////////////////////////////////////
// For second G4integration over energy
G4double G4ChanRad::SpectralXTRdEdx(G4int i, G4double energy, vector<ld1> emission, ld1 min_TT)     //нужная функция
{
  //G4int i;
  static constexpr G4int iMax = 8;
  G4double angleSum           = 0.0;
 //rad_counter_w = 0;
    ld1 sum_dwde = 0;
    ld1 coeffi = 4*pi*pow(gama, 2) / min_TT;
    //vector<ld1> radi;
   // radi.clear();
    // std::ofstream total("total_rad_new.txt");
    // for (G4int i = 0; i < 100; i++)
    // {
       // emission.w.push_back(rad_counter_w * hbar);
        //rad_counter_w += 0.01 * coeffi;
        for (G4int j = i; j < nbeam*(100); j+=101)
        {
            sum_dwde +=  emission.at(j);
        }
        //radi.push_back(sum_dwde);
        //total << emission.w.back() << " " << /*fixed << setprecision(6) <<*/ radi.back()/nbeam << "\n";
        //sum_dwde = 0;
    // }
   
    //fout89.close();
    //total.close();
	G4double result = sum_dwde;
  return result;
  // G4double lim[iMax] = { 0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 };

  // for(i = 0; i < iMax; ++i)
    // lim[i] *= fMaxThetaTR;

  // G4Integrator<G4ChanRad, G4double (G4ChanRad::*)(G4double)>
    // G4integral;

  // fEnergy = energy;
  // {
    // for(i = 0; i < iMax - 1; ++i)
    // {
      // angleSum += G4integral.Legendre96(
        // this, &G4ChanRad::SpectralAngleXTRdEdx, lim[i], lim[i + 1]);
    // }
  // }
  //return angleSum;
}

//////////////////////////////////////////////////////////////////////////////
// Check number of photons for a range of Lorentz factors from both energy
// and angular tables
void G4ChanRad::GetNumberOfPhotons()
{
  G4int iTkin;
  G4double gamma, numberE;

  std::ofstream outEn("numberE.dat", std::ios::out);
  outEn.setf(std::ios::scientific, std::ios::floatfield);

  std::ofstream outAng("numberAng.dat", std::ios::out);
  outAng.setf(std::ios::scientific, std::ios::floatfield);

  for(iTkin = 0; iTkin < fTotBin; ++iTkin)  // Lorentz factor loop
  {
    gamma =
      1.0 + (fProtonEnergyVector->GetLowEdgeEnergy(iTkin) / proton_mass_c2);
    numberE = (*(*fEnergyDistrTable)(iTkin))(0);
    if(verboseLevel > 1)
      G4cout << gamma << "\t\t" << numberE << "\t" << G4endl;
    if(verboseLevel > 0)
      outEn << gamma << "\t\t" << numberE << G4endl;
  }
  return;
}

/////////////////////////////////////////////////////////////////////////
// Returns random energy of a X-ray TR photon for given scaled kinetic energy
// of a charged particle
G4double G4ChanRad::GetXTRrandomEnergy(G4double scaledTkin, G4int iTkin)
{
  G4int iTransfer, iPlace;
  G4double transfer = 0.0, position, E1, E2, W1, W2, W;

  iPlace = iTkin - 1;

  if(iTkin == fTotBin)  // relativistic plato, try from left
  {
    position = (*(*fEnergyDistrTable)(iPlace))(0) * G4UniformRand();

    for(iTransfer = 0;; ++iTransfer)
    {
      if(position >= (*(*fEnergyDistrTable)(iPlace))(iTransfer))
        break;
    }
    transfer = GetXTRenergy(iPlace, position, iTransfer);
  }
  else
  {
    E1 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin - 1);
    E2 = fProtonEnergyVector->GetLowEdgeEnergy(iTkin);
    W  = 1.0 / (E2 - E1);
    W1 = (E2 - scaledTkin) * W;
    W2 = (scaledTkin - E1) * W;

    position = ((*(*fEnergyDistrTable)(iPlace))(0) * W1 +
                (*(*fEnergyDistrTable)(iPlace + 1))(0) * W2) *
               G4UniformRand();

    for(iTransfer = 0;; ++iTransfer)
    {
      if(position >= ((*(*fEnergyDistrTable)(iPlace))(iTransfer) *W1 +
                      (*(*fEnergyDistrTable)(iPlace + 1))(iTransfer) *W2))
        break;
    }
    transfer = GetXTRenergy(iPlace, position, iTransfer);
  }
  if(transfer < 0.0)
    transfer = 0.0;
  return transfer;
}

////////////////////////////////////////////////////////////////////////
// Returns approximate position of X-ray photon energy during random sampling
// over G4integral energy distribution
G4double G4ChanRad::GetXTRenergy(G4int iPlace, G4double, G4int iTransfer)
{
  G4double x1, x2, y1, y2, result;

  if(iTransfer == 0)
  {
    result = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);
  }
  else
  {
    y1 = (*(*fEnergyDistrTable)(iPlace))(iTransfer - 1);
    y2 = (*(*fEnergyDistrTable)(iPlace))(iTransfer);

    x1 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer - 1);
    x2 = (*fEnergyDistrTable)(iPlace)->GetLowEdgeEnergy(iTransfer);

    if(x1 == x2)
      result = x2;
    else
    {
      if(y1 == y2)
        result = x1 + (x2 - x1) * G4UniformRand();
      else
      {
        result = x1 + (x2 - x1) * G4UniformRand();
      }
    }
  }
  return result;
}
//дописано мной ниже

void G4ChanRad::GetTrajectories(vector<d> &x_path, vector<d> &time_moving, 
								vector<d> &x_velocity,
								vector<d> &z_path, vector<d> &z_velocity,
								string crystaltype, int sigen, 
								d num, int k, int l, int m, 
								array<ld1, 2*nmax+1> cpot,
								ld1 crystaltickness)
{
	//вообще константы надо подключить через файл с константами мой
	ld1 k1, k2, k3, k4, k1z, k2z, k3z, k4z, k0, k00;
    ld1 B, Bz;
	cout<<"START"<<endl;
	ofstream fout, foutx2, foutz1, foutz2;
    fout.open("table_x_path.txt");
    foutx2.open("table_x_velocity.txt");
    //foutx3.open("table_d2x.txt");
    foutz1.open("table_z_path.txt");
    foutz2.open("table_z_velocity.txt");
    //foutz3.open("table_d2z.txt");
	
	ld1 a = 0;
    ld1 t0 = crystaltickness * pow(10, -6) / c0;
    int n1 = n_steps;
    ld1 h = (t0 - a) / ld1(n1);
	num = G4int(num);
	
	//initial values
    time_moving.at(0) = 0;
    x_path.at(0)  = 1 / (ld1(nbeam) + 1) * 1 * dp / 2; // 1 here means ~i
    x_velocity.at(0)  = c0 * pow(10, 10) * sqrt(1 - 1 / pow(gama, 2)) * sin(tetta);
   
    z_path.at(0)  = 0;
    z_velocity.at(0)  = c0 * pow(10, 10) * sqrt(1 - 1 / pow(gama, 2)) * cos(tetta);
    //loop
    int s = 0;
     for (G4int i = 1; i <=nbeam; i++)
    {
        x_path.at(s) = 1 / (ld1(nbeam) + ld1(1)) * ld1(i) * dp / ld1(2); 
        x_velocity.at(s)  = c0 * pow(10, 10) * sqrt(1 - 1 / pow(gama, 2)) * sin(tetta);
		z_path.at(s)  = 0;
		z_velocity.at(s)  = c0 * pow(10, 10) * sqrt(1 - 1 / pow(gama, 2)) * cos(tetta);
        s++;
        for (G4int j = 1; j < n1; j++)
        {
            time_moving.at(s) = a + j * h;
            k0 = ducr(time_moving.at(s-1), x_path.at(s-1), cpot, dp, crystaltype); // it’s the second derivative from х, x'' = functiti()
            k1 = h * k0;

            k2 = h * ducr(time_moving.at(s-1) + h * 0.5, x_path.at(s-1) + x_velocity.at(s-1) * h * 0.5 
				+ k1 / 8 * h, cpot, dp, crystaltype);
            k3 = h * ducr(time_moving.at(s-1) + h, x_path.at(s-1) + x_velocity.at(s-1) * h 
				+ k2 * h * 0.5, cpot, dp, crystaltype);
            B = x_velocity.at(s-1) + (k1 + 2 * k2) / 6;
            //foutx3 << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << k0 << "\n"; //file
            x_path.at(s) = x_path.at(s-1) + h * B;
            x_velocity.at(s) = B + (2 * k2 + k3) / 6;
         
            k00 = functz(time_moving.at(s-1), z_path.at(s-1), z_velocity.at(s-1), x_velocity.at(s-1), k0);
            k1z = h * k00;

            //foutz3 << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << k00 << "\n";  //file
            k2z = h * functz(time_moving.at(s-1) + h * 0.5, z_path.at(s-1) + z_velocity.at(s-1) * h * 0.5 + k1z / 8 * h,
                            z_velocity.at(s-1) + k1z * 0.5, x_velocity.at(s-1), k0);

            k3z = h * functz(time_moving.at(s-1) + h * 0.5, z_path.at(s-1) + z_velocity.at(s-1) * h * 0.5 + k1z / 8 * h,
                            z_velocity.at(s-1) + k2z * 0.5, x_velocity.at(s-1), k0);

            k4z = h * functz(time_moving.at(s-1) + h, z_path.at(s-1) + z_velocity.at(s-1) * h + k3z / 2 * h,
                            z_velocity.at(s-1) + k3z, x_velocity.at(s-1), k0);

            Bz = z_velocity.at(s-1) + (k1z + k2z + k3z) / 6;
            z_path.at(s) = z_path.at(s-1) + h * Bz;
            z_velocity.at(s) = Bz + (k2z + k3z + k4z) / 6;
      
            s++;

        }
    cout<<"particle = " << i<<endl;  
    }
	cout<<"DONE"<<endl;
	
	for (int j = 0; j < s; j++)
    {
           fout << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << x_path[j] << "\n";
            foutx2 << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << x_velocity[j] << "\n";
            foutz1 << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << z_path[j] << "\n";
            foutz2 << fixed << setprecision(20) << time_moving[j] << " " << fixed << setprecision(10) << z_velocity[j] << "\n";  
    }
	
    fout.close();
    foutx2.close();
    //foutx3.close();
    //foutime_moving.close();
    foutz2.close();
    foutz1.close();
}
//tuple<ld1, array<ld1, nbeam>> 
void G4ChanRad::GetPeriods(G4int sigen,G4int k0, G4int l0, G4int n0, G4String crystaltype)
{
	// ld1 dp;
    // array<ld1, 2*nmax+1> cpot;
    // tuple<ld1, array<ld1, 2*nmax+1>> tupt = cpotss(sigen, k0, l0, n0, crystaltype);
    // tie(dp, cpot) = tupt;

    ld1 a = 0;
    ld1 b = dp / 2; 
    
    // ld1 tetta = 0;
    // crystaltickness = 20;

    // ld1 nbeam = 5;
    vector<ld1> approxim_x;
    approxim_x.clear();
    vector<ld1> evee;
    evee.clear();
	
	for (int i = 1; i <= nbeam; i++)
        evee.push_back(eve(i, cpot, dp, tetta, crystaltype));
	
	  ld1 step = b / (nbeam*10);
    ld1 ucr0 = ucr(0, cpot, dp, crystaltype);
    ld1 ucrdp_2 = ucr(b, cpot, dp, crystaltype); // dp/2
    int count = 0;
    for (int j = 0; j < nbeam; j++)
        for (int i = 0; i <nbeam*10; i++) //надо что-то написать про зависимость от числа шагов i<?
        {
            if ((ucr(step * i, cpot, dp, crystaltype) - evee[j]) * (ucr(step * (i + 1), cpot, dp, crystaltype) - evee[j]) < 0)
            {
                approxim_x.push_back((i + 0.5) * step);         //i * step + step/2
                count += 1;
            }
        }
	vector<ld1> x1, x2;
    x1.clear();
    x2.clear();
    int counter = 0;
	
	 for (int j = 0; j < nbeam; j++)
    {
            if (evee[j] < ucr0)
            {
                x1.push_back(find(approxim_x[counter], 0.000001, evee[j], cpot, dp, crystaltype));
                x2.push_back(find(approxim_x[counter + 1], 0.000001, evee[j], cpot, dp, crystaltype));
                counter += 2;
                //print(x1);
            }
            else
            {
                x1.push_back(0);
                if (ucr0 < evee[j] && evee[j] < ucrdp_2)
                {
                    x2.push_back(find(approxim_x[counter], 0.000001, evee[j], cpot, dp, crystaltype));
                    counter ++;
                }
                else if (evee[j] > ucrdp_2)
                    x2.push_back(dp / 2);
            }
        //if (counter == count)
         //   break;
    }
	
	// vector<ld1> TT;
    TT.clear();
	cout <<endl<< "x1 = ";
    print(x1);
    cout << endl << endl << endl;
    cout << "x2 = ";
    print(x2);
    cout << endl;
    for (int j = 0; j < nbeam; j++)
    {
        ld1 rrr = simpson(x1[j], x2[j], 1000, evee[j], cpot, dp, crystaltype);
        if (evee[j] < ucr0 || evee[j] > ucrdp_2)
            TT.push_back(0.5 / (c0 * pow(10, 10)) * 4 * sqrt(0.5 * m0 * gama) * rrr);
        else
            TT.push_back(1 / (c0 * pow(10, 10)) * 4 * sqrt(0.5 * m0 * gama) * rrr);
    }
	
	
    min_TT = *min_element(TT.begin(), TT.end());
 //   return {min_TT, TT};
}


//создать класс инпут параметров...

