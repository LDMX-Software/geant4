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
//
//         ------------ G4GammaConversionToMuons physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, April 2002
//
//
// 07-08-02: missprint in OR condition in DoIt : f1<0 || f1>f1_max ..etc ...
// 25-10-04: migrade to new interfaces of ParticleChange (vi)
// ---------------------------------------------------------------------------

#include "G4GammaConversionToMuons.hh"

#include "G4BetheHeitler5DModel.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4Exp.hh"
#include "G4Gamma.hh"
#include "G4Log.hh"
#include "G4LossTableManager.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Positron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

static const G4double sqrte = std::sqrt(std::exp(1.));
static const G4double PowSat = -0.88;

G4GammaConversionToMuons::G4GammaConversionToMuons(const G4String& processName,
						   G4ProcessType type)
  : G4VDiscreteProcess (processName, type),
    Mmuon(G4MuonPlus::MuonPlus()->GetPDGMass()),
    Rc(CLHEP::elm_coupling / Mmuon),
    LimitEnergy(5. * Mmuon),
    LowestEnergyLimit(2. * Mmuon),
    HighestEnergyLimit(1e12 * CLHEP::GeV),  // ok to 1e12GeV, then LPM suppression
    theGamma(G4Gamma::Gamma()),
    theMuonPlus(G4MuonPlus::MuonPlus()),
    theMuonMinus(G4MuonMinus::MuonMinus())
{
  SetProcessSubType(fGammaConversionToMuMu);
  fManager = G4LossTableManager::Instance();
  fManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4GammaConversionToMuons::~G4GammaConversionToMuons()
{
  fManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4GammaConversionToMuons::IsApplicable(const G4ParticleDefinition& part)
{
  return (&part == theGamma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GammaConversionToMuons::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  Energy5DLimit = G4EmParameters::Instance()->MaxEnergyFor5DMuPair();

  auto table = G4Material::GetMaterialTable();
  std::size_t nelm = 0;
  for (auto const& mat : *table) {
    std::size_t n = mat->GetNumberOfElements();
    nelm = std::max(nelm, n);
  }
  temp.resize(nelm, 0);

  if (Energy5DLimit > 0.0 && nullptr != f5Dmodel) {
    f5Dmodel = new G4BetheHeitler5DModel();
    f5Dmodel->SetLeptonPair(theMuonPlus, theMuonMinus);
    const std::size_t numElems = G4ProductionCutsTable::GetProductionCutsTable()->GetTableSize();
    const G4DataVector cuts(numElems);
    f5Dmodel->Initialise(&p, cuts);
  }
  PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversionToMuons::GetMeanFreePath(const G4Track& aTrack, G4double,
                                                   G4ForceCondition*)
// returns the photon mean free path in GEANT4 internal units
{
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4Material* aMaterial = aTrack.GetMaterial();
  return ComputeMeanFreePath(GammaEnergy, aMaterial);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4GammaConversionToMuons::ComputeMeanFreePath(G4double GammaEnergy,
                                              const G4Material* aMaterial)

// computes and returns the photon mean free path in GEANT4 internal units
{
  if(GammaEnergy <= LowestEnergyLimit) { return DBL_MAX; }
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

  G4double SIGMA = 0.0;
  G4double fact  = 1.0;
  G4double e = GammaEnergy;
  // low energy approximation as in Bethe-Heitler model
  if(e < LimitEnergy) {
    G4double y = (e - LowestEnergyLimit)/(LimitEnergy - LowestEnergyLimit);
    fact = y*y;
    e = LimitEnergy;
  } 

  for ( std::size_t i=0 ; i < aMaterial->GetNumberOfElements(); ++i)
  {
    SIGMA += NbOfAtomsPerVolume[i] * fact *
      ComputeCrossSectionPerAtom(e, (*theElementVector)[i]->GetZasInt());
  }
  return (SIGMA > 0.0) ? 1./SIGMA : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GammaConversionToMuons::GetCrossSectionPerAtom(
                                   const G4DynamicParticle* aDynamicGamma,
                                   const G4Element* anElement)

// gives the total cross section per atom in GEANT4 internal units
{
   return ComputeCrossSectionPerAtom(aDynamicGamma->GetKineticEnergy(),
                                     anElement->GetZasInt());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4GammaConversionToMuons::ComputeCrossSectionPerAtom(
                         G4double Egam, G4int Z)
			 
// Calculates the microscopic cross section in GEANT4 internal units.
// Total cross section parametrisation from H.Burkhardt
// It gives a good description at any energy (from 0 to 10**21 eV)
{ 
  if(Egam <= LowestEnergyLimit) { return 0.0; }

  G4NistManager* nist = G4NistManager::Instance();

  G4double PowThres, Ecor, B, Dn, Zthird, Winfty, WMedAppr, Wsatur, sigfac;

  if (Z == 1) {  // special case of Hydrogen
    B = 202.4;
    Dn = 1.49;
  }
  else {
    B = 183.;
    Dn = 1.54 * nist->GetA27(Z);
  }
  Zthird = 1. / nist->GetZ13(Z);  // Z**(-1/3)
  Winfty = B * Zthird * Mmuon / (Dn * electron_mass_c2);
  WMedAppr = 1. / (4. * Dn * sqrte * Mmuon);
  Wsatur = Winfty / WMedAppr;
  sigfac = 4. * fine_structure_const * Z * Z * Rc * Rc;
  PowThres = 1.479 + 0.00799 * Dn;
  Ecor = -18. + 4347. / (B * Zthird);

  G4double CorFuc = 1. + .04 * G4Log(1. + Ecor / Egam);
  G4double Eg =
    G4Exp(G4Log(1. - 4. * Mmuon / Egam) * PowThres)
    * G4Exp(G4Log(G4Exp(G4Log(Wsatur) * PowSat) + G4Exp(G4Log(Egam) * PowSat)) / PowSat);

  G4double CrossSection = 7. / 9. * sigfac * G4Log(1. + WMedAppr * CorFuc * Eg);
  CrossSection *= CrossSecFactor;  // increase the CrossSection by  (by default 1)
  return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4GammaConversionToMuons::SetCrossSecFactor(G4double fac)
// Set the factor to artificially increase the cross section
{
  if (fac < 0.0) return;
  CrossSecFactor = fac;
  if (verboseLevel > 1) {
    G4cout << "The cross section for GammaConversionToMuons is artificially "
           << "increased by the CrossSecFactor=" << CrossSecFactor << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VParticleChange* G4GammaConversionToMuons::PostStepDoIt(
                                                        const G4Track& aTrack,
                                                        const G4Step&  aStep)
//
// generation of gamma->mu+mu-
//
{
  aParticleChange.Initialize(aTrack);
  const G4Material* aMaterial = aTrack.GetMaterial();

  // current Gamma energy and direction, return if energy too low
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double Egam = aDynamicGamma->GetKineticEnergy();
  if (Egam <= LowestEnergyLimit) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  //
  // Kill the incident photon
  //
  aParticleChange.ProposeMomentumDirection( 0., 0., 0. ) ;
  aParticleChange.ProposeEnergy( 0. ) ;
  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;

  if (Egam <= Energy5DLimit) {
    std::vector<G4DynamicParticle*> fvect;
    f5Dmodel->SampleSecondaries(&fvect, aTrack.GetMaterialCutsCouple(), 
				aTrack.GetDynamicParticle(), 0.0, DBL_MAX);
    for(auto dp : fvect) { aParticleChange.AddSecondary(dp); }
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }  

  G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();

  // select randomly one element constituting the material
  const G4Element* anElement = SelectRandomAtom(aDynamicGamma, aMaterial);
  G4int Z = anElement->GetZasInt();
  G4NistManager* nist = G4NistManager::Instance();

  G4double xPlus = 0.5;
  G4double xMinus = 0.5;
  G4double thetaPlus = 0., thetaMinus = 0., phiHalf = 0.; // final angular variables

  // LDMX: full Tsai matrix element away from threshold
  if (Egam > LimitEnergy) {
    SampleTsai(Egam, nist->GetA27(Z), xPlus, thetaPlus, thetaMinus, phiHalf);
    xMinus = 1. - xPlus;
  } else {
    G4double B, Dn;
    G4double A027 = nist->GetA27(Z);

    if (Z == 1) {  // special case of Hydrogen
      B = 202.4;
      Dn = 1.49;
    }
    else {
      B = 183.;
      Dn = 1.54 * A027;
    }
    G4double Zthird = 1. / nist->GetZ13(Z);  // Z**(-1/3)
    G4double Winfty = B * Zthird * Mmuon / (Dn * electron_mass_c2);

    G4double C1Num = 0.138 * A027;
    G4double C1Num2 = C1Num * C1Num;
    G4double C2Term2 = electron_mass_c2 / (183. * Zthird * Mmuon);

    G4double GammaMuonInv = Mmuon / Egam;

    // generate xPlus according to the differential cross section by rejection
    G4double xmin = (Egam <= LimitEnergy) ? 0.5 : 0.5 - std::sqrt(0.25 - GammaMuonInv);
    G4double xmax = 1. - xmin;

    G4double Ds2 = (Dn * sqrte - 2.);
    G4double sBZ = sqrte * B * Zthird / electron_mass_c2;
    G4double LogWmaxInv =
      1. / G4Log(Winfty * (1. + 2. * Ds2 * GammaMuonInv) / (1. + 2. * sBZ * Mmuon * GammaMuonInv));
    G4double xPM = 0.25;

    G4int nn = 0;
    const G4int nmax = 1000;

    // sampling for Egam > LimitEnergy
    if (xmin < 0.5) {
      G4double result, W;
      do {
        xPlus = xmin + G4UniformRand() * (xmax - xmin);
        xMinus = 1. - xPlus;
        xPM = xPlus * xMinus;
        G4double del = Mmuon * Mmuon / (2. * Egam * xPM);
        W = Winfty * (1. + Ds2 * del / Mmuon) / (1. + sBZ * del);
        G4double xxp = 1. - 4. / 3. * xPM;  // the main xPlus dependence
        result = (xxp > 0.) ? xxp * G4Log(W) * LogWmaxInv : 0.0;
        if (result > 1.) {
          G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
                 << " in dSigxPlusGen, result=" << result << " > 1" << G4endl;
        }
        ++nn;
        if(nn >= nmax) { break; }
      }
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
      while (G4UniformRand() > result);
    }

    // now generate the angular variables via the auxilary variables t,psi,rho
    G4double t;
    G4double psi;
    G4double rho;

    G4double a3 = (GammaMuonInv / (2. * xPM));
    G4double a33 = a3 * a3;
    G4double f1;
    G4double b1  = 1./(4.*C1Num2);
    G4double b3  = b1*b1*b1;
    G4double a21 = a33 + b1;
  
    G4double f1_max=-(1.-xPM)*(2.*b1+(a21+a33)*G4Log(a33/a21))/(2*b3);  

    nn = 0;
    // t, psi, rho generation start  (while angle < pi)
    do {
      //generate t by the rejection method
      do { 
        ++nn;
        t=G4UniformRand();
        G4double a34=a33/(t*t);
        G4double a22 = a34 + b1;
        if(std::abs(b1)<0.0001*a34) {
          // special case of a34=a22 because of logarithm accuracy
          f1=(1.-2.*xPM+4.*xPM*t*(1.-t))/(12.*a34*a34*a34*a34);
        }
        else {
          f1=-(1.-2.*xPM+4.*xPM*t*(1.-t))*(2.*b1+(a22+a34)*G4Log(a34/a22))/(2*b3);
        }
        if (f1 < 0.0 || f1 > f1_max) {  // should never happend
          G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
          << "outside allowed range f1=" << f1
          << " is set to zero, a34 = "<< a34 << " a22 = "<<a22<<"."
          << G4endl;
          f1 = 0.0;
        }
        if(nn > nmax) { break; }
        // Loop checking, 07-Aug-2015, Vladimir Ivanchenko  
      } while ( G4UniformRand()*f1_max > f1);
      // generate psi by the rejection method
      G4double f2_max=1.-2.*xPM*(1.-4.*t*(1.-t));
      // long version
      G4double f2;
      do { 
        ++nn;
        psi=twopi*G4UniformRand();
        f2=1.-2.*xPM+4.*xPM*t*(1.-t)*(1.+std::cos(2.*psi));
        if(f2<0 || f2> f2_max) { // should never happend
          G4cout << "G4GammaConversionToMuons::PostStepDoIt WARNING:"
                 << "outside allowed range f2=" << f2 << " is set to zero" << G4endl;
          f2 = 0.0;
        }
        if(nn >= nmax) { break; }
        // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
      } while ( G4UniformRand()*f2_max > f2);

      // generate rho by direct transformation
      G4double C2Term1=GammaMuonInv/(2.*xPM*t);
      G4double C22 = C2Term1*C2Term1+C2Term2*C2Term2;
      G4double C2=4.*C22*C22/std::sqrt(xPM);
      G4double rhomax=(1./t-1.)*1.9/A027;
      G4double beta=G4Log( (C2+rhomax*rhomax*rhomax*rhomax)/C2 );
      rho=G4Exp(G4Log(C2 *( G4Exp(beta*G4UniformRand())-1. ))*0.25);

      //now get from t and psi the kinematical variables
      G4double u=std::sqrt(1./t-1.);
      G4double xiHalf=0.5*rho*std::cos(psi);
      phiHalf=0.5*rho/u*std::sin(psi);

      thetaPlus =GammaMuonInv*(u+xiHalf)/xPlus;
      thetaMinus=GammaMuonInv*(u-xiHalf)/xMinus;

      // protection against infinite loop
      if(nn > nmax) {
        if(std::abs(thetaPlus)>pi) { thetaPlus = 0.0; }
        if(std::abs(thetaMinus)>pi) { thetaMinus = 0.0; }
      }

      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while ( std::abs(thetaPlus)>pi || std::abs(thetaMinus) >pi);
  }

  // now construct the vectors
  // azimuthal symmetry, take phi0 at random between 0 and 2 pi
  G4double phi0=twopi*G4UniformRand(); 
  G4double EPlus=xPlus*Egam;
  G4double EMinus=xMinus*Egam;

  // mu+ mu- directions for gamma in z-direction
  G4ThreeVector MuPlusDirection  ( std::sin(thetaPlus) *std::cos(phi0+phiHalf),
                   std::sin(thetaPlus)  *std::sin(phi0+phiHalf), std::cos(thetaPlus) );
  G4ThreeVector MuMinusDirection (-std::sin(thetaMinus)*std::cos(phi0-phiHalf),
                  -std::sin(thetaMinus) *std::sin(phi0-phiHalf), std::cos(thetaMinus) );
  // rotate to actual gamma direction
  MuPlusDirection.rotateUz(GammaDirection);
  MuMinusDirection.rotateUz(GammaDirection);

  // create G4DynamicParticle object for the particle1
  auto aParticle1 = new G4DynamicParticle(theMuonPlus, MuPlusDirection, EPlus - Mmuon);
  aParticleChange.AddSecondary(aParticle1);
  // create G4DynamicParticle object for the particle2
  auto aParticle2 = new G4DynamicParticle(theMuonMinus, MuMinusDirection, EMinus - Mmuon);
  aParticleChange.AddSecondary(aParticle2);
  //  Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// LDMX: N. Toro, see ReleaseNotes/muonNotes.tex; units of muon mass

void G4GammaConversionToMuons::SampleTsai(G4double Egam, G4double A027,
                                          G4double& xPlus,
                                          G4double& thetaPlus,
                                          G4double& thetaMinus,
                                          G4double& phiHalf) const
{
  const G4double GammaMuonInv = Mmuon / Egam;
  const G4double gamma0 = 1. / GammaMuonInv;
  const G4double sqrtx = std::sqrt(0.25 - GammaMuonInv);
  const G4double xmin = 0.5 - sqrtx;
  const G4double xPRange = 2. * sqrtx;
  const G4int nmax = 1000;

  // Weights reach ~1.3 at TeV, below 1 at 1-10 GeV
  const G4double maxWeight = 1.5;

  // Scale of q/mMu at which the form factor sets in [PRM(6.115),muonNotes(10)]
  const G4double qFormFactor = 1. / (0.20 * A027);

  // log(rho) range scale, see muonNotes (8),(9)
  const G4double x4 = std::pow(2. * GammaMuonInv, 4.);
  const G4double aRangeNorm = G4Log((std::pow(qFormFactor, 4.) + x4) / x4);

  G4bool goodEvent = false;
  G4int nn = 0;
  do {
    ++nn;

    const G4double xP = xmin + G4UniformRand() * xPRange;
    const G4double xM = 1. - xP;
    const G4double delOverMmuon = Mmuon / (2. * Egam * xP * xM);

    // Sample t, psi and a = log(rho^4+c2^4)
    const G4double t = G4UniformRand();
    const G4double u = std::sqrt(1. / t - 1.);
    const G4double psi = CLHEP::twopi * G4UniformRand();

    const G4double rhoMax = std::min(2. * u / std::abs(std::cos(psi)),
                                     CLHEP::pi * u / std::abs(std::sin(psi)));
    const G4double c24 = std::pow(delOverMmuon / t, 4.);
    const G4double amin = G4Log(c24);
    const G4double aRange = G4Log(std::pow(rhoMax, 4.) + c24) - amin;
    const G4double a = amin + G4UniformRand() * aRange;
    const G4double rho = std::pow(G4Exp(a) - c24, 0.25);

    const G4double xiHalf = 0.5 * rho * std::cos(psi);
    const G4double phi = rho / u * std::sin(psi);
    const G4double uPlus = u + xiHalf;
    const G4double uMinus = u - xiHalf;

    // Reject angles beyond pi
    if (uPlus > xP * CLHEP::pi / GammaMuonInv ||
        uMinus > xM * CLHEP::pi / GammaMuonInv) { continue; }

    const G4double thP = GammaMuonInv * uPlus / xP;
    const G4double gammaPlus = gamma0 * xP;
    const G4double betaPlus = std::sqrt(1. - 1. / (gammaPlus * gammaPlus));
    const G4double kdotpPlus = gamma0 * gammaPlus * (1. - betaPlus * std::cos(thP));
    const G4double pPerpPlus = betaPlus * gammaPlus * std::sin(thP);

    const G4double thM = GammaMuonInv * uMinus / xM;
    const G4double gammaMinus = gamma0 * xM;
    const G4double betaMinus = std::sqrt(1. - 1. / (gammaMinus * gammaMinus));
    const G4double kdotpMinus = gamma0 * gammaMinus * (1. - betaMinus * std::cos(thM));
    const G4double pPerpMinus = betaMinus * gammaMinus * std::sin(thM);

    // Positive squared momentum transfer
    const G4double qpar = gamma0 - betaMinus * gammaMinus * std::cos(thM)
                          - betaPlus * gammaPlus * std::cos(thP);
    const G4double q2perp = (pPerpPlus - pPerpMinus) * (pPerpPlus - pPerpMinus)
                            + 2. * pPerpPlus * pPerpMinus * (1. - std::cos(phi));
    const G4double q2 = qpar * qpar + q2perp;

    // Measure and Jacobians from (a, t, psi) to (uPlus, uMinus, phi)
    const G4double integratedMeasure = xPRange * aRange / aRangeNorm;
    const G4double jacobian = G4Exp(a) / (rho * rho * rho) * rho / (t * (1. - t));

    // Tsai 1974 matrix element
    const G4double term1Plus = -(-q2 + kdotpMinus) / kdotpPlus;
    const G4double term1Minus = -(-q2 + kdotpPlus) / kdotpMinus;
    const G4double term2Plus = (-q2 / 2. + 2. * gammaMinus * gammaMinus) / (kdotpPlus * kdotpPlus);
    const G4double term2Minus = (-q2 / 2. + 2. * gammaPlus * gammaPlus) / (kdotpMinus * kdotpMinus);
    const G4double mixedTerm = 2. * ((1. + q2 / 2.) * (-2. * gammaMinus * gammaPlus - q2 / 2.)
                                     + q2 * gamma0 * gamma0 / 2.) / (kdotpMinus * kdotpPlus);
    const G4double matrixEl = std::sin(thP) * std::sin(thM) / (4. * q2 * q2)
      * (term1Plus + term1Minus + term2Plus + term2Minus + mixedTerm);

    const G4double ffSqrt = 1. / (1. + q2 / (qFormFactor * qFormFactor));
    const G4double weight = integratedMeasure * jacobian * matrixEl * ffSqrt * ffSqrt;

    if (weight < 0. || weight > maxWeight) {
      G4cerr << "G4GammaConversionToMuons::SampleTsai WARNING: weight " << weight
             << " out of bounds 0 to " << maxWeight << G4endl;
    }

    xPlus = xP;
    thetaPlus = thP;
    thetaMinus = thM;
    phiHalf = 0.5 * phi;

    if (G4UniformRand() * maxWeight < weight) { goodEvent = true; }
  } while (nn < nmax && !goodEvent);

  if (!goodEvent) {
    G4cerr << "G4GammaConversionToMuons::SampleTsai WARNING: no event accepted after "
           << nmax << " tries, using the last sampled kinematics" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element* G4GammaConversionToMuons::SelectRandomAtom(
		  const G4DynamicParticle* aDynamicGamma,
		  const G4Material* aMaterial)
{
  // select randomly 1 element within the material, invoked by PostStepDoIt

  const std::size_t NumberOfElements      = aMaterial->GetNumberOfElements();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4Element* elm = (*theElementVector)[0];

  if (NumberOfElements > 1) {
    G4double e = std::max(aDynamicGamma->GetKineticEnergy(), LimitEnergy);
    const G4double* natom = aMaterial->GetVecNbOfAtomsPerVolume();

    G4double sum = 0.;
    for (std::size_t i=0; i<NumberOfElements; ++i) {
      elm = (*theElementVector)[i];
      sum += natom[i]*ComputeCrossSectionPerAtom(e, elm->GetZasInt());
      temp[i] = sum;
    }
    sum *= G4UniformRand();
    for (std::size_t i=0; i<NumberOfElements; ++i) {
      if(sum <= temp[i]) {
        elm = (*theElementVector)[i];
        break;
      }
    }
  }
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4GammaConversionToMuons::PrintInfoDefinition()
{
  G4String comments = "gamma->mu+mu- Bethe Heitler process, SubType= ";
  G4cout << G4endl << GetProcessName() << ":  " << comments << GetProcessSubType() << G4endl;
  G4cout << "        good cross section parametrization from "
         << G4BestUnit(LowestEnergyLimit, "Energy") << " to " << HighestEnergyLimit / GeV
         << " GeV for all Z." << G4endl;
  G4cout << "        cross section factor: " << CrossSecFactor << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
