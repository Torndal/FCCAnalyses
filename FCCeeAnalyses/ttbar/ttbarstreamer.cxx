// To compile from FCCAnalyses: g++ -o FCCeeAnalyses/ttbar/streamer FCCeeAnalyses/ttbar/ttbarstreamer.cxx -I analyzers/dataframe/ -L /afs/cern.ch/user/j/jutornda/newtutorial/FCCAnalyses/install/lib -lFCCAnalyses -std=gnu++17 `root-config --cflags --libs`
//To exceute: ./FCCeeAnalyses/ttbar/streamer FCCeeAnalyses/ttbar/test.root ~/mytutorial/FCCtop/p8_ee_ttbar_ecm350.root

#include <ROOT/RDataFrame.hxx>
#include "TLorentzVector.h"
#include <TSystem.h>

// FCC event datamodel includes
#include "datamodel/ParticleData.h"
#include "datamodel/LorentzVector.h"
#include "datamodel/JetData.h"
#include "datamodel/TaggedParticleData.h"
#include "datamodel/TaggedJetData.h"
// Legacy class
#include "datamodel/FloatData.h"

#include "FCCAnalyses.h"


auto _m = fcc::ParticleData();




// Reproduce Heppy analysis
int main(int argc, char* argv[]){


#ifdef ENABLEIMPLICITMT
  ROOT::EnableImplicitMT();
#endif
  
  // fcc edm libraries
  gSystem->Load("libdatamodel.so");
  gSystem->Load("libFCCAnalyses.so");
  
  // very basic command line argument parsing
  if (argc < 3) {
    std::cout << "error: need to specify fcc data files to analyze as command line arguments" << std::endl;
    std::cout << "usage:  fccanalysis_tth_4l outfilename.root datafile1.root datafile2.root ... datafileN.root " << std::endl;
    return 1;
  }
  std::cout << "Read files... ";
  std::vector<std::string> filenames;
  
  std::string outfilename = argv[1];
  for (int i = 2; i < argc; ++i) {
    std::cout << " " << argv[i];
    filenames.push_back(argv[i]);
  }
  std::cout << std::endl;
  
  std::cout << "Creating TDataFrame ..." << std::endl;
  ROOT::RDataFrame df("events", filenames);
  
  
  std::cout << "Apply selectors and define new branches ..." << std::endl;
  auto selectors =  df
    //.Define("selected_electrons", selectParticlesPtIso(20, 0.4), {"electrons", "electronITags"})
    .Define("muons_charge", "muons.core.charge")
    .Define("muons_px","muons.core.p4.px")
    .Define("muons_py","muons.core.p4.py")
    .Define("muons_pz","muons.core.p4.pz")
    .Define("muons_E", get_e, {"muons"})
    .Define("muons_p", get_p, {"muons"})
    .Define("muons_y", get_y, {"muons"})
    .Define("muons_mctruth_px", get_px_mctruth, {"genParticles", "muonsToMC#1.index"})
    .Define("muons_mctruth_py", get_py_mctruth, {"genParticles", "muonsToMC#1.index"})
    .Define("muons_mctruth_pz", get_pz_mctruth, {"genParticles", "muonsToMC#1.index"})
    .Define("muons_mctruth_charge", get_charge_mctruth, {"genParticles", "muonsToMC#1.index"})
    .Define("muons_index", get_mcindex, {"muonsToMC#1.index"})
    
    .Define("electrons_charge", "electrons.core.charge")
    .Define("electrons_px","electrons.core.p4.px")
    .Define("electrons_py","electrons.core.p4.py")
    .Define("electrons_pz","electrons.core.p4.pz")
    .Define("electrons_E", get_e, {"electrons"})
    .Define("electrons_p", get_p, {"electrons"})
    .Define("electrons_y", get_y, {"electrons"})
    .Define("electrons_mctruth_px", get_px_mctruth, {"genParticles", "electronsToMC#1.index"})
    .Define("electrons_mctruth_py", get_py_mctruth, {"genParticles", "electronsToMC#1.index"})
    .Define("electrons_mctruth_pz", get_pz_mctruth, {"genParticles", "electronsToMC#1.index"})
    .Define("electrons_mctruth_charge", get_charge_mctruth, {"genParticles", "electronsToMC#1.index"})
    .Define("electrons_index", get_mcindex, {"electronsToMC#1.index"})

    .Define("leptons", mergeParticles, {"electrons", "muons"})
    .Define("leptons_E", get_e, {"leptons"})    
    .Define("leptons_px", get_px, {"leptons"})
    .Define("leptons_py", get_py, {"leptons"})
    .Define("leptons_pz", get_pz, {"leptons"})
    
    .Define("photons_charge", "photons.core.charge")
    .Define("photons_px","photons.core.p4.px")
    .Define("photons_py","photons.core.p4.py")
    .Define("photons_pz","photons.core.p4.pz")
    .Define("photons_E", get_e, {"photons"})
    .Define("photons_p", get_p, {"photons"})
    .Define("photons_y", get_y, {"photons"})
    .Define("photons_mctruth_px", get_px_mctruth, {"genParticles", "photonsToMC#1.index"})
    .Define("photons_mctruth_py", get_py_mctruth, {"genParticles", "photonsToMC#1.index"})
    .Define("photons_mctruth_pz", get_pz_mctruth, {"genParticles", "photonsToMC#1.index"})
    .Define("photons_mctruth_charge", get_charge_mctruth, {"genParticles", "photonsToMC#1.index"})
    .Define("photons_index", get_mcindex, {"photonsToMC#1.index"})

    
    .Define("jets_area", "jets.core.area")
    .Define("jets_mass", "jets.core.p4.mass")
    .Define("jets_px", "jets.core.p4.px")
    .Define("jets_py", "jets.core.p4.py")
    .Define("jets_pz", "jets.core.p4.pz")
    .Define("jets_Flavor", "jetsFlavor.tag")
    .Define("jets_index", get_mcindex, {"jets#0.index"})

    .Define("btags", "bTags.tag")
    .Define("btags_index", get_mcindex, {"bTags#0.index"})
    .Define("ctags", "cTags.tag")
    .Define("ctags_index", get_mcindex, {"cTags#0.index"})
    .Define("tautags", "tauTags.tag")
    .Define("tautags_index", get_mcindex, {"tauTags#0.index"})
    
    .Define("met_x", "met.position.x")
    .Define("met_y", "met.position.y")
    
    //genParticles = truth                                                                                                                                   
    .Define("genParticles_charge", "genParticles.core.charge")
    //.DefineSlotEntry("genParticles_charge", "genParticles.core.charge")
    .Define("genParticles_mass", "genParticles.core.p4.mass")
    .Define("genParticles_px", "genParticles.core.p4.px")
    .Define("genParticles_py", "genParticles.core.p4.py")
    .Define("genParticles_pz", "genParticles.core.p4.pz")
    .Define("genParticles_pdgId", "genParticles.core.pdgId")
    .Define("genParticles_status", "genParticles.core.status")
    .Define("genParticles_x_vtx", "genParticles.core.vertex.x")
    .Define("genParticles_y_vtx", "genParticles.core.vertex.y")
    .Define("genParticles_z_vtx", "genParticles.core.vertex.z")
    .Define("genParticles_index0", get_mcindex, {"genParticles#0.index"})
    .Define("genParticles_index1", get_mcindex, {"genParticles#1.index"})

    //ef = energy flow                                                                                                                                       
    .Define("efcharged_mass", "efcharged.core.p4.mass")
    .Define("efcharged_px", "efcharged.core.p4.px")
    .Define("efcharged_py", "efcharged.core.p4.py")
    .Define("efcharged_pz", "efcharged.core.p4.pz")
    .Define("efcharged_pdgId", "efcharged.core.pdgId")
    .Define("efcharged_status", "efcharged.core.status")
    .Define("efcharged_x_vtx", "efcharged.core.vertex.x")
    .Define("efcharged_y_vtx", "efcharged.core.vertex.y")
    .Define("efcharged_z_vtx", "efcharged.core.vertex.z")
    
    .Define("efphotons_mass", "efphotons.core.p4.mass")
    .Define("efphotons_px", "efphotons.core.p4.px")
    .Define("efphotons_py", "efphotons.core.p4.py")
    .Define("efphotons_pz", "efphotons.core.p4.pz")
    .Define("efphotons_pdgId", "efphotons.core.pdgId")
    .Define("efphotons_status", "efphotons.core.status")
    .Define("efphotons_x_vtx", "efphotons.core.vertex.x")
    .Define("efphotons_y_vtx", "efphotons.core.vertex.y")
    .Define("efphotons_z_vtx", "efphotons.core.vertex.z")
    
    .Define("efneutrals_mass", "efneutrals.core.p4.mass")
    .Define("efneutrals_px", "efneutrals.core.p4.px")
    .Define("efneutrals_py", "efneutrals.core.p4.py")
    .Define("efneutrals_pz", "efneutrals.core.p4.pz")
    .Define("efneutrals_pdgId", "efneutrals.core.pdgId")
    .Define("efneutrals_status", "efneutrals.core.status")
    .Define("efneutrals_x_vtx", "efneutrals.core.vertex.x")
    .Define("efneutrals_y_vtx", "efneutrals.core.vertex.y")
    .Define("efneutrals_z_vtx", "efneutrals.core.vertex.z")
    ;
  
  auto nentries = selectors.Count();
  std::cout << "Count events: " <<  *nentries << std::endl;
  std::cout << "Writing snapshot to disk ... \t" << outfilename << std::endl;
  selectors.Snapshot("events", outfilename,
		     { 
		       // fcc particles with additional infos
		       "muons_charge",
		       "muons_px",
		       "muons_py",
		       "muons_pz",
		       "muons_E",
		       "muons_p",
		       "muons_y",
		       "muons_mctruth_px",
		       "muons_mctruth_py",
		       "muons_mctruth_pz",
		       "muons_mctruth_charge",
		       "muons_index",
		       
		       "electrons_charge",
		       "electrons_px",
		       "electrons_py",
		       "electrons_pz",
		       "electrons_E",
		       "electrons_p",
		       "electrons_y",
		       "electrons_mctruth_px",
		       "electrons_mctruth_py",
		       "electrons_mctruth_pz",
		       "electrons_mctruth_charge",
		       "electrons_index",			 

		       "leptons_E",
		       "leptons_px",
		       "leptons_py",
		       "leptons_pz",
		       
		       "photons_charge",
		       "photons_px",
		       "photons_py",
		       "photons_pz",
		       "photons_E",
		       "photons_p",
		       "photons_y",
		       "photons_mctruth_px",
		       "photons_mctruth_py",
		       "photons_mctruth_pz",
		       "photons_mctruth_charge",
		       "photons_index",

		       "jets_area",
		       "jets_mass",
		       "jets_px",
		       "jets_py",
		       "jets_pz",
		       "jets_Flavor",
		       "jets_index",
		       "btags",
		       "btags_index",
		       "ctags",
		       "ctags_index",
		       "tautags",
		       "tautags_index",
		       
		       "met_x",
		       "met_y",
		       
		       "genParticles_charge",
		       "genParticles_mass",
		       "genParticles_px",
		       "genParticles_py",
		       "genParticles_pz",
		       "genParticles_pdgId",
		       "genParticles_status",
		       "genParticles_x_vtx",
		       "genParticles_y_vtx",
		       "genParticles_z_vtx",
		       "genParticles_index0",
		       "genParticles_index1",
		       
		       "efcharged_mass",
		       "efcharged_px",
		       "efcharged_py",
		       "efcharged_pz",
		       "efcharged_pdgId",
		       "efcharged_status",
		       "efcharged_x_vtx",
		       "efcharged_y_vtx",
		       "efcharged_z_vtx",
		       
		       "efphotons_mass",
		       "efphotons_px",
		       "efphotons_py",
		       "efphotons_pz",
		       "efphotons_pdgId",
		       "efphotons_status",
		       "efphotons_x_vtx",
		       "efphotons_y_vtx",
		       "efphotons_z_vtx",
		       
		       "efneutrals_mass",
		       "efneutrals_px",
		       "efneutrals_py",
		       "efneutrals_pz",
		       "efneutrals_pdgId",
		       "efneutrals_status",
		       "efneutrals_x_vtx",
		       "efneutrals_y_vtx",
		       "efneutrals_z_vtx",
		       
		     }
		     );
  
  return 0;
}
