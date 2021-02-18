import sys
import ROOT

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
ROOT.gSystem.Load("libfastjet")
ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.getMC_px
_fcc2  = ROOT.getRP2MC_p

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)
print ('fccana2  ',_fcc2)


class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        match=ROOT.getRP2MC_p_func()
        string_vec = ROOT.std.vector('string')()
        string_vec.push_back('MCRecoAssociations#0.index')
        string_vec.push_back('MCRecoAssociations#1.index')
        string_vec.push_back('ReconstructedParticles')
        string_vec.push_back('Particle')

        df2 = (self.df
               
               .Define("RP_p",          "getRP_p(ReconstructedParticles)")
               .Define("RP_px",         "getRP_px(ReconstructedParticles)")
               .Define("RP_py",         "getRP_py(ReconstructedParticles)")
               .Define("RP_pz",         "getRP_pz(ReconstructedParticles)")
               .Define("RP_charge",     "getRP_charge(ReconstructedParticles)")
               .Define("RP_mass",       "getRP_mass(ReconstructedParticles)")
               .Define("RP_e",          "getRP_e(ReconstructedParticles)")
               
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")

               .Define('RPMC_index',    "getRP2MC_index(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles)")
               #.Define('RPMC_p',        "getRP2MC_p(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define('RPMC_px',       "getRP2MC_px(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define('RPMC_py',       "getRP2MC_py(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define('RPMC_pz',       "getRP2MC_pz(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               .Define('RPMC_pdg',      "getRP2MC_pdg(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define('RPMC_charge',   "getRP2MC_charge(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define('RPMC_mass',     "getRP2MC_mass(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               .Define('RPMC_parentindex', "getMC_parentid(RPMC_index,Particle, Particle0)")
               
               .Define("MET_p",          "getRP_p(MissingET)")
               .Define("MET_px",         "getRP_px(MissingET)")
               .Define("MET_py",         "getRP_py(MissingET)")
               .Define("MET_pz",         "getRP_pz(MissingET)")
               .Define("MET_charge",     "getRP_charge(MissingET)")
               .Define("MET_mass",       "getRP_mass(MissingET)")
               .Define("MET_e",          "getRP_e(MissingET)")

               .Define("Set_lepton",     "selector_HighestEnergyLepton(ReconstructedParticles, RPMC_pdg)")
               .Define("RPlepton",       "ParticleSetCreator(ReconstructedParticles, Set_lepton)")
               .Define("RPlepton_p",     "getRP_p(RPlepton)")
               .Define("RPleptonMET",    "mergeParticles(RPlepton,MissingET)")
               .Define("RPleptonMET_invmass","RPsetInvariantMass(RPleptonMET)")
               
               .Define("Set_rest",       "selector_rest(ReconstructedParticles,Set_lepton)")
               .Define("RPrest",         "ParticleSetCreator(ReconstructedParticles, Set_rest)")
               .Define("RPrest_invmass", "RPsetInvariantMass(RPrest)")
               .Define("RPrest_px",      "getRP_px(RPrest)")
               .Define("RPrest_py",      "getRP_py(RPrest)")
               .Define("RPrest_pz",      "getRP_pz(RPrest)")
               .Define("RPrest_charge",      "getRP_charge(RPrest)")
               .Define("RPrest_e",          "getRP_e(RPrest)")
               #.Define('RPMC_p',        match,string_vec)

               .Define('EVT_thrust',      'minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVT_thrust_val',  'EVT_thrust.at(0)')

               .Define('EVTrest_thrust',      'minimize_thrust("Minuit2","Migrad")(RPrest_px, RPrest_py, RPrest_pz)')
               .Define('EVTrest_thrust_val',  'EVTrest_thrust.at(0)')

               
               #.Define('EVT_thrust',     'minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               #.Define('RP_thrustangle', 'axisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
               #.Define('EVT_thrust_x',   "EVT_thrust.at(0)")
               #.Define('EVT_thrust_y',   "EVT_thrust.at(1)")
               #.Define('EVT_thrust_z',   "EVT_thrust.at(2)")
               #.Define('EVT_thrust_val', "EVT_thrust.at(3)")
               
               #.Define('EVT_sphericity',     'minimize_sphericity("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               #.Define('EVT_sphericity_x',   "EVT_sphericity.at(0)")
               #.Define('EVT_sphericity_y',   "EVT_sphericity.at(1)")
               #.Define('EVT_sphericity_z',   "EVT_sphericity.at(2)")
               #.Define('EVT_sphericity_val', "EVT_sphericity.at(3)")
               #.Define('RP_sphericityangle', 'axisCosTheta(EVT_sphericity, RP_px, RP_py, RP_pz)')

               #.Define('RP_hemis0_mass',   "getAxisMass(0)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
               #.Define('RP_hemis1_mass',   "getAxisMass(1)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
               

               #.Define('EVTrest_thrust',     'minimize_thrust("Minuit2","Migrad")(RPrest_px, RPrest_py, RPrest_pz)')
               #.Define('RPrest_thrustangle', 'axisCosTheta(EVTrest_thrust, RPrest_px, RPrest_py, RPrest_pz)')
               #.Define('EVTrest_thrust_x',   "EVTrest_thrust.at(0)")
               #.Define('EVTrest_thrust_y',   "EVTrest_thrust.at(1)")
               #.Define('EVTrest_thrust_z',   "EVTrest_thrust.at(2)")
               #.Define('EVTrest_thrust_val', "EVTrest_thrust.at(3)")

               #.Define('RPrest_hemis0_mass',   "getAxisMass(0)(RPrest_thrustangle, RPrest_e, RPrest_px, RPrest_py, RPrest_pz)")
               #.Define('RPrest_hemis1_mass',   "getAxisMass(1)(RPrest_thrustangle, RPrest_e, RPrest_px, RPrest_py, RPrest_pz)")

               #.Define("AlgoSpher", "alg_sphericity(ReconstructedParticles)")



               )

        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "RP_p",
                "RP_px",
                "RP_py",
                "RP_pz",
                "RP_charge",
                "RP_mass",
                "RP_e",

                #"RPTRK_D0",
                #"RPTRK_Z0",

                #"RPMC_p",
                #"RPMC_px",
                #"RPMC_py",
                #"RPMC_pz",
                #"RPMC_charge",
                #"RPMC_mass",
                "RPMC_pdg",
                #"RPMC_charge",
                "RPMC_index",
                "RPMC_parentindex",

                
                "MET_p",
                "MET_px",
                "MET_py",
                "MET_pz",
                "MET_charge",
                "MET_mass",
                "MET_e",

                "RPlepton_p",
                "RPleptonMET_invmass",
                "RPrest_invmass",

                
                #"event_thrust_x",
                #"event_thrust_y",
                #"event_thrust_z",
                "EVT_thrust_val",
                #"event_thrust",
                #"event_hemis_0",
                #"event_hemis_1",

                #"eventRest_thrust_x",
                #"eventRest_thrust_y",
                #"eventRest_thrust_z",
                "EVTrest_thrust_val",
                #"eventRest_thrust",
                #"eventRest_hemis_0",
                #"eventRest_hemis_1",

                #"AlgoSpher",

                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)

# example call for standalone file
# python FCCeeAnalyses/Z_Zbb_Flavor/dataframe/analysis.py /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp/p8_ee_Ztautau_ecm91/events_012154460.root

if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    outDir = '/eos/user/j/jutornda/FCCee/newThrust/'+sys.argv[0].split('/')[1]+'/'
    #outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 0
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
