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
#_fcc  = ROOT.getMC_px
#_fcc2  = ROOT.getRP2MC_p
#Not in Namespace (works)
_fcc  = ROOT.dummyloader

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)
#print ('fccana2  ',_fcc2)


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
               .Define("MC_px",         "MCParticle::get_px(Particle)")
               .Define("MC_py",         "MCParticle::get_py(Particle)")
               .Define("MC_pz",         "MCParticle::get_pz(Particle)")
               .Define("MC_p",          "MCParticle::get_p(Particle)")
               .Define("MC_pdg",        "MCParticle::get_pdg(Particle)")
               .Define("MC_charge",     "MCParticle::get_charge(Particle)")
               .Define("MC_mass",       "MCParticle::get_mass(Particle)")
               .Define("MC_e",          "MCParticle::get_e(Particle)")
               .Define("MC_status",     "MCParticle::get_genStatus(Particle)")
               #.Define("MC_vertex_x",   "MCParticle::get_vertex_x(Particle)")
               #.Define("MC_vertex_y",   "MCParticle::get_vertex_y(Particle)")
               #.Define("MC_vertex_z",   "MCParticle::get_vertex_z(Particle)")

               .Define("RP_p",          "getRP_p(ReconstructedParticles)")
               .Define("RP_px",         "getRP_px(ReconstructedParticles)")
               .Define("RP_py",         "getRP_py(ReconstructedParticles)")
               .Define("RP_pz",         "getRP_pz(ReconstructedParticles)")
               .Define("RP_charge",     "getRP_charge(ReconstructedParticles)")
               .Define("RP_mass",       "getRP_mass(ReconstructedParticles)")
               .Define("RP_e",          "getRP_e(ReconstructedParticles)")
               
               .Define("RPTRK_D0",      "getRP2TRK_D0(ReconstructedParticles, EFlowTrack_1)")
               .Define("RPTRK_Z0",      "getRP2TRK_D0(ReconstructedParticles, EFlowTrack_1)")

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
               .Define('RPMC_parentindex', "MCParticle::get_parentid(RPMC_index,Particle, Particle0)")
               .Define('MC_daughter1', "MCParticle::getMC_daughter(0,Particle,Particle1)")
               .Define('MC_daughter2', "MCParticle::getMC_daughter(1,Particle,Particle1)")
               .Define('MC_parent1', "MCParticle::getMC_parent(0,Particle,Particle0)")
               .Define('MC_parent2', "MCParticle::getMC_parent(1,Particle,Particle0)")

               
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
               .Define("RPrest_association", "ParticleSetAssociation(ReconstructedParticles, Set_rest)")
               .Define("RPrest_invmass", "RPsetInvariantMass(RPrest)")
               .Define("RPrest_px",      "getRP_px(RPrest)")
               .Define("RPrest_py",      "getRP_py(RPrest)")
               .Define("RPrest_pz",      "getRP_pz(RPrest)")
               .Define("RPrest_charge",      "getRP_charge(RPrest)")
               .Define("RPrest_e",          "getRP_e(RPrest)")
               #.Define('RPMC_p',        match,string_vec)

               .Define("Set_finalStates", "MCParticle::selector_finalStates(Particle, Set_lepton, RPMC_index)")
               .Define("Particle_finalStates", "MCParticle::ParticleSetCreator(Particle, Set_finalStates)")
               .Define("MCfinal_px", "MCParticle::get_px(Particle_finalStates)")
               .Define("MCfinal_py", "MCParticle::get_py(Particle_finalStates)")
               .Define("MCfinal_pz", "MCParticle::get_pz(Particle_finalStates)")
               .Define("MCfinal_e",  "MCParticle::get_e(Particle_finalStates)")
               
               .Define('EVT_thrust',      'minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVT_thrust_val',  'EVT_thrust.at(0)')
               .Define('EVT_thrust_x',    'EVT_thrust.at(1)')
               .Define('EVT_thrust_x_err','EVT_thrust.at(2)')
               .Define('EVT_thrust_y',    'EVT_thrust.at(3)')
               .Define('EVT_thrust_y_err','EVT_thrust.at(4)')
               .Define('EVT_thrust_z',    'EVT_thrust.at(5)')
               .Define('EVT_thrust_z_err','EVT_thrust.at(6)')
               
               #.Define('EVT_thrust',     'minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               #.Define('RP_thrustangle', 'axisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
               #.Define('EVT_thrust_x',   "EVT_thrust.at(0)")
               #.Define('EVT_thrust_y',   "EVT_thrust.at(1)")
               #.Define('EVT_thrust_z',   "EVT_thrust.at(2)")
               #.Define('EVT_thrust_val', "EVT_thrust.at(3)")
               
               .Define('EVT_sphericity',     'minimize_sphericity("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVT_sphericity_x',   "EVT_sphericity.at(0)")
               .Define('EVT_sphericity_y',   "EVT_sphericity.at(1)")
               .Define('EVT_sphericity_z',   "EVT_sphericity.at(2)")
               .Define('EVT_sphericity_val', "EVT_sphericity.at(3)")
               .Define('RP_sphericityangle', 'axisCosTheta(EVT_sphericity, RP_px, RP_py, RP_pz)')

               #.Define('RP_hemis0_mass',   "getAxisMass(0)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
               #.Define('RP_hemis1_mass',   "getAxisMass(1)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")

               .Define('EVTrest_thrust',     'minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVTrest_thrust_val',  'EVTrest_thrust.at(0)')
               .Define('EVTrest_thrust_x',    'EVTrest_thrust.at(1)')
               .Define('EVTrest_thrust_x_err','EVTrest_thrust.at(2)')
               .Define('EVTrest_thrust_y',    'EVTrest_thrust.at(3)')
               .Define('EVTrest_thrust_y_err','EVTrest_thrust.at(4)')
               .Define('EVTrest_thrust_z',    'EVTrest_thrust.at(5)')
               .Define('EVTrest_thrust_z_err','EVTrest_thrust.at(6)')

               #.Define('EVTrest_thrust',     'minimize_thrust("Minuit2","Migrad")(RPrest_px, RPrest_py, RPrest_pz)')
               #.Define('RPrest_thrustangle', 'axisCosTheta(EVTrest_thrust, RPrest_px, RPrest_py, RPrest_pz)')
               #.Define('EVTrest_thrust_x',   "EVTrest_thrust.at(0)")
               #.Define('EVTrest_thrust_y',   "EVTrest_thrust.at(1)")
               #.Define('EVTrest_thrust_z',   "EVTrest_thrust.at(2)")
               #.Define('EVTrest_thrust_val', "EVTrest_thrust.at(3)")

               #.Define('RPrest_hemis0_mass',   "getAxisMass(0)(RPrest_thrustangle, RPrest_e, RPrest_px, RPrest_py, RPrest_pz)")
               #.Define('RPrest_hemis1_mass',   "getAxisMass(1)(RPrest_thrustangle, RPrest_e, RPrest_px, RPrest_py, RPrest_pz)")

               .Define("AlgoSpher", "alg_sphericity(ReconstructedParticles)")

               #run jet clustering with all reconstructed particles. kt_algorithm, R=1, exclusive clustering, exactly 4

               .Define("jets1_pair",           "JetClustering::clustering(1, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets1", "jets1_pair[0].second")
               .Define("jets1_association", "jets1_pair[0].first")
               .Define("jets1_px",        "JetClustering::getJet_px(jets1)")
               .Define("jets1_py",        "JetClustering::getJet_py(jets1)")
               .Define("jets1_pz",        "JetClustering::getJet_pz(jets1)")
               .Define("jets1_e",         "JetClustering::getJet_e(jets1)")

               .Define("MCjets1_pair",           "JetClustering::clustering(1, 1.0, 3, 4)(MCfinal_px, MCfinal_py, MCfinal_pz, MCfinal_e)")
               .Define("MCjets1", "MCjets1_pair[0].second")
               .Define("MCjets1_association", "MCjets1_pair[0].first")
               .Define("MCjets1_px",        "JetClustering::getJet_px(MCjets1)")
               .Define("MCjets1_py",        "JetClustering::getJet_py(MCjets1)")
               .Define("MCjets1_pz",        "JetClustering::getJet_pz(MCjets1)")
               .Define("MCjets1_e",         "JetClustering::getJet_e(MCjets1)")
               
               .Define("jets2_pair",           "JetClustering::clustering(2, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets2", "jets2_pair[0].second")
               .Define("jets2_association", "jets2_pair[0].first")
               .Define("jets2_px",        "JetClustering::getJet_px(jets2)")
               .Define("jets2_py",        "JetClustering::getJet_py(jets2)")
               .Define("jets2_pz",        "JetClustering::getJet_pz(jets2)")
               .Define("jets2_e",        "JetClustering::getJet_e(jets2)")
               
               .Define("jets3_pair",           "JetClustering::clustering(3, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets3", "jets3_pair[0].second")
               .Define("jets3_association", "jets3_pair[0].first")
               .Define("jets3_px",        "JetClustering::getJet_px(jets3)")
               .Define("jets3_py",        "JetClustering::getJet_py(jets3)")
               .Define("jets3_pz",        "JetClustering::getJet_pz(jets3)")
               .Define("jets3_e",        "JetClustering::getJet_e(jets3)")

               .Define("jets4_pair",           "JetClustering::clustering(4, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets4", "jets4_pair[0].second")
               .Define("jets4_association", "jets4_pair[0].first")
               .Define("jets4_px",        "JetClustering::getJet_px(jets4)")
               .Define("jets4_py",        "JetClustering::getJet_py(jets4)")
               .Define("jets4_pz",        "JetClustering::getJet_pz(jets4)")
               .Define("jets4_e",        "JetClustering::getJet_e(jets4)")

               .Define("jets5_pair",           "JetClustering::clustering(5, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets5", "jets5_pair[0].second")
               .Define("jets5_association", "jets5_pair[0].first")
               .Define("jets5_px",        "JetClustering::getJet_px(jets5)")
               .Define("jets5_py",        "JetClustering::getJet_py(jets5)")
               .Define("jets5_pz",        "JetClustering::getJet_pz(jets5)")
               .Define("jets5_e",        "JetClustering::getJet_e(jets5)")

               .Define("jets6_pair",           "JetClustering::clustering(6, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets6", "jets6_pair[0].second")
               .Define("jets6_association", "jets6_pair[0].first")
               .Define("jets6_px",        "JetClustering::getJet_px(jets6)")
               .Define("jets6_py",        "JetClustering::getJet_py(jets6)")
               .Define("jets6_pz",        "JetClustering::getJet_pz(jets6)")
               .Define("jets6_e",        "JetClustering::getJet_e(jets6)")

               .Define("jets7_pair",           "JetClustering::clustering(7, 1.0, 3, 4)(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("jets7", "jets7_pair[0].second")
               .Define("jets7_association", "jets7_pair[0].first")
               .Define("jets7_px",        "JetClustering::getJet_px(jets7)")
               .Define("jets7_py",        "JetClustering::getJet_py(jets7)")
               .Define("jets7_pz",        "JetClustering::getJet_pz(jets7)")
               .Define("jets7_e",        "JetClustering::getJet_e(jets7)")


               )

        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "MC_px",
                "MC_py",
                "MC_pz",
                "MC_p",
                "MC_pdg",
                "MC_charge",
                "MC_mass",
                "MC_e",
                "MC_status",
                #"MC_vertex_x",
                #"MC_vertex_y",
                #"MC_vertex_z",

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

                "MC_daughter1",
                "MC_daughter2",
                "MC_parent1",
                "MC_parent2",
                
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
                "RPrest_association",
                
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

                #"jets1_pair",
                "jets1_px",
                "jets1_py",
                "jets1_pz",
                "jets1_e",
                "jets1_association",

                "MCjets1_px",
                "MCjets1_py",
                "MCjets1_pz",
                "MCjets1_e",
                "MCjets1_association",

                "jets2_px",
                "jets2_py",
                "jets2_pz",
                "jets2_e",
                "jets2_association",

                "jets3_px",
                "jets3_py",
                "jets3_pz",
                "jets3_e",
                "jets3_association",

                "jets4_px",
                "jets4_py",
                "jets4_pz",
                "jets4_e",
                "jets4_association",

                "jets5_px",
                "jets5_py",
                "jets5_pz",
                "jets5_e",
                "jets5_association",

                "jets6_px",
                "jets6_py",
                "jets6_pz",
                "jets6_e",
                "jets6_association",

                "jets7_px",
                "jets7_py",
                "jets7_pz",
                "jets7_e",
                "jets7_association",

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
    #outDir = '/eos/user/j/jutornda/FCCee/'+sys.argv[0].split('/')[1]+'/'
    outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
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
