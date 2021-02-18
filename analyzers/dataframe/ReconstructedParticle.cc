#include "ReconstructedParticle.h"


//TOBEMOVED LATER
ResonanceBuilder::ResonanceBuilder(int arg_resonance_pdgid, float arg_resonance_mass) {m_resonance_pdgid = arg_resonance_pdgid; m_resonance_mass = arg_resonance_mass;}
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> ResonanceBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      //reso.pdg = m_resonance_pdgid;
      TLorentzVector reso_lv; 
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
    std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}



recoil::recoil(float arg_sqrts) : m_sqrts(arg_sqrts) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  recoil::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  for (auto & v1: in) {
    TLorentzVector tv1;
    tv1.SetXYZM(v1.momentum.x, v1.momentum.y, v1.momentum.z, v1.mass);
    recoil_p4 -= tv1;
  }
  auto recoil_fcc = edm4hep::ReconstructedParticleData();
  recoil_fcc.momentum.x = recoil_p4.Px();
  recoil_fcc.momentum.y = recoil_p4.Py();
  recoil_fcc.momentum.z = recoil_p4.Pz();
  recoil_fcc.mass = recoil_p4.M();
  result.push_back(recoil_fcc);
  return result;
};

ROOT::VecOps::RVec<float> getRP_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
 ROOT::VecOps::RVec<float> result;
 for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(sqrt(in[i].momentum.x * in[i].momentum.x + in[i].momentum.y * in[i].momentum.y));
 }
 return result;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> mergeParticles(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y) {
  //to be keept as ROOT::VecOps::RVec
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(x.size() + y.size());
  result.insert( result.end(), x.begin(), x.end() );
  result.insert( result.end(), y.begin(), y.end() );
  return ROOT::VecOps::RVec(result);
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getRP(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i = 0; i < index.size(); ++i) {
    if (index[i]>-1)
      result.push_back(in.at(index[i]));
    //else
    //  std::cout << "electron index negative " << index[i]<<std::endl;
  }  
  return result;
}


ROOT::VecOps::RVec<float> getRP_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.mass);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    ROOT::Math::PxPyPzMVector lv(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(lv.E());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.P());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.x);
  }
  return result;
}


ROOT::VecOps::RVec<float> getRP_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.charge);
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Rapidity());
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> getRP_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv);
  }
  return result;
}


int getRP_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x) {
  int result =  x.size();
  return result;
}

selRP_pT::selRP_pT(float arg_min_pt) : m_min_pt(arg_min_pt) {};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  selRP_pT::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (std::sqrt(std::pow(p.momentum.x,2) + std::pow(p.momentum.y,2)) > m_min_pt) {
      result.emplace_back(p);
    }
  }
  return result;
}

selRP_p::selRP_p(float arg_min_p) : m_min_p(arg_min_p) {};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  selRP_p::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (std::sqrt(std::pow(p.momentum.x,2) + std::pow(p.momentum.y,2) + std::pow(p.momentum.z,2) ) > m_min_p) {
      result.emplace_back(p);
    }
  }
  return result;
}

selRP_charge::selRP_charge(int arg_charge, bool arg_abs){m_charge = arg_charge; m_abs = arg_abs;};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  selRP_charge::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ((m_abs && abs(in[i].charge)==m_charge) || (m_charge==in[i].charge) ) {
      result.emplace_back(p);
    }
  }
  return result;
}

std::unordered_set<int> selector_HighestEnergyLepton(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float> RP2MC_pdg) {
  std::unordered_set<int> result;
  UInt_t nLepton=0;
  float max=0;
  int idx;
  for (size_t i = 0; i < in.size(); ++i) {
    if (abs(RP2MC_pdg[i])==11 || abs(RP2MC_pdg[i])==13) {
      nLepton++;
      auto & p =in[i];
      ROOT::Math::PxPyPzMVector lv(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      if (nLepton==1) {
	max=lv.E();
	idx=i;
      }
      if (lv.E() > max) {
	max=lv.E();
	idx=i;
      }
    }
  }
  result.insert(idx);
  return result;
}

std::unordered_set<int> selector_rest(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx) {
  std::unordered_set<int> result;
  for (size_t i = 0; i < in.size(); ++i) {
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) result.insert(i);
  }
  return result;
}

std::vector<edm4hep::ReconstructedParticleData> ParticleSetCreator(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx) {
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.emplace_back(p);
  }
  return result;
}

ROOT::VecOps::RVec<int> ParticleSetAssociation(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,std::unordered_set<int> idx) {
  ROOT::VecOps::RVec<int> result{};
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.push_back(i);
  }
  return result;
}

ROOT::VecOps::RVec<bool> getJet_btag(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ParticleIDData> pid, ROOT::VecOps::RVec<float> values){
  ROOT::VecOps::RVec<bool> result;
  //std::cout << "========================new event=======================" <<std::endl;
  for (size_t i = 0; i < index.size(); ++i) {
    result.push_back(values.at(pid.at(index.at(i)).parameters_begin +1));
    
    //std::cout << pid.at(index.at(i)).parameters_begin << "  ==  " << pid.at(index.at(i)).parameters_end << std::endl;
    //for (unsigned j = pid.at(index.at(i)).parameters_begin; j != pid.at(index.at(i)).parameters_end; ++j) {
    //  std::cout << " values : " << values.at(j) << std::endl;
    //}
  }
  return result;
}

int getJet_ntags(ROOT::VecOps::RVec<bool> in) {
  int result =  0;
  for (size_t i = 0; i < in.size(); ++i)
    if (in.at(i))result+=1;
  return result;
}


getAxisRP::getAxisRP(bool arg_pos): m_pos(arg_pos) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getAxisRP::operator()(ROOT::VecOps::RVec<float> angle, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i = 0; i < angle.size(); ++i) {
    if (m_pos==1 && angle.at(i)>0.) result.push_back(in.at(i));
    if (m_pos==0 && angle.at(i)<0.) result.push_back(in.at(i));;
  }
  return result;
}


float RPsetInvariantMass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    ROOT::Math::PxPyPzMVector lv(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    E+=lv.E();
    px+=p.momentum.x;
    py+=p.momentum.y;
    pz+=p.momentum.z;
  }
  float result = sqrt(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2));
  return result;
  
}
/*
float RPsetInvariantMass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    E+=tlv.E();
    px+=p.momentum.x;
    py+=p.momentum.y;
    pz+=p.momentum.z;
  }
  float result = sqrt(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2));
  return result;

}
*/


std::vector<float> alg_sphericity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {

  // ------If number of particles is less than three: special treatment ??? NOT DONE YET

  //-- Compute momentum tensor
  TMatrixD TT(3,3);
  TT.Zero();

  /*
  //--------Test vectors
  std::vector< std::array<float,3> > vect;
  //first vector
  vect.push_back({1, 0, 0});
  //second vector
  vect.push_back({0, 1, 0});
  //third vector
  vect.push_back({0, 0, 1});

  for (std::array<float, 3> p_i : vect) {
    for (int a=0;a<3;a++) {
      for (int b=0;b<3;b++) {
        TT[a][b]+=p_i[a]*p_i[b];
      }
    }
  }
  // ------- 
  */
  
  for (unsigned int i=0; i<in.size() ; i++) {
    auto & p = in[i];
    float P[3];
    P[0]=p.momentum.x;
    P[1]=p.momentum.y;
    P[2]=p.momentum.z;
    for (int I=0;I<3;I++) {
      for (int K=0;K<I;K++){
        TT[I][K]+=P[I]*P[K];
      }
    }
  }
  
  TT[0][1]=TT[1][0];
  TT[0][2]=TT[2][0];
  TT[1][2]=TT[2][1];

  //CALL VSCALE(A,ALPHA,X,N)   X(I)=A(I)*ALPHA   (I=1,..,N)
  float Alpha=1./(TT[0][0]+TT[1][1]+TT[2][2]);
  for (int I=0;I<3;I++) {
    for (int K=0;K<3;K++){
      TT[I][K]*=Alpha;
    }
  }
const TMatrixDEigen eigen(TT);
  TMatrixD eigenVal = eigen.GetEigenValues();
  TMatrixD eigenVec = eigen.GetEigenVectors();
  
  std::vector<float> EVAL;
  for (unsigned int i=0; i<3 ; i++) EVAL.push_back(eigenVal[i][i]);

  // Sort the vector in descending order
  std::vector<int> IND(3); //Saving index from sorting eigenvalues        
  size_t n(0);
  generate(begin(IND), end(IND), [&]{ return n++; });
  sort(begin(IND),end(IND),[&](int i1, int i2) { return EVAL[i1] > EVAL[i2]; } );
  
  ///--- DON<C2><B4>T KNOW HOW TO CROSS CHECK FOR IMAGINARY NUMBERS????
  // if (imaginary) continue;
  std::vector<float> major; //major axis (sphericity axis) 
  std::vector<float> semimajor; // semi-major axis
  std::vector<float> minor; //minor axis
  for (unsigned int i=0; i<3 ; i++){
    major.push_back(eigenVec[i][IND[0]]);
    semimajor.push_back(eigenVec[i][IND[1]]);
    minor.push_back(eigenVec[i][IND[2]]);
      }
  float sphericity = 1.5*(1-EVAL[IND[0]]);
  float aplanarity = 1.5*EVAL[IND[2]];
  float planarity = EVAL[IND[2]]/EVAL[IND[1]];
  /*
  std::cout <<"major axis (sphericity axis) = (" << major[0] << ", " << major[1] << ", " << major[2] << ")" << std::endl;
  std::cout <<"semi-major axis = (" << semimajor[0] << ", " << semimajor[1] << ", " << semimajor[2] << ")" << std::endl;
  std::cout <<"minor axis = (" << minor[0] << ", " << minor[1] << ", " << minor[2] << ")" << std::endl;
  std::cout << "sphericity = " << sphericity << std::endl;
  std::cout << "aplanarity = " << aplanarity << std::endl;
  std::cout << "planarity = " << planarity << std::endl;
  */
  
  return major;
}

