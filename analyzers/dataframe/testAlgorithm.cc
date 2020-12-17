//Sphericity (and later thrust) algorithms based on Aleph

#include "testAlgorithm.h"

std::vector<float> alg_sphericity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {

  // ------If number of particles is less than three: special treatment ??? NOT DONE YET

  //-- Compute momentum tensor
  TMatrixD TT(3,3);
  TT.Zero();

  //--------Test vectors
  vector< array<float,3> > vect = new vector< array<float,3> >;
  //first vector
  vect->push_back({1, 0, 0});
  //second vector
  vect->push_back({0, 1, 0});
  //third vector
  vect->push_back({0, 0, 1});
  // ------- 

  for (array<float, 3> p_i : vect) {
    for (int a=0;a<3;a++) {
      for (int b=0;b<3;b++) {
        S_tensor[a,b]+=p_i[a]*p_i[b];
      }
    }
  }

  /*for (unsigned int i=0; i<in.size() ; i++) {
    auto & p = in[i];
    float P(3);
    P(0)=p.momentum.x;
    P(1)=p.momentum.y;
    P(2)=p.momentum.z;
    for (int I=0;I<3;I++) {
      for (int k=0;K<I;K++){
	TT[I][K]+=P(I)*P(K);
      }
    }
    }*/
  TT[0][1]=TT[1][0];
  TT[0][2]=TT[2][0];
  TT[1][2]=TT[2][1];
  //CALL VSCALE(A,ALPHA,X,N)   X(I)=A(I)*ALPHA   (I=1,..,N)
  float Alpha=1./(TT[0][0]+TT[1][1]+TT[2][2]);
  for (int I=0;I<3;I++) {
    for (int k=0;K<3;K++){
      TT[I][K]*=Alpha;
    }
  }

  const TMatrixDEigen eigen(TT);
  TMatrixD eigenVal = eigen.GetEigenValues();
  TMatrixD eigenVec = eigen.GetEigenVectors();
  
  vector<float> EVAL;
  for (unsigned int i=0; i<3 ; i++) EVAL.push_back(eigenVal[i][i]);

  // Sort the vector in descending order
  vector<int> IND(3); //Saving index from sorting eigenvalues        
  size_t n(0);
  generate(begin(IND), end(IND), [&]{ return n++; });
  sort(begin(IND),end(IND),[&](int i1, int i2) { return EVAL[i1] > EVAL[i2]; } );
  
  ///--- DON´T KNOW HOW TO CROSS CHECK FOR IMAGINARY NUMBERS????
  // if (imaginary) continue;
  std::vector<float> major; //major axis (sphericity axis) 
  std::vector<float> semimajor; // semi-major axis
  std::vector<float> minor; //minor axis
  for (unsigned int i=0; i<3 ; i++){
    major.push_back(eigenVec(i,IND(0)));
    semimajor.push_back(eigenVec(i,IND(1)));
    minor.push_back(eigenVec(i,IND(2)));
      }
  float sphericity = 1.5*(1-EVAL(IND(0)));
  float aplanarity = 1.5*EVAL(IND(2));
  float planarity = EVAL(IND(2))/EVAL(IND(1));

  return major;
}
