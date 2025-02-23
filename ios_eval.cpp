#include <iostream>
#include <fstream>
#include "include.hpp"
using namespace std;


extern void eval_test(){
  int num_vec = 1;                        // number of photos for each person
  double threshold = 0.468;
  double scaler0 = 40;                         // scaler for input vector
  double scaler1 = 40;                       // Scaler for vectors in the server

  int64_t scale_threshold = (int64_t)(threshold*scaler0*scaler1);

  const int degree = 1024;//2048;                        //polynomial degree 
  const int64_t ct_modulus = 3221225473;//206158430209;//2748779069441;//576460752154525697;           //ciphertext modulus
  const int64_t modulus = 6000;//12289;                       //plaintext modulus

  //int64_t ct_modulus = 576460752154525697;
  int64_t delta = 1;
  int logq = 32;//38;//42;//60;
  //int degree = 2048;
  double var = 3.2;
  int64_t b = 2;
  int logb = 1;
  //int64_t modulus = 12289; 
  int64_t primitive_root = 3;

  //read input ct
  //vector<vector<int64_t> > ct(2,vector<int64_t>(degree,0));
  ifstream fin;
  fin.open("ciphertext.txt");
  if(!fin.is_open()){
      cout <<"Cannot open file ciphertext.txt"<<endl;
    }
  
  vector<vector<int64_t> > ct(2,vector<int64_t>(degree,0));
  for(int i = 0 ; i < 2 ; ++i){
    for(int j = 0 ; j < degree ; ++j){
	     fin >>ct[i][j];
    }
  }
  fin.close();
  cout <<"Read input ciphertext from ciphertext.txt"<<endl;

  //read encrypted vectors
  vector<vector<vector<int64_t> > > ct_server(num_vec, vector<vector<int64_t>>(2,vector<int64_t>(degree,0)));
 

  fin.open("test_enc_vec.txt");
  if(!fin.is_open()){
      cout <<"Cannot open file test_enc_vec.txt"<<endl;
    }
  for(int k = 0 ; k < num_vec ; ++k){
    for(int i = 0 ; i < 2 ; ++i){
      for(int j = 0 ; j < degree ; ++j){
        fin >>ct_server[k][i][j];
      }
    }
  }
  
  fin.close();
  cout <<"read encrypted registered vectors from test_enc_vec.txt, Number of vectors = "<<num_vec<<endl;
  
  //read RelK
  vector<vector<polynomial>> RelK(logq, vector<polynomial>(2, polynomial(degree,0)));
  fin.open("rlwe_relk.txt");
  if(!fin.is_open()){
      cout <<"Cannot open file rlwe_relk.txt"<<endl;
    }
  for (int i = 0; i < logq; ++i){
    for (int j = 0; j < 2; ++j){
      for(int k = 0 ; k < degree ; ++k){
        fin >>RelK[i][j][k];
      }
    }
  }
  fin.close();


  vector<vector<int64_t>> lwe_ip(num_vec,vector<int64_t>(degree+1,0));
  for (int i = 0; i < num_vec; ++i){
    vector<polynomial> RLWE_mul1 = ct_multiplication(ct,ct_server[i], RelK, degree, ct_modulus,modulus, 
    logq, b, logb);
    lwe_ip[i] = extract_0(RLWE_mul1,degree);
    lwe_ip[i][degree] = modq_64(lwe_ip[i][degree]-scale_threshold, ct_modulus);
  }

  ofstream fout;
  fout.open("encrypted_result.txt");
  for(int j = 0 ; j < num_vec ; ++j){
    for(int i = 0 ; i < (degree+1) ; ++i){
      fout <<lwe_ip[j][i]<<" ";
    }
    fout <<endl;
  }
  

  fout.close();
  cout <<"Stored the encrypted result in encrypted_result.txt"<<endl;

}

int main(){
  eval_test();
}
  
