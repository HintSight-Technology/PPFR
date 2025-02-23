// ================================================

//Code is designed for Facial verification 

//Date: Feb 2025

// ================================================


#include <iostream>
#include <vector>
#include <ctime>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <condition_variable>
#include <chrono>
#include <thread>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <immintrin.h>
#include <string>
#include "include.hpp"

using namespace std;

// ======================================================  LWE Encryption Scheme  ==========================================================================

//mod q
inline void modq(int64_t number, int64_t q) {
  if(number < 0){
      int64_t temp = (-1*number)/q + 1;
      number += temp*q;
    }
    if(number >= q){
      int64_t temp = number/q;
      number -= temp*q;
    }
    if(number >= q/2){
      number -= q;
    }
  return ;
}

inline void modq_poly_t(vector<int64_t> &a, int n, int64_t q) {
  
  // INPUT: polynomial a, modulus q
  // OUTPUT: (modified) a mod q, stored in a

  for (int i = 0; i < n; ++i) {

    while (a[i]<0) {
      a[i] += q;
    }

    while (a[i]>=q) {
      a[i] -= q;
    }

    if(a[i] > (q-1)/2){
      a[i] -=q;
    }

  }
  
}

inline void modq_poly_large_t(vector<int64_t> &a, int n, int64_t q) {
  
  // INPUT: polynomial a, modulus q
  // OUTPUT: (modified) a mod q, stored in a

  for (int i = 0; i < n; ++i) {

    if(a[i]<0){
      int64_t temp = -1 * a[i];
      a[i] += q*(temp/q+1);
    }

    if(a[i] >= q){
      a[i] -= q*(a[i]/q);
    }

    if(a[i] > (q-1)/2){
      a[i] -=q;
    }

  }
  
}

// Addition between polynomials
  void add_poly_t(vector<int64_t> &a, const vector<int64_t> & b, int n, int64_t q) {
  // INPUT: polynomials a and b, modulus q
  // OUTPUT: (a + b) mod q, stored in a
  for(int i = 0 ; i < n ; ++i){
    a[i] += b[i];
  }
  // mod q
  modq_poly_t(a,n, q);
}

 // Scaling
  void multi_scale_poly_t(int64_t t, vector<int64_t> &a, int n, int64_t q) {
  // INPUT: scaler t, polynomial a, modulus q
  // OUTPUT: t*a mod q, stored in a

  for (int i = 0; i < n; ++i) {
    a[i] *= t;
  }
  modq_poly_large_t(a,n, q);
}

 vector<int64_t> mul_poly_t( const vector<int64_t>& aa, 
  const vector<int64_t>& bb,int n, int64_t q) {
  vector<int64_t> c(n,0);
  for(int i = 0 ; i <n ; ++i){
    for(int j = 0 ; j < n ; ++j){
      if(i+j < n){
        c[i+j] += aa[i]*bb[j];
      }
      else{
        c[i+j-n] -= aa[i]*bb[j];
      }
    }
  }
  modq_poly_large_t(c,n,q);
  return c;
}


vector<int64_t> RLWE64_SKGen(int n) {
  vector<int64_t> x;
  for (int i = 0; i < n; ++i) {
    int64_t tempx = rand() % 3;
    x.push_back(tempx-1);
  }
  return x;
}

vector< vector<int64_t> > RLWE64_PKGen(const vector<int64_t>& sk, int64_t q){
  int len = sk.size();
  polynomial m(len,0);
  vector<vector<int64_t>> ct = RLWE64_Enc(len, q, q, 3, m, sk);

  /*

  int len = sk.size();
  vector<vector<int64_t> >ct(2, vector<int64_t>(len,0));

  //Generate random polynomial a
  for (int i = 0; i < len; ++i) {
    int64_t tempa = rand() % q-q/2;
    ct[0][i] = tempa;
  }

  //compute -(as)
  vector<int64_t> as = mul_poly_t(sk, ct[0], len, q);
  multi_scale_poly_t(-1, as,len,q);

  //Generate e and compute -as+e
  int non_zero_count = 0;
  for (int i = 0; i < len; ++i) {
    int64_t indexe=rand()%16;
      int64_t e=0;
      if(indexe == 0){
        e=1;
        non_zero_count++;
      }
      else if(indexe == 1){
        e=-1;
        non_zero_count++;
      }
      as[i] += e;
  }
  //cout <<"non zero ei = "<<non_zero_count<<endl;

  modq_poly_large_t(as,len,q);
  ct[1] = as;
  */

  return ct;



}

vector<vector<int64_t> >RLWE64_Enc_pk(const vector<int64_t>& m, 
  const vector<vector<int64_t> > & pk, int64_t p, int64_t q){
  int64_t alpha = q/p;
  int len = pk[0].size();

  vector<vector<int64_t> >ct(2, vector<int64_t>(len,0));

  //generate random polynomial u
  vector<int64_t> u(len,0);
  for (int i = 0; i < len; ++i) {
    int64_t tempa = rand() % 2;
    u[i] = tempa;
  }

  //compute ct[0] = pk[0]u+e1

  vector<int64_t> pk0u = mul_poly(pk[0], u, len, q);

  //Generate e and compute -as+e
  for (int i = 0; i < len; ++i) {
    int64_t indexe=rand()%16;
      int64_t e=0;
      if(indexe == 0){
        e=1;
      }
      else if(indexe == 1){
        e=-1;
      }
      pk0u[i] += e;
  }

  modq_poly_large(pk0u,len,q);
  ct[0] = pk0u;

  //compute ct[1] = pk[1]u+e1+alpham

  vector<int64_t> pk1u = mul_poly(pk[1], u, len, q);

  //Generate e and compute -as+e
  for (int i = 0; i < len; ++i) {
    int64_t indexe=rand()%16;
      int64_t e=0;
      if(indexe == 0){
        e=1;
      }
      else if(indexe == 1){
        e=-1;
      }
      pk1u[i] += e;
      pk1u[i] += (alpha*m[i]);
  }

  modq_poly_large(pk1u,len,q);
  ct[1] = pk1u;

  return ct;
}

vector<vector<polynomial>> extRLWE_pk(int n, int64_t q, int64_t t, int k, double var, const polynomial & m, 
  const vector<polynomial> & pk, int64_t b, int logb) {

  vector<vector<polynomial>> extRLWEct;

  //~RLWE(m,2^i)
  int64_t temp2k = b;
  polynomial mi = copy(m);
  for (int i = 0; i < (k/logb); ++i) {
   // cout <<i <<endl;
    if (i > 0) {
      multi_scale_poly(temp2k, mi,n, q);
    }
    vector<polynomial> c1i = RLWE64_Enc_pk(mi, pk, t, q);
   // cout <<"enc "<<c1i.size()<<endl;
    extRLWEct.push_back(c1i);
  }

  return extRLWEct;
}

extern vector<vector<polynomial>> RelK_Gen_pk(const vector<polynomial> & pk, const polynomial & s, int N, int64_t ct_modulus, int64_t pt_modulus, int logq, 
  double var, int64_t b, int logb){
  polynomial s2 = NTTMul(s,s,N,pt_modulus);
  print_polynomial(s2);
  return extRLWE_pk(N,ct_modulus,pt_modulus,logq,var,s2,pk,b,logb);
}

int main(){
  srand(time(0));


  //set n,q,t for RLWE scheme
  const int n1 = 1024;//2048;                        //polynomial degree 
  const int64_t q1 = 3221225473;//206158430209;//2748779069441;//576460752154525697;           //ciphertext modulus
  const int p1 = 6000;//12289;                       //plaintext modulus

  int logq = 32;//38;//42;//60;
  double var = 3.2;
  int64_t b = 2;
  int logb = 1;
  int64_t modulus = 6000;//12289; 
  int64_t primitive_root = 3;

  vector < int64_t > x = RLWE64_SKGen(n1);
  ofstream fout;
  fout.open("rlwe_sk.txt");
  for(int i=0;i<n1;++i){
    fout<<x[i]<<" ";
  }
  fout.close();
  cout<<"Generate and store RLWE secret key."<<endl;


  vector<vector<int64_t> > pk = RLWE64_PKGen(x, q1);
  fout.open("rlwe_pk.txt");
  for(int i=0 ; i < 2 ; ++i){
    for (int j = 0; j < n1; ++j){
      fout << pk[i][j]<<" ";
    }
    fout <<endl;
  }
  fout.close();
  cout<<"Generate and store RLWE public key."<<endl;


  vector<vector<polynomial>> RelK=RelK_Gen(x,n1,q1,q1,logq,var,b,logb);


  fout.open("rlwe_relk.txt");
  for(int k = 0 ; k < logq ; ++k){
    for(int i=0 ; i < 2 ; ++i){
      for (int j = 0; j < n1; ++j){
        fout << RelK[k][i][j]<<" ";
      }
    fout <<endl;
    }
  }
  
  fout.close();

  cout <<"Generate and store evaluation key. "<<endl;
  cout <<endl;



  return 0;

}


