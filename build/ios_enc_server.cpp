// ================================================

//Code is designed for Facial verification in IOS

//Date: Feb 2025

// ================================================


#include <iostream>
#include <vector>
#include<unistd.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;



// Test set 


int scaler0 = 40;                         // scaler for input vector

int scaler1 = 40;                       // Scaler for Weights 

int num_vec = 1;                        // number of photos for each person

vector<vector<int>>  test_vector(num_vec);

// ======================================================  Neural Network  ==========================================================================

//read test vectors and labels
void read_vector() {
  ifstream infile;
  // plain weight path: size 512*2
  infile.open("test_vec.txt");
  //double sum_check = 0;
  for(int i = 0 ; i < num_vec ; ++i){
    for (int j = 0; j < 512; ++j) {
      double t;
      infile >> t;
      //sum_check += t*t;
      t = t * scaler0;
      int tt=(int)t;
      if (t > 0 && t - tt > 0.5) {
        test_vector[i].push_back( tt + 1);
      }
      else if (t < 0 && tt - t>0.5) {
        test_vector[i].push_back( tt - 1);
      }
      else {
        test_vector[i].push_back( tt);
      }
    }
  }
  
  infile.close();
  //cout <<"l2 - norm ^2 = " << sum_check << endl;
}

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


inline void modq_poly(vector<int64_t> &a, int n, int64_t q) {
  
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

inline void modq_poly_large(vector<int64_t> &a, int n, int64_t q) {
  
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
 void add_poly(vector<int64_t> &a, const vector<int64_t> & b, int n, int64_t q) {
  // INPUT: polynomials a and b, modulus q
  // OUTPUT: (a + b) mod q, stored in a
  for(int i = 0 ; i < n ; ++i){
    a[i] += b[i];
  }
  // mod q
  modq_poly_large(a,n, q);
}

 // Scaling
 void multi_scale_poly(int64_t t, vector<int64_t> &a, int n, int64_t q) {
  // INPUT: scaler t, polynomial a, modulus q
  // OUTPUT: t*a mod q, stored in a

  for (int i = 0; i < n; ++i) {
    a[i] *= t;
  }
  modq_poly_large(a,n, q);
}

 vector<int64_t> mul_poly( const vector<int64_t>& aa, 
  const vector<int64_t>& bb,int n, int64_t q) {
  vector<int64_t> c(n,0);
  for(int i = 0 ; i <n ; ++i){
    for(int j = 0 ; j < n ; ++j){
      if(i+j < n){
        c[i+j] += aa[i]*bb[j];
        modq(c[i+j],q);
      }
      else{
        c[i+j-n] -= aa[i]*bb[j];
        modq(c[i+j-n],q);
      }
    }
  }
  modq_poly_large(c,n,q);
  return c;
}

vector<vector<int64_t> >RLWE64_Enc_sk(const vector<int64_t>& m, 
  const vector<int64_t>  & s, int64_t p, int64_t q){
  int64_t alpha = q/p;
  int len = s.size();

  vector<vector<int64_t> >ct(2, vector<int64_t>(len,0));

  //generate random polynomial u
  vector<int64_t> u(len,0);
  for (int i = 0; i < len; ++i) {
    int64_t tempa = rand() % 2;
    u[i] = tempa;
  }

  //compute ct[0] = u
  ct[0] = u;


  //compute ct[1] = -as+alpha m+e

  vector<int64_t> as = mul_poly(s, u, len, q);
  multi_scale_poly(-1, as,len,q);

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
      as[i] += e;
      as[i] += (alpha*m[i]);
  }

  modq_poly_large(as,len,q);
  ct[1] = as;

  return ct;
}


vector<vector<int64_t> >RLWE64_Enc(const vector<int64_t>& m, 
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

vector<int64_t> RLWE64_Dec(const vector<vector<int64_t> >& ct, const vector<int64_t> & sk,
  int64_t p, int64_t q){
  int64_t alpha = q/p;
  int len = ct[0].size();

  vector<int64_t> pt = mul_poly(ct[0],sk,len,q);
  add_poly(pt, ct[1], len, q);

 // for (int i = 0; i < len; ++i){
 //   pt[i] = pt[i]/alpha;
 // }

  double tq = (double)p/(double)q;

  //(t/q)*as
  for(int i = 0 ; i < len ; ++i){
   // as[i] *= t;
    //double tq = (double)t/(double)q;
    double temp = tq*(double)pt[i] ;
    int64_t temp2 = (int64_t)temp;
    if(temp - (double)temp2 >= 0.5 && temp > 0){
      temp2 ++;
    }
    else if((double)temp2 - temp >= 0.5 && temp < 0){
        temp2 --;
    }
    pt[i] = temp2;
  //  if(as[i] < 0){
   //   as[i] += t;
   // }
  }


  return pt;


}

int main(){

  srand(time(0));
  read_vector();

  //set n,q,t for RLWE scheme
  const int n1 = 1024;//2048;                        //polynomial degree 
  const int64_t q1 = 3221225473;//206158430209;//2748779069441;//576460752154525697;           //ciphertext modulus
  const int p1 = 6000;//12289;                       //plaintext modulus

  //network parameter
  const int l1 = 512;

  //read lwe secret key
  ifstream fin;

  /*
  vector<int64_t> x(n1,0);
  fin.open("rlwe_sk.txt");
  for(int i=0;i<n1;++i){
    fin>>x[i];
  }
  fin.close();
  cout <<"Read RLWE secret key (for test). "<<endl;
  */

  //read lwe public key
  fin.open("rlwe_pk.txt");
  vector< vector<int64_t> > pk(2, vector<int64_t>(n1,0));
  for(int i=0 ; i < 2 ; ++i){
    for (int j = 0; j < n1; ++j){
      fin>> pk[i][j];
    }
  }
  fin.close();
  cout <<"Read RLWE public key. "<<endl;

/*
  //for test
  for (int i = 0; i < l1; ++i){
    cout <<test_vector_1[i]<<" ";
  }
  for (int i = 0; i < l1; ++i){
    cout <<test_vector_2[i]<<" ";
  }
  */

  vector<vector<int64_t>> input(num_vec, vector<int64_t>(n1,0));
  for (int i = 0; i < num_vec; ++i){
    input[i][0] = test_vector[i][0];
    for (int j = 1 ; j < l1 ; ++j){
      input[i][n1-j] = -1 * test_vector[i][j];
    }
  }

  ofstream fout;
  fout.open("test_enc_vec.txt");
  
  for (int k = 0; k < num_vec; ++k){
    vector<vector<int64_t> > ct1=RLWE64_Enc(input[k],pk,p1,q1); 
    for (int i = 0; i < 2; ++i){
      for(int j=0;j<n1;++j){
        fout<<ct1[i][j]<<" ";
      }
      fout << endl;
      //for test
      //cout <<"( "<<test_vector[i]<<", "<<Decrypt(q1,p1,n1,ct[i],x)<<") ";
    }
  }
  

  fout.close();
  cout<<"Encrypt the plain vectors and store the encrypted vectors. Number of vectors = "<<num_vec<<endl;


  return 0;

}


