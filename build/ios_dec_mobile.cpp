// ================================================

//Code is designed for Facial verification in IOS

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
//#include <immintrin.h>
#include <string>


using namespace std;

int scaler0 = 40;                         // scaler for input vector

int scaler1 = 40;                       // Scaler for Weights 




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


// LWE Decryption
int64_t Decrypt(int64_t q, int64_t p, int n, const vector<int64_t>& c, const vector<int>& x) {
  // INPUT: modulus q, dimension n, LWE ciphertext c = (vec(a),b), LWE key x
  // OUTPUT: Decryption (b + <a,x>) mod q
  int64_t alpha = q/p;
  int64_t ip2 = c[n];
  for (int i = 0; i < n; ++i) {
    ip2 += (c[i]*(int64_t)x[i]);
      modq(ip2, q);
  }

  while(ip2 < 0){
      ip2 += q;
    }
    int64_t temp2 = (ip2 + alpha/2) % q;
    temp2 /= alpha;
    if (temp2 > p/2){
      temp2 -= p;
    }
    return temp2;

}


string ios_enc_vec_dec(){

  // ======================== Initialization and Information Output ==========================

  srand(time(0));
  //read_vector();
  int num_vec = 1;                        // number of photos for each person

  //set n,q,t for RLWE scheme
  const int n1 = 1024;//2048;                        //polynomial degree 
  const int64_t q1 = 3221225473;//206158430209;//2748779069441;//576460752154525697;           //ciphertext modulus
  const int p1 = 6000;//12289;                       //plaintext modulus

  //network parameter
  const int l1 = 512;


  //read LWE secret key
  ifstream fin;
  fin.open("rlwe_sk.txt");
  vector<int> x(n1);
  for (int i=0;i<n1;++i){
    fin >> x[i];
  }
  fin.close();
  cout <<"read LWE secret key."<<endl;

  //read encrypted result
  vector<vector<int64_t>> ct_ip(num_vec, vector<int64_t>(n1+1,0));

  fin.open("encrypted_result.txt");
  if(!fin.is_open()){
    cout <<"Cannot open file encrypted_result.txt"<<endl;
  }
  for(int k = 0 ; k < num_vec ; ++k){
    for(int i = 0 ; i < (n1+1) ; ++i){
      fin >>ct_ip[k][i];
    }
  }
  
  fin.close();
  cout <<"Read encrypted result from encrypted_result.txt."<<endl;

  int num_yes = 0;
  int num_no = 0;

  for (int i = 0; i < num_vec; ++i){
    int64_t dec_ip1 = Decrypt(q1,p1,n1,ct_ip[i],x);
    cout <<"Decrypt result = "<<dec_ip1<<endl;
    if (dec_ip1 >= 0){
      num_yes ++;
    }
    else{
      num_no++;
    }
  }
  
  cout <<"verification result: ";
  if(num_yes >= num_no){
    string a = "yes";
    cout <<a <<endl;
    return a;
  }
  else{
    string a = "no";
    cout <<a <<endl;
    return a;
  }
    


}


int main(){
  string b = ios_enc_vec_dec();
}





