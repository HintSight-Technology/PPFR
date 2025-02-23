# N2HE-FacialVerification-Demo

## Prerequisites
- [OpenSSL](https://www.openssl.org/)  3.2.1 or later
- [hexl](https://github.com/intel/hexl) v1.2.5 or later

## Installation
Supported platforms: Linux / macOS with Intel CPU.  

```
cd build
cmake ..
make
```

## Usage
```
./init       // generate keys
./enc        // encrypt the features
./enc_server // encrypt the model
./eval       // homomorphically inference the encrypted model, using the encrypted features
./dec        // decrypt the result
```
