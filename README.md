# Pivot-tangent Method
This repo contains the code for pivot-tangent method - A method to polynomially approximate inverse sqrt function for application in homomorphic encryption scheme like CKKS. The file `pivot_tangent_unencrypted.py` can be used to find different parameters under different setttings in the pivot-tangent method. The folder `inv_sqr_comp` contains the implementation of pivot-tangent method in CKKS homomorphic scheme using the [Microsoft SEAL library](https://github.com/microsoft/SEAL). The file `inv_sqrt_comp.cpp` contains the comparision of pivot-tangent method and constrained linear regression approach used in the [HPCA algorithm](https://github.com/pandasamanvaya/Homomorphic_PCA). To run the example, clone the repo and then:-
```
cd inv_sqr_comp
cmake .
make
./comp
```
