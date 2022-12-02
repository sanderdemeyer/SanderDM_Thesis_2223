# MatlabTrack

## Installation
### Prerequisites
In order to get started, you'll need the following:
- MATLAB (> R2020b)
- Parallel Computing Toolbox
- C compiler

You can then extract the code into any folder, making sure to add it to the MATLAB path. This is done either manually (right click folder -> Add to Path -> Selected Folders and Subfolders), or by navigating to the MatlabTrack directory and running ```doPath```.
Additionally, for the initial installation you should run the command ```GetMD5``` once, which will compile the required mex files.

### Settings
Note that usage of this package will create a folder ```LocalStorage``` in the current working directory, in order to store cached data. The default location of this folder can be changed in ```Mediators/@Mediator/Mediator.m```, line 61.

## Getting Started
You can now start using the MatlabTrack. The cornerstone of this package are the classes ```@TensorNone```, ```@TensorSymm``` and ```@Tensor6j```, the first defining non-symmetric tensors, while the latter are two separate implementations of symmetric tensors.

### Creating tensors
Creating non-symmetric tensors is done in the following way:
```matlabsession
>> dims = [2 3 2];
>> t = TensorNone.Random(dims)

t = 

  TensorNone with properties:

     var: [2×3×2 double]
    dims: [2 3 2]
    legs: 3
>> t.var

ans(:,:,1) =

  -0.1924 - 0.8045i  -0.7648 + 0.8351i  -1.4224 + 0.2157i
   0.8886 + 0.6966i  -1.4023 - 0.2437i   0.4882 - 1.1658i


ans(:,:,2) =

  -0.1774 - 1.1480i   1.4193 + 0.7223i   0.1978 - 0.6669i
  -0.1961 + 0.1049i   0.2916 + 2.5855i   1.5877 + 0.1873i
```
This creates a three-leg tensor of dimension 2 x 3 x 2, with normally distributed pseudorandom complex numbers. For more control over the initialisation, see also ```TensorNone.New```, ```TensorNone.Ones``` and ```TensorNone.Zeros```.

These objects behave like tensors, and as such support the following (incomplete) list of operations:
### Permutations
```matlabsession
>> t2 = Permute(t, [2 1 3])

t2 = 

  TensorNone with properties:

     var: [3×2×2 double]
    dims: [3 2 2]
    legs: 3
>> t2.var

ans(:,:,1) =

  -0.1924 - 0.8045i   0.8886 + 0.6966i
  -0.7648 + 0.8351i  -1.4023 - 0.2437i
  -1.4224 + 0.2157i   0.4882 - 1.1658i


ans(:,:,2) =

  -0.1774 - 1.1480i  -0.1961 + 0.1049i
   1.4193 + 0.7223i   0.2916 + 2.5855i
   0.1978 - 0.6669i   1.5877 + 0.1873i
```
### Contractions
The contraction routine uses ncon-like syntax, with arbitrary amount of tensors:
```matlabsession
>> t1 = TensorNone.Random([2 3 2]);
>> t2 = TensorNone.Random([3 4]);
>> t3 = Contract({t1 t2}, {[-1 1 -2] [1 -3]})

t3 = 

  TensorNone with properties:

     var: [2×2×4 double]
    dims: [2 2 4]
    legs: 3
>> t4 = Contract({t1 t2 t3}, {[1 2 3] [2 -1] [1 3 -2]})

t4 = 

  TensorNone with properties:

     var: [4×4 double]
    dims: [4 4]
    legs: 2
```
### Factorisations
