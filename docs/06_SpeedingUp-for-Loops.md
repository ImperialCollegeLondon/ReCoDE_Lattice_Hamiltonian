# Speeding up ```for``` loops: the ```get_rate_mat``` function

The purpose of the ```get_rate_mat``` function (lines 477-523 of ```lattice.py```) is to calculate the elements of the rate matrix which describes the rate of population transfer between the excited states of the lattice. Each element of the rate matrix takes the form 
```math
k_{\alpha \beta}(\omega_{\alpha \beta}) = \left[ \sum_{k}|c_{kk}^{(\alpha)}|^{2}|c_{kk}^{(\beta)}|^{2} + \sum_{i\neq j}|c_{ij}^{(\alpha)}|^{2}|c_{ij}^{(\beta)}|^{2} + \sum_{j\neq i}|c_{ij}^{(\alpha)}|^{2}|c_{ij}^{(\beta)}|^{2} \right] C(\omega_{\alpha \beta})
```
as is described in greater detail [here](02_FindingSteadyStatePopulations.md). Looking at this summation, we see that each term in the rate matrix requires a summation over all the basis states of the system (basis states are indexed using the latin alphabet, and eigenstates using the greek alphabet). For a $N\times N$ lattice, the number of basis states is $N^{4}$. In addition to this, the rate matrix itself contains $N^{2}$ elements as we calculate rates between every possible pair of eigenstates. Thus, the calculation of the rate matrix scales as $N^{8}$ meaning that it limits the speed of the code for all but the smallest lattice sizes. This means that we need this function to run as quickly as possible in order to minimise the overall runtime of the code and, in this file, we will explore some of the techniques we have used to achieve this. 

## 1) Vectorize functions

Vectorizing functions is always a good first step to speed up python code. Simply put, a vectorized function is one which can accept an iterable as its input parameter, rather than only accepting a single element at a time from an interable. For example, the fact that ```sum``` is a vectorized function means that you can replace the code
```
total = 0
for i in range(1000):
  total += i
```
with 
```
total = sum(range(1000))
```
In general, it is always good to try and use pre-exsisting functions from libraries such as numpy or scipy as these functions are typically vectorized and also use C to do much of the heavy lifting, making them much faster than equivalent functions written purely in Python. For example, in lines 515-516, we use the ```np.matmul``` function to perform the summation over all of the lattice's basis states.

However, in some cases, there will not be a pre-existing function which does what you need and then it is best to write your own, vectorized function. We have done this to calculate the values of the correlation function ($C(\omega_{\alpha \beta})$ in the expression for $k_{\alpha \beta}(\omega_{\alpha \beta})$ above). Our functions are saved in the ```redfield.py``` file and we call them in line 509 of ```get_rate_mat``` via the ```correlation_function_real_part``` function. However, to take advantage of the functions in ```redfield.py```, we first need to make a numpy array containing the difference in energy between each possible pair of eigenstates (i.e., the $\omega_{\alpha \beta}$) which we do using...

## 2) List Comprehension

To make a numpy array containing the values of $\omega_{\alpha \beta}$, we combine the numpy function ```np.fromiter``` (see the documentation [here](https://numpy.org/doc/stable/reference/generated/numpy.fromiter.html)) with a generator written in the form of a *list comprehension*. A list comprehension is a way of rewriting ```for``` loops which is particularly useful if the logic involved in each iteration of the for loop is simple enough to be encapsulated in a single line. For example, the code:
```
a = [1,2,3,4,5]
b = []
for i in range(len(a)):
  b.append(a[i]**2)
```
can be rewritten as:
```
a = [1,2,3,4,5]
b = [i**2 for i in a]
```

## 3) Convert nested for loops into a single for loop

Nested ```for``` loops are ones which involve looping over two (or more) indicies e.g.:
```
for i in range(10):
  for j in range(10):
    ...
```
In many situations, it is possible to convert a nested ```for``` loop into a single ```for``` loop by combining the indices into a single lcontaining pairs of indices using a package like [itertools](https://docs.python.org/3/library/itertools.html). 

In the case of ```get_rate_matrix```, the indices needed are those which label the row and column of the element of the rate matrix being calculated (i.e., the values of $\alpha$ and $\beta$ in $k_{\alpha \beta}$). Fortunately, numpy has a built in function (```np.tril_indices_from```, see the [documentation](https://numpy.org/doc/stable/reference/generated/numpy.tril_indices_from.html)) which can generate an iterable of the indices describing the elements making up a lower (or upper) triangular matrix. We use this to remove the double for loop which we would otherwise need to iterate over all the possible pairs of eigenstates. 

