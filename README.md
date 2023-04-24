# Reed-Solomon RS(n=255,k=223) encoding and decoding using Berlekamp-Welch algorithm


## Encoding
The encoding process is simple, it works by puting the original information as coefficients of a polynomial of degree k-1 and then evaluating it at n distinct points.



## Decoding
For more informations on how the Berlekamp-Welch decoding algorithm works please refer to the wikipedia page.

> [Berlekamp-Welch decoding algorithm](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Welch_algorithm)


## Usage

### Compiling

```console
$ clang -o rs_encode rs_encode.c galois/galois.c -O3
$ clang -o rs_decode rs_decode.c galois/galois.c -O3
```

### Encoding

To encode text directly
```console
$ ./rs_encode -t "text_to_encode"
```

To encode a text file
```console
$ ./rs_encode -f "file_to_encode"
```


### Decoding

To decode an encoded file
```console
$ ./rs_decode "file_to_decode"
```


## Limitations

For now the encoding program can only encode a bloc of $k = 223$ bytes. 

## Credit
- [Fast Galois Field Arithmetic Library in C/C++ ](http://web.eecs.utk.edu/~jplank/plank/papers/CS-07-593/) by [James S. Plank](http://web.eecs.utk.edu/~jplank/)

- [Gaussian elimination algorithm in Galois Field](https://github.com/edisonyangyang/gf) by [edisonyangyang](edisonyangyang)

- [Polynomial-Division](https://github.com/Remonhanyz/Polynomial-Division) by [Remonhanyz](https://github.com/Remonhanyz)