all:
	clang -o rs_decode rs_decode.c galois/galois.c -O3
	clang -o rs_encode rs_encode.c galois/galois.c -O3
