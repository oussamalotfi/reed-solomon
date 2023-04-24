#include "common.h"
#include "galois.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



void gf_guass_RowM_ColN(unsigned char coef[], int coef_row, int coef_col,
                        unsigned char b[], unsigned char x[]) {
  int i, j, k;
  unsigned char mi;
  unsigned char mx, tmp;
  unsigned char a[N][N] = {0};

  for (i = 0; i < coef_row; i++) {
    for (j = 0; j < coef_col; j++)
      a[i][j] = coef[i * coef_col + j];
  }

  for (i = 0; i < coef_col - 1; i++) {
    // 找主元素
    for (j = i + 1, mi = i, mx = a[i][i]; j < coef_row; j++) {
      if (a[j][i] > mx) {
        mi = j;
        mx = a[j][i];
      }
    }

    // 交换两行
    if (i < mi) {
      tmp = b[i];
      b[i] = b[mi];
      b[mi] = tmp;
      for (j = i; j < coef_col; j++) {
        tmp = a[i][j];
        a[i][j] = a[mi][j];
        a[mi][j] = tmp;
      }
    }

    // 高斯消元
    for (j = i + 1; j < coef_row; j++) {
      if (a[i][i] == 0)
        // printf("a[i=%d][i=%d] == %d\n", i, i, a[i][i]);
        ;
      else {
        tmp = galois_single_divide(a[j][i], a[i][i], W);
        b[j] = (b[j] ^ galois_single_multiply(b[i], tmp, W));
        for (k = i; k < coef_col; k++)
          a[j][k] = a[j][k] ^ galois_single_multiply(a[i][k], tmp, W);
      }
    }
  }
  // 求解方程
  // 判断是否为0，即是否能偶求解
  if (a[coef_col - 1][coef_col - 1] == 0) {
    // printf("a[coef_col - 1=%d][coef_col - 1=%d] == %d\n", i, i, a[coef_col -
    // 1][coef_col - 1]);
    ;
  } else {
    x[coef_col - 1] =
        galois_single_divide(b[coef_col - 1], a[coef_col - 1][coef_col - 1], W);
    for (i = coef_col - 2; i >= 0; i--) {
      x[i] = b[i];
      for (j = i + 1; j < coef_col; j++)
        x[i] = x[i] ^ galois_single_multiply(a[i][j], x[j], W);
      x[i] = galois_single_divide(x[i], a[i][i], W);
    }
    // printf("success\n");
  }
}

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("Usage:\n%s \"file_to_decode\"\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /*
   * Received Codeword
   */

  unsigned char received_codeword[N] = {0};

  FILE *input;
  input = fopen(argv[1], "rb");
  if (input == NULL) {
    printf("ERROR: Could not open the file: %s\n", argv[1]);
    return (EXIT_FAILURE);
  }
  fread(received_codeword, sizeof(received_codeword), 1, input);
  printf("File read successfully\n");
  fclose(input);

  /*
   * Generating the matrix A and the vector b
   */
  int A[N][N] = {0};
  int b[N] = {0};

  for (int i = 0; i < N; i++) {
    A[i][0] = received_codeword[i];
    for (int j = 1; j < e; j++) {
      // evaluating the first e-1 columns of the matrix A
      A[i][j] = galois_single_multiply(A[i][j - 1], eval_points[i], W);
    }

    A[i][e] = 1;
    for (int j = e + 1; j < N; j++) {
      // evaluating the rest (e to N columns) of the columns of the matrix A
      A[i][j] = galois_single_multiply(A[i][j - 1], eval_points[i], W);
    }

    // evaluating the vector B
    b[i] = galois_single_multiply(A[i][e - 1], eval_points[i], W);
  }

  // /*
  //  * Generating the Augmented Matrix  augA = A|b
  //  */

  // int augA[N][N + 1];
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N + 1; j++) {
  //     if (j < N)
  //       augA[i][j] = A[i][j];
  //     else
  //       augA[i][j] = b[i];
  //   }
  // }

  // unsigned int c;
  // unsigned int x[N];
  // for (int j = 0; j < N; j++) {
  //   for (int i = 0; i < N; i++) {
  //     // printf("%d,%d\n",j,i);
  //     if (i != j) {
  //       // c = galois_single_divide(augA[i][j], augA[j][j], W);
  //       c = galois_single_multiply(augA[i][j], galois_inverse(augA[j][j], W),
  //                                  W);
  //       for (int k = 0; k < N + 1; k++) {
  //         augA[i][k] = augA[i][k] ^ galois_single_multiply(c, augA[j][k], W);
  //       }
  //     }
  //   }
  // }
  //   for (int i = 0; i < N; i++) {
  //     // c = galois_single_divide(augA[i][N], augA[i][i], W);

  //   x[i] = galois_single_multiply(augA[i][N], galois_inverse(augA[i][i], W),
  //   W);
  // }



  /*
   * Solving the linear system Ax = B using
   * Gaussian elimination
   */

  unsigned char copyx[N] = {0};
  unsigned char copyA[N][N], copyb[N];

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      copyA[i][j] = A[i][j];
    }
    copyb[i] = b[i];
  }

  gf_guass_RowM_ColN(*copyA, N, N, copyb, copyx);
  int x[N];
  for (int i = 0; i < N; i++) {
    x[i] = copyx[i];
  }

  // testing the solution
  int sum = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      sum ^= galois_single_multiply(copyA[i][j], x[j], W);
    }
    if (sum != b[i]) {
      printf("The sum doesn't match, there was an error while solving Ax=b\n");
      exit(EXIT_FAILURE);
    }
    sum = 0;
  }

  /*
   * Defining E and Q
   */

  int E[e + 1] = {0};
  E[e] = 1;
  int Q[N - e] = {0};

  memcpy(E, x, sizeof(int) * e);
  memcpy(Q, x + e, sizeof(int) * (N - e));



  /*
   * Dividing Q by E
   * The quotient is stored in total_coef
   */

  int nom_deg = N - e - 1, den_deg = e;

  int nom_coef[N - e - 1 + 1];
  memcpy(nom_coef, Q, sizeof(Q));
  int den_coef[e + 1];
  memcpy(den_coef, E, sizeof(E));

  // doing long division
  int total_coef[nom_deg - den_deg + 1];
  int j;
  int temp_coef[e + 1];

  for (int i = 0; i < nom_deg - den_deg + 1; i++) {
    // divide
    total_coef[i] =
        galois_single_divide(nom_coef[nom_deg - i], den_coef[den_deg], W);
    // multiply
    for (j = 0; j < den_deg + 1; j++) {
      temp_coef[den_deg - j] =
          galois_single_multiply(total_coef[i], den_coef[den_deg - j], W);
    }
    // subtract
    for (j = 0; j < den_deg + 1; j++) {
      nom_coef[nom_deg - j - i] =
          (nom_coef[nom_deg - j - i] ^ temp_coef[den_deg - j]);
    }
  }



  // /*
  //  * The following commented lines is to print the result
  //  * of the Q/E division
  //  */

  // printf("\nthe solution is:\n");
  // for (int i = 0; i < nom_deg - den_deg + 1; i++) {
  //   printf("(%d)X^%d+\n", total_coef[i], nom_deg - den_deg - i);
  // }

  // dealing with remainder
  //------------------------
  // printf("[");
  // for (int i = 0; i < nom_deg + 1; i++) {

  //   if (nom_coef[nom_deg - i] != 0.0 && nom_deg - i != 0) {
  //     printf("(%d)X^%d+", nom_coef[nom_deg - i], nom_deg - i);
  //   } else if (nom_coef[nom_deg - i] != 0.0 && nom_deg - i == 0) {
  //     printf("(%0d)X^%d", nom_coef[nom_deg - i], nom_deg - i);
  //     printf("]/[");
  //     for (j = 0; j < den_deg + 1; j++) {
  //       printf("(%d)X^%d+", den_coef[den_deg - j], den_deg - j);
  //     }
  //   }
  // }

  // printf("0]\n");



  /*
   * Print Decoded Message
   */

  char decoded_message[K];
  int zero_index = 0;
  printf("Decoded message:\n");
  printf("\n-----------------\n");
  for (int i = 0; i <= K; i++) {
    decoded_message[i] = total_coef[K - 1 - i];
    if (decoded_message[i] == 0) {
      zero_index = i;
      break;
    }
  }
  for (int i = 0; i < K; i++) {
    printf("%d\n", decoded_message[i]);
  }
  printf("zero index is: %d\n", zero_index);

  printf("\n-----------------\n");



  /*
   * Export to a file
   */

  FILE *output;
  char filename[30];
  strcat(filename, "dec_");
  sprintf(filename + 4, "%lu", (unsigned long)time(NULL));
  strcat(filename, ".txt");
  output = fopen(filename, "wb");

  if (output == NULL) {
    printf("ERROR: Could not open the file:\n");
    return (EXIT_FAILURE);
  }

  fwrite(decoded_message, zero_index * sizeof(char), 1, output);
  printf("Message decoded successfully to %s\n", filename);
  fclose(output);

  return 0;
}