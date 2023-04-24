#include "common.h"
#include "galois.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv) {

  if (argc != 3) {
    printf("Usage: \n%s -t \"text_to_encode\"\n%s -f "
           "\"file_to_encode\"\n",
           argv[0], argv[0]);
    exit(EXIT_FAILURE);
  }

  /*
   * Getting data from user input
   */

  char text[K] = {0};

  // check for flags in the arguments
  if (strcmp(argv[1], "-t") == 0) {
    // -t is the flag for text mode
    strcpy(text, argv[2]);
  } else if (strcmp(argv[1], "-f") == 0) {
    // -f is the flag for file mode
    FILE *input;
    input = fopen(argv[2], "rb");
    if (input == NULL) {
      printf("ERROR: Could not open the file: %s\n", argv[2]);
      return (EXIT_FAILURE);
    }
    fread(text, sizeof(text), 1, input);
    printf("File read successfully\n");
    fclose(input);

  } else {
    printf("Please provide a valid flag\n");
    exit(EXIT_FAILURE);
  }

  if (sizeof(text) / sizeof(char) > K) {
    printf(
        "The size of the text (=%ld) is greater than the length of the message "
        "block(K=223 bytes).\nPlease provide a shorter message.\n",
        sizeof(text) / sizeof(char));
    exit(EXIT_FAILURE);
  }

  /*
   * Generating the message
   */

  u_int8_t message[K] = {0};
  memcpy(message, text, sizeof(text));
  printf("Message to encode: %s\n", text);

  /*
   * Encoding
   */

  // define the codeword
  u_int8_t codeword[N] = {0};

  // encoding using Horner polynomial evaluation
  for (int i = 0; i < N; i++) {
    codeword[i] = message[K - 1];
    for (int j = K - 2; j >= 0; j--) {
      codeword[i] =
          galois_single_multiply(codeword[i], eval_points[i], W) ^ message[j];
    }
  }

  /*
   * Exporting to a file
   */

  FILE *output;
  char filename[30];
  strcat(filename, "enc_");
  sprintf(filename + 4, "%lu", (unsigned long)time(NULL));
  strcat(filename, ".rs");
  output = fopen(filename, "wb");

  if (output == NULL) {
    printf("ERROR: Could not open the file:\n");
    return (EXIT_FAILURE);
  }

  fwrite(codeword, sizeof(codeword), 1, output);
  printf("Message encoded successfully to %s\n", filename);

  fclose(output);

  return 0;
}