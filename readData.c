#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Returns 0 on error, 1 on success */
/* In param: game file */
/* Out param: nothing for now */
/* Side effects: prints out the information */
int readGame(char * gameData) {
  FILE * fp;
  char line[60];
  if (fp = fopen(gameData, "r")) {
    int nActions1 = 1, nActions2 = 1;
    float payoff1 = 0., payoff2 = 0.;

    /* Get no. actions for player 1 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player1: %d", &nActions1);
    
    /* Get no. actions for player 2 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player2: %d", &nActions2);

    printf("P1: %d, P2: %d\n", nActions1, nActions2);

    /* Get all of the payoffs */
    while (fgets(line, sizeof(line), fp)) {
      sscanf(line, "%f %f", &payoff1, &payoff2);
      printf("%f, %f\n", payoff1, payoff2);
    }
    return 1;
  } else {
    return 0;  
  }
}

int main(int argc, char * argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Incorrect command line args");
    return 1;
  } else {
    readGame(argv[1]);
    return 0;
  }
}

