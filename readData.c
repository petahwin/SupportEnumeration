#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

extern int nActions1, nActions2;
extern float * payoffsA, * payoffsB;

/* Returns 0 on error, 1 on success */
/* In param: game file */
/* Populates constant valued global variables
 * including the payoff matrices and set of actions */
int readGame(char * gameData) {
  FILE * fp;
  char line[60];
  if ( (fp = fopen(gameData, "r")) ) {
//    int nActions1 = 1, nActions2 = 1;
    float payoff1 = 0., payoff2 = 0.;

    /* Get no. actions for player 1 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player1: %d", &nActions1);
    
    /* Get no. actions for player 2 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player2: %d", &nActions2);

    payoffsA = malloc(nActions1 * nActions2 * sizeof(float));
    payoffsB = malloc(nActions1 * nActions2 * sizeof(float));

    int i = 0;
    /* Get all of the payoffs */
    while (fgets(line, sizeof(line), fp)) {
      sscanf(line, "%f %f", &payoff1, &payoff2);
      payoffsA[i] = payoff1, payoffsB[i] = payoff2;
      ++i;
    }
    fclose(fp);
    return 1;
  } else {
    return 0;  
  }
}

