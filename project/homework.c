#include "fem.h"

double *femElasticitySolve(femProblem *theProblem) {

  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  renumberMesh(theGeometry);
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;

  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int map[4], mapX[4], mapY[4];

  int nLocal = theMesh->nLocalNode;

  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double g = theProblem->g;
  double **A = theSystem->A;
  double *B = theSystem->B;

  //
  //  A faire :-)
  //
  int *theConstrainedNodes = theProblem->constrainedNodes;
  for (int i = 0; i < theSystem->size; i++) {
    if (theConstrainedNodes[i] != -1) {
      double value = theProblem->conditions[theConstrainedNodes[i]]->value;
      femFullSystemConstrain(theSystem, i, value);
    }
  }

  return femFullSystemEliminate(theSystem);
}


void renumberMesh(femGeo *theGeometry) {
  if (theGeometry == NULL || theGeometry->theNodes == NULL) {
      Error("theGeometry ou theNodes non initialisé");
  }

  int nNodes = theGeometry->theNodes->nNodes;

  /* --- 1. Construction de la liste d’adjacence --- */
  // On utilise la connectivité des éléments si disponible, sinon celle des arêtes.
  femMesh *mesh = (theGeometry->theElements) ? theGeometry->theElements : theGeometry->theEdges;
  if (!mesh) {
      Error("Aucune connectivité (éléments/ arêtes) disponible pour la renumérotation RCM.");
  }
  int nElem = mesh->nElem;
  int nLocal = mesh->nLocalNode;
  int *elem = mesh->elem;

  // Allocation de la structure temporaire pour chaque nœud
  typedef struct {
      int count;
      int capacity;
      int *neighbors;
  } NodeNeighbors;
  NodeNeighbors *adj = malloc(nNodes * sizeof(NodeNeighbors));
  if (adj == NULL) {
      Error("Allocation mémoire échouée pour 'adj'");
  }
  for (int i = 0; i < nNodes; i++) {
      adj[i].count = 0;
      adj[i].capacity = 4; // capacité initiale
      adj[i].neighbors = malloc(adj[i].capacity * sizeof(int));
      if (adj[i].neighbors == NULL) {
          Error("Allocation mémoire échouée pour 'adj[i].neighbors'");
      }
  }

  // Pour chaque élément, ajouter une arête entre chaque paire de nœuds
  for (int e = 0; e < nElem; e++) {
      for (int i = 0; i < nLocal; i++) {
          for (int j = i + 1; j < nLocal; j++) {
              int ni = elem[e * nLocal + i];
              int nj = elem[e * nLocal + j];
              // Ajout de nj comme voisin de ni si non déjà présent
              int exists = 0;
              for (int k = 0; k < adj[ni].count; k++) {
                  if (adj[ni].neighbors[k] == nj) { exists = 1; break; }
              }
              if (!exists) {
                  if (adj[ni].count == adj[ni].capacity) {
                      adj[ni].capacity *= 2;
                      adj[ni].neighbors = realloc(adj[ni].neighbors, adj[ni].capacity * sizeof(int));
                  }
                  adj[ni].neighbors[adj[ni].count++] = nj;
              }
              // Ajout de ni comme voisin de nj
              exists = 0;
              for (int k = 0; k < adj[nj].count; k++) {
                  if (adj[nj].neighbors[k] == ni) { exists = 1; break; }
              }
              if (!exists) {
                  if (adj[nj].count == adj[nj].capacity) {
                      adj[nj].capacity *= 2;
                      adj[nj].neighbors = realloc(adj[nj].neighbors, adj[nj].capacity * sizeof(int));
                  }
                  adj[nj].neighbors[adj[nj].count++] = ni;
              }
          }
      }
  }

  /* --- 2. Calcul de l’ordre RCM --- */
  // Allocation d'un tableau pour l'ordre RCM et d'un tableau de marqueurs de visite.
  int *rcm_order = malloc(nNodes * sizeof(int));
  int *visited = calloc(nNodes, sizeof(int));
  if (rcm_order == NULL || visited == NULL) {
      Error("Allocation mémoire échouée pour rcm_order ou visited");
  }
  int rcm_index = 0;

  // Parcourir tous les nœuds pour traiter chaque composante connexe.
  for (int start = 0; start < nNodes; start++) {
      if (!visited[start]) {
          // Pour cette composante, on choisit le nœud de départ (ici, start).
          int root = start;
          // Utilisation d'une file pour le parcours en largeur (BFS)
          int *queue = malloc(nNodes * sizeof(int));
          if (queue == NULL) { Error("Allocation mémoire échouée pour la file"); }
          int q_start = 0, q_end = 0;
          queue[q_end++] = root;
          visited[root] = 1;
          while (q_start < q_end) {
              int curr = queue[q_start++];
              rcm_order[rcm_index++] = curr;

              // Récupérer les voisins non visités de curr
              int numNeighbors = 0;
              int *temp = malloc(adj[curr].count * sizeof(int));
              if (temp == NULL) { Error("Allocation mémoire échouée pour temp"); }
              for (int i = 0; i < adj[curr].count; i++) {
                  int nb = adj[curr].neighbors[i];
                  if (!visited[nb]) {
                      temp[numNeighbors++] = nb;
                      visited[nb] = 1;
                  }
              }
              // Trier les voisins par degré croissant
              for (int i = 0; i < numNeighbors - 1; i++) {
                  for (int j = i + 1; j < numNeighbors; j++) {
                      if (adj[temp[i]].count > adj[temp[j]].count) {
                          int swap = temp[i];
                          temp[i] = temp[j];
                          temp[j] = swap;
                      }
                  }
              }
              // Ajouter les voisins triés à la file
              for (int i = 0; i < numNeighbors; i++) {
                  queue[q_end++] = temp[i];
              }
              free(temp);
          }
          free(queue);
      }
  }
  free(visited);

  // Inverser l'ordre obtenu pour obtenir l'ordre RCM
  for (int i = 0; i < nNodes / 2; i++) {
      int tmp = rcm_order[i];
      rcm_order[i] = rcm_order[nNodes - i - 1];
      rcm_order[nNodes - i - 1] = tmp;
  }

  /* --- 3. Construction du mapping à partir de l'ordre RCM --- */
  // mapping[old_index] = new_index
  int *mapping = malloc(nNodes * sizeof(int));
  if (mapping == NULL) {
      Error("Allocation mémoire échouée pour mapping");
  }
  for (int i = 0; i < nNodes; i++) {
      mapping[rcm_order[i]] = i;
  }

  /* --- 4. Réorganisation des tableaux de coordonnées --- */
  double *newX = malloc(nNodes * sizeof(double));
  double *newY = malloc(nNodes * sizeof(double));
  if (newX == NULL || newY == NULL) {
      Error("Allocation mémoire échouée pour les nouveaux tableaux de coordonnées");
  }
  for (int i = 0; i < nNodes; i++) {
      newX[i] = theGeometry->theNodes->X[rcm_order[i]];
      newY[i] = theGeometry->theNodes->Y[rcm_order[i]];
  }
  free(theGeometry->theNodes->X);
  free(theGeometry->theNodes->Y);
  theGeometry->theNodes->X = newX;
  theGeometry->theNodes->Y = newY;

  /* --- 5. Mise à jour de la connectivité --- */
  // Pour les arêtes
  femMesh *edges = theGeometry->theEdges;
  if (edges) {
      for (int i = 0; i < edges->nElem; i++) {
          for (int j = 0; j < edges->nLocalNode; j++) {
              int oldIndex = edges->elem[i * edges->nLocalNode + j];
              edges->elem[i * edges->nLocalNode + j] = mapping[oldIndex];
          }
      }
  }
  // Pour les éléments
  if (theGeometry->theElements) {
      femMesh *elems = theGeometry->theElements;
      for (int i = 0; i < elems->nElem; i++) {
          for (int j = 0; j < elems->nLocalNode; j++) {
              int oldIndex = elems->elem[i * elems->nLocalNode + j];
              elems->elem[i * elems->nLocalNode + j] = mapping[oldIndex];
          }
      }
  }

  /* --- 6. Libération des ressources temporaires --- */
  for (int i = 0; i < nNodes; i++) {
      free(adj[i].neighbors);
  }
  free(adj);
  free(rcm_order);
  free(mapping);
}


double *femFullSystemEliminate(femFullSystem *mySystem) {
  double **A, *B, factor;
  int i, j, k, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  /* Gauss elimination */

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    for (i = k + 1; i < size; i++) {
      factor = A[i][k] / A[k][k];
      for (j = k + 1; j < size; j++)
        A[i][j] = A[i][j] - A[k][j] * factor;
      B[i] = B[i] - B[k] * factor;
    }
  }

  /* Back-substitution */

  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    for (j = i + 1; j < size; j++)
      factor += A[i][j] * B[j];
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (mySystem->B);
}

double *femFullSystemEliminateBand(femFullSystem *mySystem) {
  double **A = mySystem->A;
  double *B = mySystem->B;
  int size = mySystem->size;
  int band = (int)ComputeBand(mySystem); // Compute the half-bandwidth

  // Forward elimination (only within the band)
  for (int k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    // Only update rows within the band (i from k+1 to k+band)
    int max_i = (k + band < size ? k + band : size - 1);
    for (int i = k + 1; i <= max_i; i++) {
      double factor = A[i][k] / A[k][k];
      // Update only the entries in the band (j from k+1 to k+band)
      int max_j = (k + band < size ? k + band : size - 1);
      for (int j = k + 1; j <= max_j; j++) {
        A[i][j] -= factor * A[k][j];
      }
      B[i] -= factor * B[k];
      // Optionally set A[i][k] to 0 explicitly
      A[i][k] = 0.0;
    }
  }

  // Back-substitution (again using the band information)
  for (int i = size - 1; i >= 0; i--) {
    double sum = 0.0;
    int max_j = (i + band < size ? i + band : size - 1);
    for (int j = i + 1; j <= max_j; j++) {
      sum += A[i][j] * B[j];
    }
    B[i] = (B[i] - sum) / A[i][i];
  }

  return B;
}


femBandSystem *femBandSystemCreate(femFullSystem *mySystem) {
  femBandSystem *femBandSystemCreate = malloc(sizeof(femBandSystem));
  femBandSystemCreate->size = mySystem->size;
  femBandSystemCreate->band = ComputeBand(mySystem);
}
double ComputeBand(femFullSystem *mySystem) {
  double **A;
  int size;
  int band = 0;

  A = mySystem->A;
  size = mySystem->size;

  // Loop over each element of the matrix A.
  for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
          // Check if the element is nonzero
          if (A[i][j] != 0.0) {
              int current_band = abs(i - j);
              if (current_band > band) {
                  band = current_band;
              }
          }
      }
  }
  return (double)band;
}


void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) {
  double **A, *B;
  int i, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    B[i] -= myValue * A[i][myNode];
    A[i][myNode] = 0;
  }

  for (i = 0; i < size; i++)
    A[myNode][i] = 0;

  A[myNode][myNode] = 1;
  B[myNode] = myValue;
}
