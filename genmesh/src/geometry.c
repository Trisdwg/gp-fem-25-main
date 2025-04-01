#include "fem.h"
#include <gmshc.h>

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {

  return 1.0;
}

void geoMeshGenerate(double lc) {
  femGeo *theGeometry = geoGetGeometry();
  theGeometry->LxPlate = 1.0;
  theGeometry->LyPlate = 1.0;
  theGeometry->h = 1.0 * 0.05;
  double ri = 18.85;
  double ro = 25.0;
  double curv = 1.0/(1.0*ri);
  double curvRatio = 0.4583333333333333;
  int nTeeth = 32;
  theGeometry->Rinner = ri;
  theGeometry->Router = ro;
  theGeometry->curvature = curv;
  theGeometry->curvatureRatio = curvRatio;
  theGeometry->elementType = FEM_TRIANGLE;

  geoSetSizeCallback(geoSize);

  int ierr;

  // Create disk structure
  double alpha = curvRatio * (2.0*M_PI/3.0);
  if(ri*sin(alpha/2.0) > 1/curv)
  {
    printf("Error: curvRatio is too high\n");
  }
  double beta = 2.0*asin(ri*curv*sin(alpha/2.0));
  double dist = ri*cos(alpha/2.0) + cos(beta/2.0)/curv;
  int idOutDisk = gmshModelOccAddDisk(0,0,0,ro,ro,-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int idInDisk = gmshModelOccAddDisk(0,0,0,ri,ri,-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int idCurvDisk  = gmshModelOccAddDisk(0,dist,0,1/curv,1/curv,-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int outDisk[] = {2, idOutDisk};
  int inDisk[] = {2, idInDisk};
  int curvDisk[] = {2, idCurvDisk};

  gmshModelOccCut(inDisk, 2,curvDisk, 2, NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr);
  ErrorGmsh(ierr);
  for(int k = 0; k < 2; k++)
  {
    gmshModelOccRotate(curvDisk, 2, 0, 0, 0, 0, 0, 1, 2*M_PI/3, &ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(inDisk, 2,curvDisk, 2, NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr);
    ErrorGmsh(ierr);
  }
  gmshModelOccRemove(curvDisk, 2, 0, &ierr);
  ErrorGmsh(ierr);
  int *holder;
  size_t holderdim;
  gmshModelOccCut(outDisk, 2,inDisk, 2, &holder,&holderdim,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);

  // Create tooth structure
  double toothLength = 3.0;
  double toothWidth = 1.5;
  double toothAngle = 2*M_PI*5.0/360;
  int idTooth = gmshModelOccAddDisk(ro,0,0,toothLength,toothWidth,-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int tooth[] = {2, idTooth};
  int idTorus = gmshModelOccAddDisk(0,0,0,2*ro,2*ro,-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int torus[] = {2, idTorus};
  int idInTorus = gmshModelOccAddDisk(0,0,0,ro+(toothLength/1.25),ro+(toothLength/1.25),-1,NULL,NULL,NULL,NULL, &ierr);
  ErrorGmsh(ierr);
  int inTorus[] = {2, idInTorus};
  gmshModelOccCut(torus, 2,inTorus, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);
  int *endTooth;
  size_t endToothdim;
  gmshModelOccCut(tooth, 2,torus, 2, &endTooth,&endToothdim,NULL,NULL,NULL,-1,1,1,&ierr);
  ErrorGmsh(ierr);

  for(int t = 0; t < nTeeth; t++)
  {
    // printf("t = %d\n", t);
    gmshModelOccRotate(endTooth, endToothdim, 0, 0, 0, 0, 0, 1, (2*M_PI/nTeeth), &ierr);
    ErrorGmsh(ierr);
    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);
    int *temp;
    size_t tempdim;
    gmshModelOccFuse(holder, holderdim,endTooth, endToothdim, &temp,&tempdim,NULL,NULL,NULL,-1,1,0,&ierr);
    ErrorGmsh(ierr);
    holder = temp;
    holderdim = tempdim;
  }
  gmshModelOccRemove(endTooth, endToothdim, 0, &ierr);
  ErrorGmsh(ierr);
  // gmshModelOccRemove(holder, 2, 0, &ierr);
  // ErrorGmsh(ierr);

  gmshModelOccSynchronize(&ierr);

  // Use a frontal delaunay algorithm
  gmshOptionSetNumber("Mesh.Algorithm", 6, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);

  return;
}

void geoMeshGenerateGeo(void) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_TRIANGLE;

  geoSetSizeCallback(geoSize);

  /*
  4 ------------------ 3
  |                    |
  |                    |
  5 ------- 6          |
             \         |
              )        |
             /         |
  8 ------- 7          |
  |                    |
  |                    |
  1 ------------------ 2
  */

  int ierr;
  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;
  double r = w / 4;
  double lc = theGeometry->h;

  int p1 = gmshModelGeoAddPoint(-w / 2, -h / 2, 0., lc, 1, &ierr);
  int p2 = gmshModelGeoAddPoint(w / 2, -h / 2, 0., lc, 2, &ierr);
  int p3 = gmshModelGeoAddPoint(w / 2, h / 2, 0., lc, 3, &ierr);
  int p4 = gmshModelGeoAddPoint(-w / 2, h / 2, 0., lc, 4, &ierr);
  int p5 = gmshModelGeoAddPoint(-w / 2, r, 0., lc, 5, &ierr);
  int p6 = gmshModelGeoAddPoint(0., r, 0., lc, 6, &ierr);
  int p7 = gmshModelGeoAddPoint(0., -r, 0., lc, 7, &ierr);
  int p8 = gmshModelGeoAddPoint(-w / 2, -r, 0., lc, 8, &ierr);
  int p9 = gmshModelGeoAddPoint(0., 0., 0., lc, 9, &ierr); // center of circle

  int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
  int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
  int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
  int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
  int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
  int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
  int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
  int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

  int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed
  int c1[] = {1};
  c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);
  int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
  gmshModelGeoSynchronize(&ierr);

  gmshOptionSetNumber("Mesh.Algorithm", 3, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);
}

void geoMeshGenerateGeoFile(const char *filename) {
  femGeo *theGeometry = geoGetGeometry();
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);

  gmshOptionSetNumber("Mesh.Algorithm", 3, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);
  return;
}

void geoMeshGenerateMshFile(const char *filename) {
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  return;
}
