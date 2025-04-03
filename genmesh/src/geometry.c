#include "fem.h"
#include "gmshc.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {
  femGeo *theGeometry = geoGetGeometry();
  double dc = 3.5;
  theGeometry->dhCenter = dc;
  double hc = 0.5;
  theGeometry->hCenter = hc;
  double forcePosx = 50.0;
  theGeometry->forcePositionX = forcePosx;
  double forcePosy = 0.0;
  theGeometry->forcePositionY = forcePosy;
  double forceR = 30.0;
  theGeometry->forceRadius = forceR;
  double dt = 10.0;
  theGeometry->dhTooth = dt;
  double ht = 0.5;
  theGeometry->hTooth = ht;
  double h = 5.0;
  theGeometry->h = h;
  double curvRatio = theGeometry->curvatureRatio;
  double curv = theGeometry->curvature;
  double rin = theGeometry->Rinner;
  double rout = theGeometry->Router;
  double toothLength = theGeometry->toothL;
  
  double distcenter = sqrt((y)*(y) + (x)*(x));
  double htemp = h;
  if(distcenter <= rin+dc)
  {
    double alpha = curvRatio * (2.0*M_PI/3.0);
    double beta = 2.0*asin(rin*curv*sin(alpha/2.0));
    double dist = rin*cos(alpha/2.0) + cos(beta/2.0)/curv;
    double distInner1 = sqrt((y-dist)*(y-dist) + (x)*(x));
    double distInner2 = sqrt((y-dist*cos(2.0*M_PI/3.0))*(y-dist*cos(2.0*M_PI/3.0)) + (x-dist*sin(2.0*M_PI/3.0))*(x-dist*sin(2.0*M_PI/3.0)));
    double distInner3 = sqrt((y-dist*cos(4.0*M_PI/3.0))*(y-dist*cos(4.0*M_PI/3.0)) + (x-dist*sin(4.0*M_PI/3.0))*(x-dist*sin(4.0*M_PI/3.0)));
    double distEdge = distcenter - rin;

    double tolerance = 1.0;
    if(distInner1 <= 1/curv+tolerance){
      double theta = asin(x/distInner1);
      if(theta < alpha/2.0 && theta > -alpha/2.0){
        distEdge = 1/curv - distInner1;
      }
    }
    if(distInner2 <= 1/curv+tolerance){
      double theta = M_PI/3.0 - atan((dist*sin(2.0*M_PI/3.0)-x)/(y-dist*cos(2.0*M_PI/3.0)));
      if(theta < alpha/2.0 && theta > -alpha/2.0){
        distEdge = 1/curv - distInner2;
      }
    }
    if(distInner3 <= 1/curv+tolerance){
      double theta = -M_PI/3.0 - atan((dist*sin(4.0*M_PI/3.0)-x)/(y-dist*cos(4.0*M_PI/3.0)));
      if(theta < alpha/2.0 && theta > -alpha/2.0){
        distEdge = 1/curv - distInner3;
      }
    }
    if(distEdge <= dc)
    {
      htemp = hc + 3*(h-hc)*((distEdge/dc)*(distEdge/dc)) + 2*(hc-h)*((distEdge/dc)*(distEdge/dc)*(distEdge/dc));
      h = fmin(htemp, h);
    }
  }
  
  htemp = h;
  double distForce = sqrt((y-forcePosy)*(y-forcePosy) + (x-forcePosx)*(x-forcePosx));
  if(distForce <= forceR)
  {
    double distEdge = rout+toothLength-distcenter;
    if(distEdge <= dt)
    {
      htemp = ht + 3*(h-hc)*((distEdge/dt)*(distEdge/dt)) + 2*(hc-h)*((distEdge/dt)*(distEdge/dt)*(distEdge/dt));
      h = fmin(htemp, h);
    }

  }

  return h;

}

void geoMeshGenerate(double lc) {
  femGeo *theGeometry = geoGetGeometry();
  // theGeometry->LxPlate = 1.0;
  // theGeometry->LyPlate = 1.0;
  // theGeometry->h = 1.0 * 0.05;
  double ri = 18.85;
  double ro = 25.0;
  double curv = 1.0/(1.0*ri);
  double curvRatio = 0.4583333333333333;
  int nTeeth = 32;
  double toothLength = 3.0;
  double toothWidth = 1.5;
  theGeometry->Rinner = ri;
  theGeometry->Router = ro;
  theGeometry->curvature = curv;
  theGeometry->curvatureRatio = curvRatio;
  theGeometry->elementType = FEM_TRIANGLE;
  theGeometry->toothL = toothLength;
  theGeometry->toothW = toothWidth;

  geoSetSizeCallback(geoSize);

  int ierr;

  // Create disk structure
  double alpha = curvRatio * (2.0*M_PI/3.0);
  if(ri*sin(alpha/2.0) > 1/curv)
  {
    printf("Error: curvRatio is too high\n");
    return;
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
  geoMeshGenerate(0.1);
  return;
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

void geoAssembleDomains(void){

  
  femGeo *theGeometry = geoGetGeometry();
  double rin = theGeometry->Rinner;
  double rout = theGeometry->Router;
  double limit = rin + (rout-rin)/2.0;
  double forceR = theGeometry->forceRadius;
  double forcePosx = theGeometry->forcePositionX;
  double forcePosy = theGeometry->forcePositionY;
  // printf("Limit: %f\n", limit);

  int domainAppartenance[theGeometry->nDomains];
  // printf("Number of domains: %d\n", theGeometry->nDomains);
  int innerNElem = 0;
  int freeNElem = 0;
  int forceNElem = 0;

  for(int i = 0; i < theGeometry->nDomains; i++)
  {
    domainAppartenance[i] = 0;
    femDomain *currentDomain = theGeometry->theDomains[i];
    int edge = currentDomain->elem[0];
    int nodeA = currentDomain->mesh->elem[2*edge];
    int nodeB = currentDomain->mesh->elem[2*edge+1];
    double xA = theGeometry->theNodes->X[nodeA];
    double yA = theGeometry->theNodes->Y[nodeA];
    double xB = theGeometry->theNodes->X[nodeB];
    double yB = theGeometry->theNodes->Y[nodeB];
    double x = xA + (xB-xA)/2.0;
    double y = yA + (yB-yA)/2.0;
    double dist = sqrt((x)*(x) + (y)*(y));
    double distForce = sqrt((x-forcePosx)*(x-forcePosx) + (y-forcePosy)*(y-forcePosy));
    // printf("Domain %i : %f\n", i, dist);
    if(dist < limit){
      domainAppartenance[i] = 1;
      innerNElem += currentDomain->nElem;
    }
    else if(distForce < forceR && xB < xA && fabs(yB-yA) < 0.35){
      domainAppartenance[i] = 2;
      forceNElem += currentDomain->nElem;
    }
    else{
      freeNElem += currentDomain->nElem;
    }
  }
  //newdomains
  femDomain *freeDomain = malloc(sizeof(femDomain));
  // freeDomain->name = "Free";
  freeDomain->mesh = theGeometry->theEdges;
  freeDomain->nElem = freeNElem;
  freeDomain->elem = malloc(sizeof(int) * freeNElem);
  femDomain *innerDomain = malloc(sizeof(femDomain));
  // freeDomain->name = "Inner";
  innerDomain->mesh = theGeometry->theEdges;
  innerDomain->nElem = innerNElem;
  innerDomain->elem = malloc(sizeof(int) * innerNElem);
  femDomain *forceDomain = malloc(sizeof(femDomain));
  forceDomain->mesh = theGeometry->theEdges;
  forceDomain->nElem = forceNElem;
  forceDomain->elem = malloc(sizeof(int) * forceNElem);
  int freeIndex = 0;
  int innerIndex = 0;
  int forceIndex = 0;
  for(int i = 0; i < theGeometry->nDomains; i++)
  {
    // printf("Domain %d: %d\n", i, domainAppartenance[i]);  
    femDomain *currentDomain = theGeometry->theDomains[i];
    if(domainAppartenance[i] == 1){
      for(int j = 0; j < currentDomain->nElem; j++){
        innerDomain->elem[innerIndex] = currentDomain->elem[j];
        innerIndex++;
      }
    }
    else if(domainAppartenance[i] == 2){
      for(int j = 0; j < currentDomain->nElem; j++){
        forceDomain->elem[forceIndex] = currentDomain->elem[j];
        forceIndex++;
      }
    }
    else{
      for(int j = 0; j < currentDomain->nElem; j++){
        freeDomain->elem[freeIndex] = currentDomain->elem[j];
        freeIndex++;
      }
    }
    currentDomain->nElem = 0;
    free(currentDomain->elem);
    currentDomain->mesh = NULL;
    free(currentDomain);
  }
  theGeometry->theDomains = realloc(theGeometry->theDomains, sizeof(femDomain *) * 3);
  theGeometry->theDomains[0] = freeDomain;
  theGeometry->theDomains[1] = innerDomain;
  theGeometry->theDomains[2] = forceDomain;
  theGeometry->nDomains = 3;

  return;
}