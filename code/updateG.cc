#include "system.h"
#include "mbsim/joint.h"
#include "mbsim/contact.h"
#include "mbsim/contours/point.h"
#include "mbsim/contours/plane.h"
#include "mbsim/constitutive_laws.h"
#include "mbsim/utils/rotarymatrices.h"
#include "mbsim/environment.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/spineextrusion.h>
#include <openmbvcppinterface/cuboid.h>
#include <openmbvcppinterface/polygonpoint.h>
#endif

#include <mbsim/utils/eps.h>
#include <mbsim/utils/stopwatch.h>

using namespace MBSimFlexibleBody;
using namespace MBSim;
using namespace fmatvec;
using namespace std;

Mat cs2Mat(cs* sparseMat) {
  Mat newMat(sparseMat->m, sparseMat->n, INIT, 0.);

  if (sparseMat->nz < 0) {
    for (int j = 0; j < sparseMat->n; j++) {
      for (int p = sparseMat->p[j]; p < sparseMat->p[j + 1]; p++) {
        int row = sparseMat->i[p];
        int col = j;
        double entry = sparseMat->x ? sparseMat->x[p] : 1;
        newMat(row, col) = entry;
      }
    }
  }
  else {
    for (int i = 0; i < sparseMat->nz; i++) {
      int row = sparseMat->i[i];
      int col = sparseMat->p[i];
      double entry = sparseMat->x ? sparseMat->x[i] : 1;
      newMat(row, col) = entry;
    }
  }

  return newMat;
}

void System::updateG(double t, int j) {
  
/* Original Implementation */
//    cout << "DynamicSystemSolver::updateG(): W[" << j<< "]: "<< W[j] << endl;
//    cout << "DynamicSystemSolver::updateG(): LLM[" << j<< "]: "<< LLM[j] << endl;
//    cout << "DynamicSystemSolver::updateG(): V[" << j<< "]: "<< V[j] << endl;

//  if (j == 0) {
//    static int ITER = 0;
//    std::ostringstream strs;
//    strs << ITER++;
//    std::string str = strs.str();
//    ofstream Ofile(str.c_str());
//    Ofile << LLM[j];
//    Ofile.close();
//  }

//     G << SqrMat(W[j].T() * slvLLFac(LLM[j], V[j])); //TODO@zhan: this is the line!
//     ofstream WOfile("G_original.txt");
//     WOfile << G;
//     WOfile.close();

/* new Implementation with Sparse Matrix technique */
  StopWatch sw;
  sw.start();
  cs *cs_L_LLM, *cs_Wj;

  int n = LLM[j].cols();

  cs_L_LLM = compressLLM_LToCsparse(t, j); // change the Matrix from the FMATVEC format to Csparse format

//  cs_print(cs_L_LLM, 0);

// construct the w in the triplet form here and using cs_triplet (const cs *T) to transform
// or directly transform the w
  cs_Wj = compressWToCsparse(j);

  int * xi = (int *) cs_malloc(2 * n, sizeof(int));

  Mat y(n, W[j].cols(), INIT);
  double * yele = y.operator ()();
  int ylda = y.ldim();

// solve Ly=w
  for (int k = 0; k < W[j].cols(); k++) {
    cs_splsolve (cs_L_LLM, cs_Wj, k, xi, yele + ylda * k , 0); // yele + ylda * k is the pointer to the kth column of y.ele. 
                                                               // Directly stores the data to the pointer of y in order to
                                                               // avoid to create temporary vector x.
//    double * x = (double *) calloc(n, sizeof(double));
//    cs_splsolve(cs_L_LLM, cs_Wj, k, xi, x, 0);
//    for (int iter = 0; iter < n; iter++) {
//      y(iter, k) = x[iter];
//    }
//    free(x);
  }

/* calculate G and copy to Gs 
 * Todo: Do we need sparse matrix multiplication for G = y'*y, when G is a dense matrix. check by matlab. (Done!)
 * Result: the matlab test shows that with full matrix format is faster, this is because that y is stored in compressed-column format and it is hard to access the the column of y'.
 * maybe it is better to store the y into the format compressed-column based on fmatvec and then implement the multiplication, it is different than the normal sparse matrix multiplication
 * if G is symmetric type, only the elements in the lower triangular is stored in column wise.
 * if it is accessed in row wise order in the upper triangular part, is will be re-mapped to the lower part,
 * so ele can be still one by one in the order how they stored in memory.
 **/
  G.resize(y.cols(), NONINIT);

  for (int i = 0; i < y.cols(); i++) {
    Vec temp = y.col(i); // temp can be kept in the fast memmory part
    for (int j = i; j < y.cols(); j++) {
      double val = scalarProduct(temp, y.col(j));
      G(i, j) = val;
      G(j, i) = val;
    }
  }

/* test if the block matrix technique can improve the performance */ 
//    Mat temp(2, n , NONINIT);
//    SqrMat2 val(NONINIT);
//    for (int i = 0; i < y.cols(); i = i + 2) {
//      temp = y(Index(0, n-1), Index(i, i + 1)).T(); // temp can be kept in the fast memmory part
//    for (int j = 0; j < y.cols(); j = j + 2) {
//      G(Index(i,i+1), Index(j,j+1))= temp * y(Index(0, n-1), Index(j, j + 1));
////      G(Index(i,i+1), Index(j,j+1)) = val;
////      G(Index(j,j+1), Index(i,i+1)) = val;
//    }
//  }

/* write the result in to text file for validation and compare the performance of two implementations*/
  double t_cs = sw.stop(true);
  cout << nrm1(W[j] - cs2Mat(cs_Wj)) << endl;
  cout << W[j] - cs2Mat(cs_Wj) << endl;
  cout << nrm1(LLM[j] - cs2Mat(cs_L_LLM)) << endl;

  ofstream Ofile("W_reference.txt");
  Ofile << W[j];
  Ofile.close();
  Ofile.open("W_cs.txt");
  Ofile << cs2Mat(cs_Wj);
  Ofile.close();
  Ofile.open("W_diff.txt");
  Ofile << W[j] - cs2Mat(cs_Wj);
  Ofile.close();

  Ofile.open("LLM_reference.txt");
  Ofile << LLM[j];
  Ofile.close();
  Ofile.open("LLM_cs.txt");
  Ofile << cs2Mat(cs_L_LLM);
  Ofile.close();
  Ofile.open("LLM_diff.txt");
  Ofile << LLM[j] - cs2Mat(cs_L_LLM);
  Ofile.close();

  ofstream ycsOfile("y_cs.txt");
  ycsOfile << y;
  ycsOfile.close();
  sw.start();
  SqrMat Greference(W[j].T() * slvLLFac(LLM[j], V[j]));
  double t_ref = sw.stop(true);
  cout << "t_cs = " << t_cs << " s   | t_ref =" << t_ref << "s   | dt = " << (t_ref - t_cs) / t_ref * 100 << " %" << endl;

  ofstream GOfile("G_cs.txt");
  GOfile << G;
  GOfile.close();
  ofstream GOreffile("G_reference.txt");
  GOreffile << Greference;
  GOreffile.close();
  SqrMat Gcmp(G - Greference);
  ofstream GOcmpfile("G_diff.txt");
  GOcmpfile << Gcmp;
  GOcmpfile.close();
  cout << nrm1(Gcmp) << endl;

  //todo: free the allocated memory by the cSparse
  cs_spfree(cs_Wj);
  cs_spfree(cs_L_LLM);
  free(xi);

  if (checkGSize)
    ;  // Gs.resize();
  else if (Gs.cols() != G.size()) {
    static double facSizeGs = 1;
    if (G.size() > limitGSize && fabs(facSizeGs - 1) < epsroot())
      facSizeGs = double(countElements(G)) / double(G.size() * G.size()) * 1.5;
    Gs.resize(G.size(), int(G.size() * G.size() * facSizeGs));
  }
  Gs << G;
}

cs * System::compressWToCsparse(int j) {

  int m, n, nz, *Cp, *Ci, counter;;
  double *Cx;;
  cs *C;
  double EPSILON = 1e-17; /*todo: a better epsilion? */

  m = W[j].rows();
  n = W[j].cols();
  nz = 7 * m; // todo: find a method to determine the value of nz

  C = cs_spalloc(m, n, nz, 1, 0); /* allocate memory for the sparse matrix C in Triplet formate, the last input value is 0!*/
  if (!C)
    return (cs_done(C, 0, 0, 0)); /* out of memory */

  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  counter = 0;
  for (vector<Link*>::iterator i = linkSetValued.begin(); i != linkSetValued.end(); ++i) {
    LinkMechanics * i_temp = dynamic_cast<LinkMechanics*>(*i);
    for (int col = (**i).getlaInd(); col < (**i).getlaInd() + (**i).getlaSize(); col++) {
      Cp[col] = counter;
      if ((**i).getType() == "Joint") {  // for joint
        const size_t Noframes = (*i_temp).getFrame().size();
        for (size_t partner = 0; partner < Noframes; partner++) {
          int lowerRow = (*i_temp).getFrame()[partner]->gethInd(j);
          int upperRow = (*i_temp).getFrame()[partner]->gethInd(j) + (*i_temp).getFrame()[partner]->gethSize(j);

          for (int row = lowerRow; row < upperRow; row++) {
            double entry = W[j](row, col);
            if (fabs(entry) > EPSILON){
              Ci[counter] = row;
              Cx[counter] = entry;
              counter++;
            }
          }
        }
      }
      else if ((**i).getType() == "Contact") {  // for contour
        const size_t NoContacts = (*i_temp).getContour().size();
        for (size_t partner = 0; partner < NoContacts; partner++) {
          int lowerRow = (*i_temp).getContour()[partner]->gethInd(j);
          int upperRow = (*i_temp).getContour()[partner]->gethInd(j) + (*i_temp).getContour()[partner]->gethSize(j);

          for (int row = lowerRow; row < upperRow; row++) {
            double entry = W[j](row, col);
            if (fabs(entry) > EPSILON){
              Ci[counter] = row;
              Cx[counter] = entry;
              counter++;
            }
          }
        }
      }

      else {
        throw MBSimError("Not implemented!");
      }
    }
  }

  Cp[n] = counter;

  // free the allocated space
  return C;

}

cs * System::compressLLM_LToCsparse(double t, int j) {

  int n, nz, *Cp, *Ci, counter;;
  double *Cx;
  cs *LLM_L_cs;
  double EPSILON = 1e-17; /*todo: a better epsilion? */

  n = LLM[j].cols();

  nz = 7 * n;  // cs_entry will dynamic extend (double) the memory size by realloc(), which will consider whether there are still enough space after the original position.

  // creat matrix
  LLM_L_cs = cs_spalloc(n, n, nz, 1, 0); /* allocate memory for the sparse matrix C in Triplet formate, the last input value is 0!*/
  if (!LLM_L_cs)
    return (cs_done(LLM_L_cs, 0, 0, 0)); /* out of memory */

  Cp = LLM_L_cs->p;
  Ci = LLM_L_cs->i;
  Cx = LLM_L_cs->x;

  //create the three arrays
  counter = 0;
  // Run loop over all sub systems
  for (int i = 0; i < (int) dynamicsystem.size(); i++) {
    int uInd = dynamicsystem[i]->getuInd(j);
    int uSize = dynamicsystem[i]->getuSize(j);

    for (int col = uInd; col < uInd + uSize; col++) {
      Cp[col] = counter;
      for (int row = col; row < uInd + uSize; row++) {
        double entry = LLM[j](row, col);
        if (fabs(entry) > EPSILON){
          Ci[counter] = row;
          Cx[counter] = entry;
          counter++;
        }
      }
    }

  }

  // Run loop over all objects
  for (int i = 0; i < (int) object.size(); i++) {
    int uInd = object[i]->getuInd(j);
    int uSize = object[i]->getuSize(j);

    for (int col = uInd; col < uInd + uSize; col++) {
      Cp[col] = counter;
      for (int row = col; row < uInd + uSize; row++) {
        double entry = LLM[j](row, col);
        if (fabs(entry) > EPSILON){
          Ci[counter] = row;
          Cx[counter] = entry;
          counter++;
        }
      }
    }
  }

  Cp[n] = counter;

  return LLM_L_cs;

}

System::System(const string &projectName) :
    DynamicSystemSolver(projectName) {

  // acceleration of gravity
  Vec grav(3, INIT, 0.);
  grav(1) = -9.81;
  MBSimEnvironment::getInstance()->setAccelerationOfGravity(grav);

  // input flexible ring
  double l0 = 1.; // length ring
  double E = 2.5e9; // E-Modul alu
  double rho = 2.5e3; // density alu
  int elements = 32; // number of finite elements
  double b0 = 0.02; // width
  double A = b0 * b0; // cross-section area
  double I = 1. / 12. * b0 * b0 * b0 * b0; // moment inertia

  // input infty-norm balls (cuboids)
  int nBalls = 80; // number of balls
  double mass = 0.025; // mass of ball

  // flexible ring
  rod = new FlexibleBody1s21RCM("Rod", false);
  rod->setLength(l0);
  rod->setEModul(E);
  rod->setCrossSectionalArea(A);
  rod->setMomentInertia(I);
  rod->setDensity(rho);
  rod->setFrameOfReference(this->getFrame("I"));
  rod->setNumberElements(elements);
  rod->initRelaxed(M_PI / 2.);
  this->addObject(rod);

#ifdef HAVE_OPENMBVCPPINTERFACE
  OpenMBV::SpineExtrusion *cuboid=new OpenMBV::SpineExtrusion;
  cuboid->setNumberOfSpinePoints(elements*4); // resolution of visualisation
  cuboid->setStaticColor(0.5);// color in (minimalColorValue, maximalColorValue)
  cuboid->setScaleFactor(1.);// orthotropic scaling of cross section
  vector<OpenMBV::PolygonPoint*> *rectangle = new vector<OpenMBV::PolygonPoint*>;// clockwise ordering, no doubling for closure
  OpenMBV::PolygonPoint *corner1 = new OpenMBV::PolygonPoint(b0*0.5,b0*0.5,1);
  rectangle->push_back(corner1);
  OpenMBV::PolygonPoint *corner2 = new OpenMBV::PolygonPoint(b0*0.5,-b0*0.5,1);
  rectangle->push_back(corner2);
  OpenMBV::PolygonPoint *corner3 = new OpenMBV::PolygonPoint(-b0*0.5,-b0*0.5,1);
  rectangle->push_back(corner3);
  OpenMBV::PolygonPoint *corner4 = new OpenMBV::PolygonPoint(-b0*0.5,b0*0.5,1);
  rectangle->push_back(corner4);

  cuboid->setContour(rectangle);
  rod->setOpenMBVSpineExtrusion(cuboid);
#endif

  // balls
  assert(nBalls > 1);
  double d = 7. * l0 / (8. * nBalls); // thickness
  double b = b0 * 1.5; // height / width

  for (int i = 0; i < nBalls; i++) {
    stringstream name;
    name << "Element_" << i;
    RigidBody *ball = new RigidBody(name.str());
    balls.push_back(ball);
    balls[i]->setFrameOfReference(this->getFrame("I"));
    balls[i]->setFrameForKinematics(balls[i]->getFrame("C"));
    balls[i]->setTranslation(new LinearTranslation("[1,0;0,1;0,0]"));
    balls[i]->setRotation(new RotationAboutFixedAxis(Vec("[0;0;1]")));
    balls[i]->setMass(mass);
    SymMat Theta(3, INIT, 0.);
    Theta(0, 0) = 1. / 6. * mass * b * b;
    Theta(1, 1) = 1. / 12. * mass * (d * d + b * b);
    Theta(2, 2) = 1. / 12. * mass * (d * d + b * b);
    balls[i]->setInertiaTensor(Theta);
    this->addObject(balls[i]);

    Point *pt = new Point("COG");
    balls[i]->addContour(pt, Vec(3, INIT, 0.), SqrMat(3, EYE), balls[i]->getFrame("C"));

    Point *tP = new Point("topPoint");
    balls[i]->addContour(tP, d * Vec("[0.5;0;0]") + b * Vec("[0;0.5;0]"), SqrMat(3, EYE), balls[i]->getFrame("C"));

    Point *bP = new Point("bottomPoint");
    balls[i]->addContour(bP, d * Vec("[0.5;0;0]") - b * Vec("[0;0.5;0]"), SqrMat(3, EYE), balls[i]->getFrame("C"));

    Plane *plane = new Plane("Plane");
    SqrMat trafoPlane(3, INIT, 0.);
    trafoPlane(0, 0) = -1.;
    trafoPlane(1, 1) = 1.;
    trafoPlane(2, 2) = -1.;
    balls[i]->addContour(plane, -d * Vec("[0.5;0;0]"), trafoPlane, balls[i]->getFrame("C"));

#ifdef HAVE_OPENMBVCPPINTERFACE
    OpenMBV::Cuboid *cube=new OpenMBV::Cuboid;
    cube->setLength(d,b,b);
    cube->setStaticColor(1.);
    balls[i]->setOpenMBVRigidBody(cube);
#endif
  }

  //Set balls to correct position
  FlexibleBody1s21RCM * rodInfo = new FlexibleBody1s21RCM("InfoRod", false);

  rodInfo->setq0(rod->getq());
  rodInfo->setu0(rod->getu());
  rodInfo->setNumberElements(rod->getNumberElements());
  rodInfo->setLength(rod->getLength());
  rodInfo->setFrameOfReference(rod->getFrameOfReference());

  rodInfo->initInfo();
  rodInfo->updateStateDependentVariables(0.);

  for (unsigned int i = 0; i < balls.size(); i++) {
    Vec q0(3, INIT, 0.);
    double xL = fmod(i * rodInfo->getLength() / balls.size() + rodInfo->getLength() * 0.25, rodInfo->getLength());
    ContourPointData cp;
    cp.getContourParameterType() = CONTINUUM;
    cp.getLagrangeParameterPosition() = Vec(1, INIT, xL);

    rodInfo->updateKinematicsForFrame(cp, position_cosy);
    q0(0) = cp.getFrameOfReference().getPosition()(0);
    q0(1) = cp.getFrameOfReference().getPosition()(1);
    q0(2) = -AIK2Cardan(cp.getFrameOfReference().getOrientation())(2) + M_PI * 0.5;
    balls[i]->setInitialGeneralizedPosition(q0);
  }

  delete rodInfo;

  // inertial ball constraint
  this->addFrame("BearingFrame", l0 / (2 * M_PI) * Vec("[0;1;0]"), SqrMat(3, EYE), this->getFrame("I"));
  Joint *joint = new Joint("BearingJoint");
  joint->setForceDirection(Mat("[1,0;0,1;0,0]"));
  joint->setForceLaw(new BilateralConstraint);
  joint->setImpactForceLaw(new BilateralImpact);
  joint->connect(this->getFrame("BearingFrame"), balls[0]->getFrame("C"));
  this->addLink(joint);

  // constraints balls on flexible band
  for (int i = 0; i < nBalls; i++) {
    Contact *contact = new Contact("Band_" + balls[i]->getName());
    contact->setContactForceLaw(new BilateralConstraint);
    contact->setContactImpactLaw(new BilateralImpact);
    contact->connect(balls[i]->getContour("COG"), rod->getContour("Contour1sFlexible"));
    contact->enableOpenMBVContactPoints(0.01);
    this->addLink(contact);
  }

  // inner-ball contacts
  for (int i = 0; i < nBalls; i++) {
    stringstream namet, nameb;
    namet << "ContactTop_" << i;
    nameb << "ContactBot_" << i;
    Contact *ctrt = new Contact(namet.str());
    Contact *ctrb = new Contact(nameb.str());
    ctrt->setContactForceLaw(new UnilateralConstraint);
    ctrt->setContactImpactLaw(new UnilateralNewtonImpact(0.));
    ctrb->setContactForceLaw(new UnilateralConstraint);
    ctrb->setContactImpactLaw(new UnilateralNewtonImpact(0.));
    if (i == nBalls - 1) {
      ctrt->connect(balls[0]->getContour("topPoint"), balls[i]->getContour("Plane"));
      ctrb->connect(balls[0]->getContour("bottomPoint"), balls[i]->getContour("Plane"));
    }
    else {
      ctrt->connect(balls[i + 1]->getContour("topPoint"), balls[i]->getContour("Plane"));
      ctrb->connect(balls[i + 1]->getContour("bottomPoint"), balls[i]->getContour("Plane"));
    }
    this->addLink(ctrt);
    this->addLink(ctrb);
  }
}

/* other methods to compress W matrix into Csparse format
//cs * System::compressWToCsparse(int j) {
//
//  int m, n, nz; // row, column, nz;//, I1;
//  cs *C;
//  double EPSILON = 1e-17; /*todo: a better epsilion? */
//
//  m = W[j].rows();
//  n = W[j].cols();
//  nz = 5 * m; // todo: find a method to determine the value of nz
//
//  C = cs_spalloc(m, n, nz, 1, 1); /* allocate memory for the sparse matrix C in Triplet formate, the last input value is 1!*/
//  if (!C)
//    return (cs_done(C, 0, 0, 0)); /* out of memory */
//
//  // todo: need to refine when considering the sub DynamicSystem
//  // convert the upper part of W matrix (which corresponding to the DOFs of flexible belt) into cs triplet format
////  I1 = object[0]->getuInd(j) + object[0]->getuSize(j) - 1;
////  for (column = 0; column < n; column++) {
////    for (row = 0; row <= I1; row++) {  // i = 0: scan the whole rectangle area.
////      double entry = W[j](row, column);
////      if (fabs(entry) > EPSILON)
////        cs_entry(C, row, column, entry);  //REMARK: why won't you use this function in the compress_LLM-function?
////
////    }
////  }
//
//  // convert the bottom part of W matrix into cs triplet format
////  int LinkPartner1Bottom = 0;
////  int LinkPartner2Upper = 0;
//  // only consider the contour in the bottom part of W matrix
////  int LinkPartner1B = 0;
////  int LinkPartner2U = 0;
//
//  for (vector<Link*>::iterator i = linkSetValued.begin(); i != linkSetValued.end(); ++i) {
//    LinkMechanics * i_temp = dynamic_cast<LinkMechanics*>(*i);
//    for (int col = (**i).getlaInd(); col < (**i).getlaInd() + (**i).getlaSize(); col++) {
//      if ((**i).getType() == "Joint") {  // for joint
//        const size_t Noframes = (*i_temp).getFrame().size();
//        for (size_t partner = 0; partner < Noframes; partner++) {
//          int lowerRow = (*i_temp).getFrame()[partner]->gethInd(j);
//          int upperRow = (*i_temp).getFrame()[partner]->gethInd(j) + (*i_temp).getFrame()[partner]->gethSize(j);
//
//          for (int row = lowerRow; row < upperRow; row++) {
//            double entry = W[j](row, col);
//            if (fabs(entry) > EPSILON)
//              cs_entry(C, row, col, entry);
//          }
//        }
////        LinkPartner1Bottom =;
////        LinkPartner2Upper =;
//        // only consider the frame in the bottom part of W matrix
////        LinkPartner1B = max((*i_temp).getFrame()[0]->gethInd(j), (*i_temp).getFrame()[1]->gethInd(j));
////        LinkPartner2U = max((*i_temp).getFrame()[0]->gethInd(j) + (*i_temp).getFrame()[0]->getJacobianOfTranslation(j).cols(), (*i_temp).getFrame()[1]->gethInd(j) + (*i_temp).getFrame()[1]->getJacobianOfTranslation(j).cols()) - 1;
//      }
//      else if ((**i).getType() == "Contact") {  // for contour
////        LinkPartner1Bottom = (**i).getlaInd();
////        LinkPartner2Upper = (**i).getlaInd() + (**i).getlaSize() - 1;
//      // only consider the contour in the bottom part of W matrix
//
//        const size_t NoContacts = (*i_temp).getContour().size();
//        for (size_t partner = 0; partner < NoContacts; partner++) {
//          int lowerRow = (*i_temp).getContour()[partner]->gethInd(j);
//          int upperRow = (*i_temp).getContour()[partner]->gethInd(j) + (*i_temp).getContour()[partner]->gethSize(j);
//
//          for (int row = lowerRow; row < upperRow; row++) {
//            double entry = W[j](row, col);
//            if (fabs(entry) > EPSILON)
//              cs_entry(C, row, col, entry);
//          }
//        }
////        LinkPartner1B = max((*i_temp).getContour()[0]->gethInd(j), (*i_temp).getContour()[1]->gethInd(j));
////        LinkPartner2U = max((*i_temp).getContour()[0]->gethInd(j) + (*i_temp).getContour()[0]->gethSize(j), (*i_temp).getContour()[1]->gethInd(j) + (*i_temp).getContour()[1]->gethSize(j)) - 1;
//      }
//
//      else {
//        throw MBSimError("Not implemented!");
//      }
//
////      for (column = LinkPartner1Bottom; column <= LinkPartner2Upper; column++) {
////        for (row = LinkPartner1B; row <= LinkPartner2U; row++) {
////          double entry = W[j](row, column);
////          if (fabs(entry) > EPSILON)
////            cs_entry(C, row, col, entry);
////        }
////      }
//    }
//  }
//
//  // compress triplet format into Compress column format
////  cs_print(C, 0);
//  cs * cs_Wj = cs_triplet(C);
////  cs_print(cs_Wj, 0);
//
//  // free the allocated space
//  cs_spfree(C);
//
//  return cs_Wj;
//
//}
//
//
//
//
//cs * System::compressLLM_LToCsparse(double t, int j) {
//
//  int n, nz;
//  cs *LLM_L_csTriplet;
//  double EPSILON = 1e-17; /*todo: a better epsilion? */
//
//  n = LLM[j].cols();
//
//  // check for maximal non-zero-size: right now there are three different ways, when compare the time, the percentage is not reliable for comparing because the ref_time varies a little.
//
//  // this method allocates about 32 times larger memory than needed, and has to calculate the nz at every time
////  for (size_t i = 0; i < dynamicsystem.size(); i++) {
////    nz += dynamicsystem[i]->getuSize(j) * dynamicsystem[i]->getuSize(j);
////  }
////
////  for (size_t i = 0; i < object.size(); i++) {
////    nz += object[i]->getuSize(j) * object[i]->getuSize(j);
////  }
//
//  // this one gives appropriate memory size as needed. nz is only calculate once, but has branches.  this one can be still simplified, is j always equto 0
//// //  if(nz0 == 0 || nz1 == 0){
////  static int nz0, nz1;
////  if(t < 2e-6) { // 2e-6 should be timestep size // todo: find a better way to do this
////    if(j == 0)
////      nz0 = countElementsLT(LLM[j]) + 50;
////    else if (j == 1)
////      nz1 =  countElementsLT(LLM[j]) + 50;
////  }
////  nz = (j==0) ? nz0 : nz1;
//////  cout << "j = " << j << endl; // j is always 0 ??
//////  cout << "nz = " << nz << endl;
//
//  // this one is the simplest one, can avoid extend the memory in most case and does not waste too much memory
//  nz = 7 * n;  // cs_entry will dynamic extend (double) the memory size by realloc(), which will consider whether there are still enough space after the original position.
//
//  // creat matrix
//  LLM_L_csTriplet = cs_spalloc(n, n, nz, 1, 1); /* allocate memory for the sparse matrix C in Triplet formate, the last input value is 1!*/
//  if (!LLM_L_csTriplet)
//    return (cs_done(LLM_L_csTriplet, 0, 0, 0)); /* out of memory */
//
//  // Run loop over all sub systems
//  for (int i = 0; i < (int) dynamicsystem.size(); i++) {
//    int uInd = dynamicsystem[i]->getuInd(j);
//    int uSize = dynamicsystem[i]->getuSize(j);
//
//    for (int row = uInd; row < uInd + uSize; row++) {
//      for (int col = row; col >= uInd; col--) {
////      for (int col = 0; col < uInd + uSize; col++) {
//        double entry = LLM[j](row, col);
//        if (fabs(entry) > EPSILON)
//          cs_entry(LLM_L_csTriplet, row, col, entry);
//      }
//    }
//  }
//
//  // Run loop over all objects
//  for (int i = 0; i < (int) object.size(); i++) {
//    int uInd = object[i]->getuInd(j);
//    int uSize = object[i]->getuSize(j);
//
//    for (int row = uInd; row < uInd + uSize; row++) {
//      for (int col = row; col >= uInd; col--) {
////      for (int col = 0; col < uInd + uSize; col++) {
//        double entry = LLM[j](row, col);
//        if (fabs(entry) > EPSILON)
//          cs_entry(LLM_L_csTriplet, row, col, entry);
//      }
//    }
//  }
//
//  // compress triplet format into Compress column format
////  cs_print(LLM_L_csTriplet, 0);
////  cout << "LLM_L_csTriplet.nzmax = "<< LLM_L_csTriplet->nzmax << endl;
//  cs * LLM_L_cs = cs_triplet(LLM_L_csTriplet);
////  cs_print(LLM_L_cs, 0);
////  cout << "LLM_L_cs.nzmax = "<< LLM_L_cs->nzmax << endl;
//
//
//  // free the allocated space
//  cs_spfree(LLM_L_csTriplet);
//
//  return LLM_L_cs;
//
//}
