#ifndef _WERKZEUGMASCHINE_H
#define _WERKZEUGMASCHINE_H

#include "mbsim/dynamic_system_solver.h"
#include "mbsimFlexibleBody/flexible_body/flexible_body_1s_21_rcm.h"
#include "mbsim/rigid_body.h"
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
    #include <fmatvec/csparse.h>
#ifdef __cplusplus
}
#endif

class System : public MBSim::DynamicSystemSolver {
  public:
    System(const std::string &projectName);

    /**
     * \brief compress the W matrix into csparse compressed-column form
     * \param time
     */
    cs * compressWToCsparse(int j=0);

    /**
     * \brief compress the W matrix into csparse compressed-column form
     * \param time
     */
    cs * compressLLM_LToCsparse(int j=0);

    virtual void updateG(double t, int i=0);

  private:
    /** flexible ring */
    MBSimFlexibleBody::FlexibleBody1s21RCM *rod;

    /** vector of balls */
    std::vector<MBSim::RigidBody*> balls;
};

#endif

