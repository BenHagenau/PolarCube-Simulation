#include <stdio.h>
#include <string.h>
#include <math.h>
#include "GPS/vector3D.h"
#include "svd.h"
#include "mex.h"

/* --------------SPACECRAFT PARAMETERS----------------*/

double Gs[][3] = {{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
double Gt[][3] = {{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
double Js[]    = {1.1249e-05, 1.1249e-05, 1.1249e-05};

double SC_I[][3] = {{0.05220503,  -0.00026736,  0.000475987},
    {-0.00026736, 0.052387143,  -0.000342469},
    {0.000475987, -0.000342469, 0.004468478}};

double I_VSCMG[][3] = {{2.2596e-05, 0, 0},
    {0, 2.2596e-05, 0},
    {0, 0, 2.2596e-05}};


double SV_OmegaNom[] = {0.663145596216231, 0.663145596216231, 0.663145596216231};

/* --------------SIMULATION CONDITIONS----------------*/

// int contFilt = 0;
double eps = 1e-14;
double I[3][3];
/*----------------PROTOTYPING-------------------------*/
void ctrl_law(double errMRP[3], double z[3], double omega[3], double errOmega[3], double omegaRef[3],
        double omegaDotRef[3], double L[3], double SV_Omega[3], double B_body[3],
        double K, double P[3][3], double Km, double Ki[3][3], double deadband, double *LrT, double *tau_mag);

/*--------------FUNCTIONS USED BY CTRL_LAW------------*/

// computes us (Lr) from equation 15, ignoring the Gs matrix which comes back in at the final RW torque sum
// where us = Lr + u_s^* + delta u
static void Lr_RW(double Lr[3], double errMRP[3], double z[3], double omega[3], double errOmega[3],
        double omegaRef[3], double omegaDotRef[3], double L[3], double SV_Omega[3], double B_body[3],
        double K, double P[3][3], double Km, double Ki[3][3])
{

    double m_tmp[3][3];
    double m1_tmp[3][3];
    /* double m2_tmp[3][3]; */

    double v_tmp[3];
    double v1_tmp[3];
    double v2_tmp[3];

    // part 1 => -I*(omegaDotRef - tilde(omega)*omegaRef)
    double first[3];

    tilde(omega, m_tmp);
    Mdot(m_tmp, omegaRef, v_tmp);
    sub(omegaDotRef, v_tmp, v_tmp);
    Mdot(I, v_tmp, first); // first must be subtracted

    // part 2 => K*errMRP + P*errOmega + P*Ki*z
    double second[3];

    mult(K, errMRP, v_tmp);
    Mdot(P, errOmega, v1_tmp);
    MdotM(P, Ki, m_tmp);
    Mdot(m_tmp, z, v2_tmp);
    add(v_tmp, v1_tmp, second);
    add(second, v2_tmp, second);

    // part 3 => - (tilde(omegaRef) - tilde(Ki*z))*(I*omega + Gs*hs) + L
    double third[3];

    tilde(omegaRef, m_tmp);
    Mdot(Ki, z, v_tmp);
    tilde(v_tmp, m1_tmp);
    Msub(m_tmp, m1_tmp, m_tmp);

    Mdot(I, omega, v_tmp);

    // Gs is symmetric in this case, still doing it
    transpose(Gs, m1_tmp); // hs = diag(Js)*(Omega + Gs'*omega)
    Mdot(m1_tmp, omega, v1_tmp);
    add(SV_Omega, v1_tmp, v1_tmp);
    diag(Js, m1_tmp);
    Mdot(m1_tmp, v1_tmp, v2_tmp);

    Mdot(Gs, v2_tmp, v1_tmp);

    add(v_tmp, v1_tmp, v_tmp);
    Mdot(m_tmp, v_tmp, v1_tmp);
    sub(L, v1_tmp, third);

    // all together
    add(second, third, Lr);
    sub(Lr, first, Lr);
}

static void pinv(double m[3][3])
{
    double v[3][3];
    double s[3];
    double m_tmp[3][3];
    double m1_tmp[3][3];

    dsvd(m, 3, 3, s, v); // SVD -- u goes into m, s and v are as expected

    memset(m_tmp, 0, sizeof(double [3][3]));

    // create the "inverse" diagonal for S
    size_t i;
    for (i = 0; i < 3; i++)
        m_tmp[i][i] = (fabs(s[i]) < eps) ? 0 : 1/s[i];

    // V*(S"inv")*U'
    MdotM(v, m_tmp, m1_tmp);
    transpose(m, m_tmp);
    MdotM(m1_tmp, m_tmp, m);
}

static void mom_manage(double uStar[3], double delta_u[3], double tau_mag[3], double SV_Omega[3], double B_body[3], double Km)
{
    double m_tmp[3][3];
    double m1_tmp[3][3];
    double m2_tmp[3][3];

    double v_tmp[3];
    double v1_tmp[3];

    // uStar = -Km*diag(Js)*delOmega
    tilde(B_body, m2_tmp); // B_body skew is saved in m2
    diag(Js, m_tmp);
    sub(SV_Omega, SV_OmegaNom, v_tmp);
    Mdot(m_tmp, v_tmp, v1_tmp);
    mult(-1*Km, v1_tmp, uStar);

    // muStar = -pinv(tilde(B_body)*Gt)*Gs*uStar
    MdotM(m2_tmp, Gt, m_tmp);
    pinv(m_tmp); // inverse is saved in m_tmp
    MdotM(m_tmp, Gs, m1_tmp);
    Mdot(m1_tmp, uStar, v_tmp);
    mult(-1, v_tmp, v_tmp); // muStar is now in v_tmp

    // tau_mag = -tilde(B_body)*Gt*muStar
    MdotMT(m2_tmp, Gt, m_tmp);
    Mdot(m_tmp, v_tmp, tau_mag);
    mult(-1, tau_mag, tau_mag);

    // delta_u = inv(Gs)\(tau_mag - Gs*uStar)
    // the inv(Gs)\ ... is suposed to be a minimum norm inverse but basically evaluates to Gs
    Mdot(Gs, uStar, v_tmp);
    sub(tau_mag, v_tmp, v_tmp);
    Mdot(Gs, v_tmp, delta_u);
}


/*----------------CONTROL LAW FUNCTION----------------*/
void ctrl_law(double errMRP[3], double z[3], double omega[3], double errOmega[3],
        double omegaRef[3], double omegaDotRef[3], double L[3],
        double SV_Omega[3], double B_body[3], double K, double P[3][3],
        double Km, double Ki[3][3], double deadband, double *LrT, double *tau_mag)
{    
    Madd(SC_I, I_VSCMG, I);
    double Lr[3];
    double uStar[3];
    double delta_u[3];
    double v_tmp[3];
    static size_t i;
   
    // compute basic RW torques
    Lr_RW(Lr, errMRP, z, omega, errOmega, omegaRef, omegaDotRef, L, SV_Omega, B_body, K, P, Km, Ki);
    // compute the torquer torques and the rest of the wheel torques for momentum management
    mom_manage(uStar, delta_u, tau_mag, SV_Omega, B_body, Km);

    // LrT = Lr + uStar + delta_u -- total wheel torques
    add(Lr, uStar, v_tmp);
    add(v_tmp, delta_u, v_tmp);
    /* Mdot(Gs, v_tmp, LrT); // this one should be correct but thats not how it is in the sim*/
    equal(v_tmp, LrT);

    // check to see if output should be filtered
    // simple running average
    /*if (contFilt) {
        add(LrT_prev[0], LrT_prev[1], v_tmp);
        add(LrT, v_tmp, LrT);
        mult(1/3.0, LrT, LrT);
    }*/

    // update previous LrT values
    //equal(LrT, LrT_prev[i]);
    //i = !i;

    if (norm(LrT) < deadband)
        setZero(LrT);
}

// -------------------BUILDING MEX GATEWAY FUNCTION TO LAUNCH FROM MATLAB---------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    if(nrhs != 14) {
        mexErrMsgIdAndTxt("MyToolbox:ctrl_law:nrhs","14 inputs required.");
    }
    
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:ctrl_law:nlhs", "Two outputs required");
    }
    
    // Create pointers to matlab inputs
    double *errMPR = mxGetPr(prhs[0]);
    double *z = mxGetPr(prhs[1]);
    double *Omega = mxGetPr(prhs[2]);
    double *errOmega = mxGetPr(prhs[3]);
    double *omegaRef = mxGetPr(prhs[4]);
    double *omegaDotRef = mxGetPr(prhs[5]);
    double *L = mxGetPr(prhs[6]);
    double *SV_Omega = mxGetPr(prhs[7]); // wheel rates
    double *B_body = mxGetPr(prhs[8]);
    double *K        = mxGetPr(prhs[9]);
    double *P_mex   = mxGetPr(prhs[10]);
    double *Km       = mxGetPr(prhs[11]);
    double *Ki_mex = mxGetPr(prhs[12]);
    double *deadband = mxGetPr(prhs[13]);
    
    double P[3][3] = {{0}};
    double Ki[3][3] = {{0}};
    
    int i;
    for (i = 0; i < 3; i++){
        P[i][i] = P_mex[i];
        Ki[i][i] = Ki_mex[i];
    }
       
    // Create output pointers
    plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL);
    
    // Get pointers to the output matrices
    double *outMatrix1 = mxGetPr(plhs[0]);
    double *outMatrix2 = mxGetPr(plhs[1]);
    
    // Launch the c function
    ctrl_law(errMPR, z, Omega, errOmega, omegaRef, omegaDotRef, L,
            SV_Omega, B_body, *K, P, *Km, Ki, *deadband, outMatrix1, outMatrix2);
}
