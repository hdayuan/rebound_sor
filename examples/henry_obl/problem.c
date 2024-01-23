#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const double MEARTH = 3.003e-6;
const double REARTH = 4.263e-5;
const double LHAT[3] = {0.0, 0.0, 1.0};

void print_array_sq(double a[3][3]) {
    printf("[[%.3f, %.3f, %.3f],\n[%.3f, %.3f, %.3f],\n[%.3f, %.3f, %.3f]]\n", a[0][0],a[0][1],a[0][2],a[1][0],a[1][1],a[1][2],a[2][0],a[2][1],a[2][2]);
}

void print_array_l(double a[12]) {
    printf("[%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f,  %.3f, %.3f, %.3f]\n", a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11]);
}

double radians(double degrees) {
    return degrees * M_PI / 180.0;
}

double ts_dot(double x[3], double y[3]) {
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

double mag(double v[3]) {
    return sqrt(ts_dot(v, v));
}

void ts_cross(double x[3], double y[3], double z[3]) {
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = -x[0] * y[2] + x[2] * y[0];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

double norm(double v[3]) {
    return mag(v);
}

void assign(double x[], double y[], int length, int xoffset, int yoffset) {
    /*
    the ODE interface doesn't let us use slices for assignment, so let's
    abstract this out
    */
    int i;
    for (i=0; i<length; i++) {
        x[i + xoffset] = y[i + yoffset];
    }
}

// This function stores transpose
// of A[][] in B[][]
void transpose(double A[3][3], double B[3][3])
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            // Assigns the transpose of element A[j][i] to
            // B[i][j]
            B[i][j] = A[j][i];
}

void matmul(double A[3][3], double B[3], double C[3]) {
    int i;
    for (i = 0; i < 3; i++) {
        C[i] = ts_dot(A[i], B);
    }
}

void euler_to_body(double q1, double q2, double q3, double rot_arr[3][3]) {
    /* we use the ZXZ rotation convention */
    double c1 = cos(q1);
    double s1 = sin(q1);
    double c2 = cos(q2);
    double s2 = sin(q2);
    double c3 = cos(q3);
    double s3 = sin(q3);
    rot_arr[0][0] = c1 * c3 - c2 * s1 * s3;
    rot_arr[0][1] = -c1 * s3 - c2 * c3 * s1;
    rot_arr[0][2] = s1 * s2;
    rot_arr[1][0] = c3 * s1 + c1 * c2 * s3;
    rot_arr[1][1] = c1 * c2 * c3 - s1 * s3;
    rot_arr[1][2] = -c1 * s2;
    rot_arr[2][0] = s2 * s3;
    rot_arr[2][1] = c3 * s2;
    rot_arr[2][2] = c2;
}

/**/

void dydt_rb(struct reb_ode* const ode, double* const yDot, const double* const y, const double t) {

    double k=0.331;
    double J2=1e-3;
    /* double triax=1e-5; */
    double triax=0;
    double Q_tide=1e2;
    double k2=1.5;
    double ms=1;
    double rp=2 * REARTH;
    double mp=4 * MEARTH;
    double tidal_dt;
    int i;
    double s_norm, si, sj, sk;
    double shat_ijk[3];
    double inv_rot[3][3];
    double rotmat[3][3];
    double ihat[3];
    double jhat[3];
    double khat[3];
    double svec_ijk[3];
    double svec_xyz[3];
    double n;
    double a = 0.4;
    double I = k * mp * rp*rp;
    double J = I * (1 + triax);
    double K = I * (1 + J2);
    double year;
    double rhovec[3];
    double l_cross_r[3];
    double s_cross_r[3];
    double prefac_rb;
    double prefac_tide;
    double torque_rb[3];
    double torque_tide_xyz[3];
    double torque_tide[3];
    double torque[3];
    double rho_cross_r[3];

    assign(ihat, y, 3, 0, 0);
    assign(jhat, y, 3, 0, 3);
    assign(khat, y, 3, 0, 6);
    assign(svec_ijk, y, 3, 0, 9);

    s_norm = norm(svec_ijk);
    si = svec_ijk[0]; sj = svec_ijk[1]; sk = svec_ijk[2];
    shat_ijk[0] = si / s_norm; shat_ijk[1] = sj / s_norm; shat_ijk[2] = sk / s_norm;

    for (i=0;i<3;i++) {
        inv_rot[0][i] = ihat[i];
        inv_rot[1][i] = jhat[i];
        inv_rot[2][i] = khat[i];
    }
    transpose(inv_rot, rotmat);
    matmul(rotmat, svec_ijk, svec_xyz);

    // struct reb_orbit o = reb_tools_particle_to_orbit(ode->r->G, ode->r->particles[1], ode->r->particles[0]);
    struct reb_particle* const ps = ode->r->particles;

    double rvec[3] = {ps[1].x - ps[0].x, ps[1].y - ps[0].y, ps[1].z - ps[0].z};

    // a = o.a;
    year = pow(a, 1.5);
    n = 2 * M_PI / year;
    tidal_dt = atan(1 / Q_tide) / (2 * n);

    ts_cross(LHAT, rvec, l_cross_r);
    ts_cross(svec_xyz, rvec, s_cross_r);

    for (i=0;i<3;i++) {
        rhovec[i] = rvec[i] - (n * l_cross_r[i] * tidal_dt) + (s_cross_r[i] * tidal_dt);
    }

    prefac_rb = 3 * ode->r->G * ms / pow(a,5);
    prefac_tide = 3 * k2 * ode->r->G * pow(ms,2) * pow(rp,5) / pow(a,10);

    torque_rb[0] = prefac_rb * (K - J) * ts_dot(rvec, jhat) * ts_dot(rvec, khat);
    torque_rb[1] = prefac_rb * (I - K) * ts_dot(rvec, khat) * ts_dot(rvec, ihat);
    torque_rb[2] = prefac_rb * (J - I) * ts_dot(rvec, ihat) * ts_dot(rvec, jhat);
    ts_cross(rhovec, rvec, rho_cross_r);
    for (i=0;i<3;i++) {
        torque_tide_xyz[i] = prefac_tide * ts_dot(rvec, rhovec) * rho_cross_r[i];
    }
    matmul(inv_rot, torque_tide_xyz, torque_tide);
    for (i=0;i<3;i++) {
        torque[i] = torque_rb[i] + torque_tide[i];
    }

    double sidot = (torque[0] - (K - J) * sj * sk) / I;
    double sjdot = (torque[1] - (I - K) * sk * si) / J;
    double skdot = (torque[2] - (J - I) * si * sj) / K;

    double s_cross_i[3];
    double s_cross_j[3];
    double s_cross_k[3];
    ts_cross(svec_xyz, ihat, s_cross_i);
    ts_cross(svec_xyz, jhat, s_cross_j);
    ts_cross(svec_xyz, khat, s_cross_k);
    assign(yDot, s_cross_i, 3, 0, 0);
    assign(yDot, s_cross_j, 3, 3, 0);
    assign(yDot, s_cross_k, 3, 6, 0);
    double sdot[3] = {sidot, sjdot, skdot};
    assign(yDot, sdot, 3, 9, 0);
}

int main(int argc, char* argv[]) {

    // args: omega/n, theta (deg), final integration time, data filename
    if (argc != 5){
        printf("Incorrect number of command-line args.\n");
        exit(1);
    }
    // spin params
    double omega_by_n=atof(argv[1]);
    double obl0=radians(atof(argv[2]));
    double tf = atof(argv[3]); // in orbital periods
    char* data_fn = argv[4];
    double beta0=0.0;
    double ax0=0.0;

    // orbital params
    double a=0.4;
    double ms=1;
    double mp=4 * MEARTH;
    // tide params
    // integration params
    double year = pow(a, 1.5);
    double tol=1e-8; // integrator tolerance
    double reb_dt = 1e-2; // integrator dt, in years
    double out_dt=(int) tf / 1000; // in orbital periods; outputs 1000 data points
    double s0[3][3];
    double w0[3];

    double n = 2 * M_PI / year;
    int i, j;

    struct reb_simulation* r = reb_create_simulation();
    r->G = 39.476926421373; // to get units in AU, yr, M_sun

    reb_add_fmt(r, "m", ms);                // star
    reb_add_fmt(r, "m a e", mp, a, 0.0); // planet
    reb_move_to_com(r);

    r->integrator = REB_INTEGRATOR_BS;  // Bulirsch-Stoer integrator
    r->ri_bs.eps_rel = tol;            // Relative tolerance
    r->ri_bs.eps_abs = tol;            // Absolute tolerance
    r->dt = reb_dt;
    // r->ri_bs.max_dt = 10*reb_dt;


    struct reb_ode* ho = reb_create_ode(r,12);   // Add an ODE with 2 dimensions
    ho->derivatives = dydt_rb;              // Right hand side of the ODE
    struct reb_particle* const ps = r->particles;

    // initial conditions
    euler_to_body(0, obl0, ax0, s0);
    w0[0] = 0; w0[1]= sin(beta0)*n*omega_by_n; w0[2] = cos(beta0)*n*omega_by_n;
    assign(ho->y, s0[0], 3, 0, 0);
    assign(ho->y, s0[1], 3, 3, 0);
    assign(ho->y, s0[2], 3, 6, 0);
    assign(ho->y, w0, 3, 9, 0);

    int nt = round(tf / out_dt) + 1;
    int nv = 16;
    double* y_rb = malloc(nv*nt*sizeof(double)); // t, i[3], j[3], k[3], s[3], r[3]
    if (y_rb == NULL) {
        reb_free_ode(ho);
        reb_free_simulation(r);
        return 1;
    }
    clock_t start = clock();

    for (i=0; i<nt; i++) {
        reb_integrate(r, out_dt * i*year);
        // printf("Integrated %f, %f\n", r->t, ho->y[0]);
        // fflush(stdout);

        y_rb[i*nv] = r->t;
        for (j=1; j < 13; j++) {
            y_rb[i*nv + j] = ho->y[j-1];
        }
        y_rb[i*nv + 13] = ps[1].x - ps[0].x;
        y_rb[i*nv + 14] = ps[1].y - ps[0].y;
        y_rb[i*nv + 15] = ps[1].z - ps[0].z;
    }
    int tot_secs = (int) ((clock() - start) / CLOCKS_PER_SEC);
    int hours = (int) (tot_secs / 3600);
    int mins = (int) ((tot_secs % 3600) / 60);
    int secs = (tot_secs % 3600) % 60;
    printf("Integration took %dh %dm %ds\n", hours, mins, secs);

    FILE *f = fopen(data_fn, "wb");
    fwrite(y_rb, sizeof(double), nv*nt, f);
    fclose(f);

    free(y_rb);
    reb_free_ode(ho);
    reb_free_simulation(r);

    return 0;
}
