#include <math.h>
#include "rhs_values.c"

# define Power(x, y) (pow((double) (x), (double) (y)))

static void
get_dofs_g(DOF **dofs_Psi, DOF **dofs_g)
{
    INT i, n, data_count;
    FLOAT *p_Psi[6], *p_g[6];
    FLOAT Psi11, Psi12, Psi13, Psi22, Psi23, Psi33, det;
    FLOAT g11, g12, g13, g22, g23, g33;

    for(i=0; i<6; i++){
        p_Psi[i] = DofData(dofs_Psi[i + 4]);
        p_g[i] = DofData(dofs_g[i]);
    }

    data_count = DofGetDataCount(dofs_g[0]);
    for(n=0; n<data_count; n++){
        Psi11 = *(p_Psi[0]);
        Psi12 = *(p_Psi[1]);
        Psi13 = *(p_Psi[2]);
        Psi22 = *(p_Psi[3]);
        Psi23 = *(p_Psi[4]);
        Psi33 = *(p_Psi[5]);

        det = -(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 
                - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33;

        g11 = (-Power(Psi23,2) + Psi22*Psi33)/det;
        g12 = (Psi13*Psi23 - Psi12*Psi33)/det;
        g13 = (-(Psi13*Psi22) + Psi12*Psi23)/det;
        g22 = (-Power(Psi13,2) + Psi11*Psi33)/det;
        g23 = (Psi12*Psi13 - Psi11*Psi23)/det;
        g33 = (-Power(Psi12,2) + Psi11*Psi22)/det;

        *p_g[0] = g11;
        *p_g[1] = g12;
        *p_g[2] = g13;
        *p_g[3] = g22;
        *p_g[4] = g23;
        *p_g[5] = g33;

        for(i=0; i<6; i++){
            p_Psi[i]++;
            p_g[i]++;
        }
    } 
}


static void
get_dofs_N(DOF **dofs_Psi, DOF **dofs_g, DOF **dofs_N)
{
    INT i, n, data_count;
    FLOAT *p_g[6], *p_Psi[4], *p_N[4];
    FLOAT Psi00, Psi01, Psi02, Psi03;
    FLOAT g11, g12, g13, g22, g23, g33;
    FLOAT N, N1, N2, N3;

    for(i=0; i<4; i++){
        p_Psi[i] = DofData(dofs_Psi[i]);
        p_N[i] = DofData(dofs_N[i]);
    }

    for(i=0; i<6; i++){
        p_g[i] = DofData(dofs_g[i]);
    }

    data_count = DofGetDataCount(dofs_g[0]);
    for(n=0; n<data_count;n++){
        Psi00 = *(p_Psi[0]);
        Psi01 = *(p_Psi[1]);
        Psi02 = *(p_Psi[2]);
        Psi03 = *(p_Psi[3]);

        g11 = *(p_g[0]);
        g12 = *(p_g[1]);
        g13 = *(p_g[2]);
        g22 = *(p_g[3]);
        g23 = *(p_g[4]);
        g33 = *(p_g[5]);

        N = Power(-Psi00 + g11*Power(Psi01,2) + 2*g12*Psi01*Psi02 
            + g22*Power(Psi02,2) + 2*g13*Psi01*Psi03 + 2*g23*Psi02*Psi03 
            + g33*Power(Psi03,2),0.5);
        N1 = g11*Psi01 + g12*Psi02 + g13*Psi03;
        N2 = g12*Psi01 + g22*Psi02 + g23*Psi03;
        N3 = g13*Psi01 + g23*Psi02 + g33*Psi03;

        *p_N[0] = N;
        *p_N[1] = N1;
        *p_N[2] = N2;
        *p_N[3] = N3;

        for(i=0; i<4; i++){
            p_Psi[i]++;
            p_N[i]++;
        }
        for(i=0; i<6; i++){
            p_g[i]++;
        }
    }
}

static void
get_dofs_auxi(DOF **dofs_var, DOF **dofs_g, DOF **dofs_N, DOF **dofs_rhs)
{
    DOF **dofs_Psi = dofs_var;
    get_dofs_g(dofs_Psi, dofs_g);
    get_dofs_N(dofs_Psi, dofs_g, dofs_N);
    INT i, n, data_count;
    FLOAT *p_var[50], *p_rhs[50], values_var[50], values_rhs[50];
    for(INT i=0; i < 50; i++){
        p_var[i] = DofData(dofs_var[i]);
        p_rhs[i] = DofData(dofs_rhs[i]); 
    }

    //evaluate dofs values at every data point
    data_count = DofGetDataCount(dofs_var[0]);
    for(n=0;n<data_count;n++){
        //compute rhs values at data point
        for(i=0; i < 50; i++){
            values_var[i] = *(p_var[i]);
        }

        rhs_values(values_var, values_rhs); 

        for(i=0; i<50; i++){
            *p_rhs[i] = values_rhs[i];
            p_rhs[i]++;
            p_var[i]++; 
        }
    }
}
