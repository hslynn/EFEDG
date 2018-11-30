#include <math.h>

#define gamma0 1
#define gamma1 -1
#define gamma2 1

#define Power(x,y) (pow((double) (x), (double) (y)))

static void 
get_values_src(FLOAT *values_var, FLOAT *values)
{ 
    FLOAT invPsi00, invPsi01, invPsi02, invPsi03, invPsi11, invPsi12, invPsi13, invPsi22, invPsi23, invPsi33;
    FLOAT g11, g12, g13, g22, g23, g33;
    FLOAT N, N1, N2, N3;
    FLOAT t0, t1, t2, t3;
    FLOAT Gamma000,Gamma001,Gamma002,Gamma003,Gamma011,Gamma012,Gamma013,Gamma022,Gamma023,Gamma033;
    FLOAT Gamma100,Gamma101,Gamma102,Gamma103,Gamma111,Gamma112,Gamma113,Gamma122,Gamma123,Gamma133;
    FLOAT Gamma200,Gamma201,Gamma202,Gamma203,Gamma211,Gamma212,Gamma213,Gamma222,Gamma223,Gamma233;
    FLOAT Gamma300,Gamma301,Gamma302,Gamma303,Gamma311,Gamma312,Gamma313,Gamma322,Gamma323,Gamma333;
    FLOAT vecGamma0, vecGamma1, vecGamma2, vecGamma3;
    
    FLOAT *Psi = values_var, *Pi = values_var + 10, *Phi = values_var + 30; 
    FLOAT Psi00 = Psi[0], Psi01 = Psi[1], Psi02 = Psi[2], Psi03 = Psi[3], Psi11 = Psi[4];
    FLOAT Psi12 = Psi[5], Psi13 = Psi[6], Psi22 = Psi[7], Psi23 = Psi[8], Psi33 = Psi[9]; 
    
    FLOAT Pi00 = Pi[0], Pi01 = Pi[1], Pi02 = Pi[2], Pi03 = Pi[3], Pi11 = Pi[4];
    FLOAT Pi12 = Pi[5], Pi13 = Pi[6], Pi22 = Pi[7], Pi23 = Pi[8], Pi33 = Pi[9];
    
    FLOAT Phi000, Phi001, Phi002, Phi003, Phi011;//To be determined below
    FLOAT Phi012, Phi013, Phi022, Phi023, Phi033;

    FLOAT Phi100 = Phi[0], Phi101 = Phi[1], Phi102 = Phi[2], Phi103 = Phi[3], Phi111 = Phi[4];
    FLOAT Phi112 = Phi[5], Phi113 = Phi[6], Phi122 = Phi[7], Phi123 = Phi[8], Phi133 = Phi[9]; 
                                                                                     
    FLOAT Phi200 = Phi[10], Phi201 = Phi[11], Phi202 = Phi[12], Phi203 = Phi[13], Phi211 = Phi[14];
    FLOAT Phi212 = Phi[15], Phi213 = Phi[16], Phi222 = Phi[17], Phi223 = Phi[18], Phi233 = Phi[19];  
                                                                                     
    FLOAT Phi300 = Phi[20], Phi301 = Phi[21], Phi302 = Phi[22], Phi303 = Phi[23], Phi311 = Phi[24];
    FLOAT Phi312 = Phi[25], Phi313 = Phi[26], Phi322 = Phi[27], Phi323 = Phi[28], Phi333 = Phi[29]; 
    
    FLOAT srcPsi00,srcPsi01,srcPsi02,srcPsi03,srcPsi11,srcPsi12,srcPsi13,srcPsi22,srcPsi23,srcPsi33;
    
    FLOAT srcPi00, t1srcPi00, t2srcPi00, t3srcPi00, t4srcPi00, t5srcPi00, t6srcPi00;
    FLOAT srcPi01, t1srcPi01, t2srcPi01, t3srcPi01, t4srcPi01, t5srcPi01, t6srcPi01;
    FLOAT srcPi02, t1srcPi02, t2srcPi02, t3srcPi02, t4srcPi02, t5srcPi02, t6srcPi02;
    FLOAT srcPi03, t1srcPi03, t2srcPi03, t3srcPi03, t4srcPi03, t5srcPi03, t6srcPi03;
    FLOAT srcPi11, t1srcPi11, t2srcPi11, t3srcPi11, t4srcPi11, t5srcPi11, t6srcPi11;
    FLOAT srcPi12, t1srcPi12, t2srcPi12, t3srcPi12, t4srcPi12, t5srcPi12, t6srcPi12;
    FLOAT srcPi13, t1srcPi13, t2srcPi13, t3srcPi13, t4srcPi13, t5srcPi13, t6srcPi13;
    FLOAT srcPi22, t1srcPi22, t2srcPi22, t3srcPi22, t4srcPi22, t5srcPi22, t6srcPi22;
    FLOAT srcPi23, t1srcPi23, t2srcPi23, t3srcPi23, t4srcPi23, t5srcPi23, t6srcPi23;
    FLOAT srcPi33, t1srcPi33, t2srcPi33, t3srcPi33, t4srcPi33, t5srcPi33, t6srcPi33;
    
    FLOAT srcPhi100, t1srcPhi100, t2srcPhi100, t3srcPhi100;
    FLOAT srcPhi101, t1srcPhi101, t2srcPhi101, t3srcPhi101;
    FLOAT srcPhi102, t1srcPhi102, t2srcPhi102, t3srcPhi102;
    FLOAT srcPhi103, t1srcPhi103, t2srcPhi103, t3srcPhi103;
    FLOAT srcPhi111, t1srcPhi111, t2srcPhi111, t3srcPhi111;
    FLOAT srcPhi112, t1srcPhi112, t2srcPhi112, t3srcPhi112;
    FLOAT srcPhi113, t1srcPhi113, t2srcPhi113, t3srcPhi113;
    FLOAT srcPhi122, t1srcPhi122, t2srcPhi122, t3srcPhi122;
    FLOAT srcPhi123, t1srcPhi123, t2srcPhi123, t3srcPhi123;
    FLOAT srcPhi133, t1srcPhi133, t2srcPhi133, t3srcPhi133;
                                                   
    FLOAT srcPhi200, t1srcPhi200, t2srcPhi200, t3srcPhi200;
    FLOAT srcPhi201, t1srcPhi201, t2srcPhi201, t3srcPhi201;
    FLOAT srcPhi202, t1srcPhi202, t2srcPhi202, t3srcPhi202;
    FLOAT srcPhi203, t1srcPhi203, t2srcPhi203, t3srcPhi203;
    FLOAT srcPhi211, t1srcPhi211, t2srcPhi211, t3srcPhi211;
    FLOAT srcPhi212, t1srcPhi212, t2srcPhi212, t3srcPhi212;
    FLOAT srcPhi213, t1srcPhi213, t2srcPhi213, t3srcPhi213;
    FLOAT srcPhi222, t1srcPhi222, t2srcPhi222, t3srcPhi222;
    FLOAT srcPhi223, t1srcPhi223, t2srcPhi223, t3srcPhi223;
    FLOAT srcPhi233, t1srcPhi233, t2srcPhi233, t3srcPhi233;
                                                   
    FLOAT srcPhi300, t1srcPhi300, t2srcPhi300, t3srcPhi300;
    FLOAT srcPhi301, t1srcPhi301, t2srcPhi301, t3srcPhi301;
    FLOAT srcPhi302, t1srcPhi302, t2srcPhi302, t3srcPhi302;
    FLOAT srcPhi303, t1srcPhi303, t2srcPhi303, t3srcPhi303;
    FLOAT srcPhi311, t1srcPhi311, t2srcPhi311, t3srcPhi311;
    FLOAT srcPhi312, t1srcPhi312, t2srcPhi312, t3srcPhi312;
    FLOAT srcPhi313, t1srcPhi313, t2srcPhi313, t3srcPhi313;
    FLOAT srcPhi322, t1srcPhi322, t2srcPhi322, t3srcPhi322;
    FLOAT srcPhi323, t1srcPhi323, t2srcPhi323, t3srcPhi323;
    FLOAT srcPhi333, t1srcPhi333, t2srcPhi333, t3srcPhi333;
    
    //Gauge functions H^a
    FLOAT H0 = 0;
    FLOAT H1 = 0;
    FLOAT H2 = 0;
    FLOAT H3 = 0;
    
    FLOAT dH00 = 0;
    FLOAT dH01 = 0;
    FLOAT dH02 = 0;
    FLOAT dH03 = 0;
    
    FLOAT dH10 = 0;
    FLOAT dH11 = 0;
    FLOAT dH12 = 0;
    FLOAT dH13 = 0;
    
    FLOAT dH20 = 0;
    FLOAT dH21 = 0;
    FLOAT dH22 = 0;
    FLOAT dH23 = 0;
    
    FLOAT dH30 = 0;
    FLOAT dH31 = 0;
    FLOAT dH32 = 0;
    FLOAT dH33 = 0;
        
    //inverse spacetime metric
    invPsi00 = (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi01 = (Psi03*Psi13*Psi22 - Psi03*Psi12*Psi23 - Psi02*Psi13*Psi23 + Psi01*Power(Psi23,2) + 
           Psi02*Psi12*Psi33 - Psi01*Psi22*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi02 = (-(Psi03*Psi12*Psi13) + Psi02*Power(Psi13,2) + Psi03*Psi11*Psi23 - Psi01*Psi13*Psi23 - 
           Psi02*Psi11*Psi33 + Psi01*Psi12*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi03 = (Psi03*Power(Psi12,2) - Psi02*Psi12*Psi13 - Psi03*Psi11*Psi22 + Psi01*Psi13*Psi22 + 
           Psi02*Psi11*Psi23 - Psi01*Psi12*Psi23)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi11 = (-(Power(Psi03,2)*Psi22) + 2*Psi02*Psi03*Psi23 - Psi00*Power(Psi23,2) - Power(Psi02,2)*Psi33 + 
           Psi00*Psi22*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi12 = (Power(Psi03,2)*Psi12 - Psi02*Psi03*Psi13 - Psi01*Psi03*Psi23 + Psi00*Psi13*Psi23 + 
           Psi01*Psi02*Psi33 - Psi00*Psi12*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi13 = (-(Psi02*Psi03*Psi12) + Power(Psi02,2)*Psi13 + Psi01*Psi03*Psi22 - Psi00*Psi13*Psi22 - 
           Psi01*Psi02*Psi23 + Psi00*Psi12*Psi23)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi22 = (-(Power(Psi03,2)*Psi11) + 2*Psi01*Psi03*Psi13 - Psi00*Power(Psi13,2) - Power(Psi01,2)*Psi33 + 
           Psi00*Psi11*Psi33)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi23 = (Psi02*Psi03*Psi11 - Psi01*Psi03*Psi12 - Psi01*Psi02*Psi13 + Psi00*Psi12*Psi13 + 
           Power(Psi01,2)*Psi23 - Psi00*Psi11*Psi23)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
    
    invPsi33 = (-(Power(Psi02,2)*Psi11) + 2*Psi01*Psi02*Psi12 - Psi00*Power(Psi12,2) - Power(Psi01,2)*Psi22 + 
           Psi00*Psi11*Psi22)/
         (Power(Psi03,2)*Power(Psi12,2) - 2*Psi02*Psi03*Psi12*Psi13 + Power(Psi02,2)*Power(Psi13,2) - Power(Psi03,2)*Psi11*Psi22 + 
           2*Psi01*Psi03*Psi13*Psi22 - Psi00*Power(Psi13,2)*Psi22 + 2*Psi02*Psi03*Psi11*Psi23 - 2*Psi01*Psi03*Psi12*Psi23 - 
           2*Psi01*Psi02*Psi13*Psi23 + 2*Psi00*Psi12*Psi13*Psi23 + Power(Psi01,2)*Power(Psi23,2) - Psi00*Psi11*Power(Psi23,2) - 
           Power(Psi02,2)*Psi11*Psi33 + 2*Psi01*Psi02*Psi12*Psi33 - Psi00*Power(Psi12,2)*Psi33 - Power(Psi01,2)*Psi22*Psi33 + 
           Psi00*Psi11*Psi22*Psi33);
     
    //inverse space metric
    g11 = (-Power(Psi23,2) + Psi22*Psi33)/
            (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);

    g12 = (Psi13*Psi23 - Psi12*Psi33)/
         (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);
    
    g13 = (-(Psi13*Psi22) + Psi12*Psi23)/
         (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);
    
    g22 = (-Power(Psi13,2) + Psi11*Psi33)/
         (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);
    
    g23 = (Psi12*Psi13 - Psi11*Psi23)/
         (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);
    
    g33 = (-Power(Psi12,2) + Psi11*Psi22)/
         (-(Power(Psi13,2)*Psi22) + 2*Psi12*Psi13*Psi23 - Psi11*Power(Psi23,2) - Power(Psi12,2)*Psi33 + Psi11*Psi22*Psi33);
    
    //Shift
    N1 = g11*Psi01 + g12*Psi02 + g13*Psi03;
    N2 = g12*Psi01 + g22*Psi02 + g23*Psi03;
    N3 = g13*Psi01 + g23*Psi02 + g33*Psi03;
    
    //Lapse
    N = Power(-Psi00 + g11*Power(Psi01,2) + 2*g12*Psi01*Psi02 + g22*Power(Psi02,2) + 2*g13*Psi01*Psi03 + 2*g23*Psi02*Psi03 + 
        g33*Power(Psi03,2),0.5);
    
    //unit normal vector t^a
    t0 = 1/N;
    t1 = -N1/N;
    t2 = -N2/N;
    t3 = -N3/N;

    //time derivative of Psi, i.e, Phi0ab

    Phi000 = N1*Phi100 + N2*Phi200 + N3*Phi300 - N*Pi00;
    Phi001 = N1*Phi101 + N2*Phi201 + N3*Phi301 - N*Pi01;
    Phi002 = N1*Phi102 + N2*Phi202 + N3*Phi302 - N*Pi02;
    Phi003 = N1*Phi103 + N2*Phi203 + N3*Phi303 - N*Pi03;
    Phi011 = N1*Phi111 + N2*Phi211 + N3*Phi311 - N*Pi11;
    Phi012 = N1*Phi112 + N2*Phi212 + N3*Phi312 - N*Pi12;
    Phi013 = N1*Phi113 + N2*Phi213 + N3*Phi313 - N*Pi13;
    Phi022 = N1*Phi122 + N2*Phi222 + N3*Phi322 - N*Pi22;
    Phi023 = N1*Phi123 + N2*Phi223 + N3*Phi323 - N*Pi23;
    Phi033 = N1*Phi133 + N2*Phi233 + N3*Phi333 - N*Pi33;
         
    //Christoffel Symbles
    Gamma000 = 0.5*Phi000;
    Gamma001 = 0.5*Phi100;
    Gamma002 = 0.5*Phi200;
    Gamma003 = 0.5*Phi300;
    Gamma011 = 0.5*(-Phi011 + 2*Phi101);
    Gamma012 = 0.5*(-Phi012 + Phi102 + Phi201);
    Gamma013 = 0.5*(-Phi013 + Phi103 + Phi301);
    Gamma022 = 0.5*(-Phi022 + 2*Phi202);
    Gamma023 = 0.5*(-Phi023 + Phi203 + Phi302);
    Gamma033 = 0.5*(-Phi033 + 2*Phi303);
    
    Gamma100 = 0.5*(2*Phi001 - Phi100);
    Gamma101 = 0.5*Phi011;
    Gamma102 = 0.5*(Phi012 - Phi102 + Phi201);
    Gamma103 = 0.5*(Phi013 - Phi103 + Phi301);
    Gamma111 = 0.5*Phi111;
    Gamma112 = 0.5*Phi211;
    Gamma113 = 0.5*Phi311;
    Gamma122 = 0.5*(-Phi122 + 2*Phi212);
    Gamma123 = 0.5*(-Phi123 + Phi213 + Phi312);
    Gamma133 = 0.5*(-Phi133 + 2*Phi313);
    
    Gamma200 = 0.5*(2*Phi002 - Phi200);
    Gamma201 = 0.5*(Phi012 + Phi102 - Phi201);
    Gamma202 = 0.5*Phi022;
    Gamma203 = 0.5*(Phi023 - Phi203 + Phi302);
    Gamma211 = 0.5*(2*Phi112 - Phi211); 
    Gamma212 = 0.5*Phi122;
    Gamma213 = 0.5*(Phi123 - Phi213 + Phi312);
    Gamma222 = 0.5*Phi222;
    Gamma223 = 0.5*Phi322;
    Gamma233 = 0.5*(-Phi233 + 2*Phi323);
     
    Gamma300 = 0.5*(2*Phi003 - Phi300);
    Gamma301 = 0.5*(Phi013 + Phi103 - Phi301);
    Gamma302 = 0.5*(Phi023 + Phi203 - Phi302);
    Gamma303 = 0.5*Phi033;
    Gamma311 = 0.5*(2*Phi113 - Phi311);
    Gamma312 = 0.5*(Phi123 + Phi213 - Phi312);
    Gamma313 = 0.5*Phi133;
    Gamma322 = 0.5*(2*Phi223 - Phi322);
    Gamma323 = 0.5*Phi233;
    Gamma333 = 0.5*Phi333;
 
    //Gamma contracted with metric
    vecGamma0 = Gamma000*invPsi00 + Gamma100*invPsi00 + Gamma200*invPsi00 + Gamma300*invPsi00 + 2*Gamma001*invPsi01 + 
                 2*Gamma101*invPsi01 + 2*Gamma201*invPsi01 + 2*Gamma301*invPsi01 + 2*Gamma002*invPsi02 + 2*Gamma102*invPsi02 + 
                 2*Gamma202*invPsi02 + 2*Gamma302*invPsi02 + 2*Gamma003*invPsi03 + 2*Gamma103*invPsi03 + 2*Gamma203*invPsi03 + 
                 2*Gamma303*invPsi03 + Gamma011*invPsi11 + Gamma111*invPsi11 + Gamma211*invPsi11 + Gamma311*invPsi11 + 2*Gamma012*invPsi12 + 
                 2*Gamma112*invPsi12 + 2*Gamma212*invPsi12 + 2*Gamma312*invPsi12 + 2*Gamma013*invPsi13 + 2*Gamma113*invPsi13 + 
                 2*Gamma213*invPsi13 + 2*Gamma313*invPsi13 + Gamma022*invPsi22 + Gamma122*invPsi22 + Gamma222*invPsi22 + Gamma322*invPsi22 + 
                 2*Gamma023*invPsi23 + 2*Gamma123*invPsi23 + 2*Gamma223*invPsi23 + 2*Gamma323*invPsi23 + Gamma033*invPsi33 + 
                 Gamma133*invPsi33 + Gamma233*invPsi33 + Gamma333*invPsi33;
    vecGamma1 = Gamma000*invPsi00 + Gamma100*invPsi00 + Gamma200*invPsi00 + Gamma300*invPsi00 + 2*Gamma001*invPsi01 + 2*Gamma101*invPsi01 + 
                2*Gamma201*invPsi01 + 2*Gamma301*invPsi01 + 2*Gamma002*invPsi02 + 2*Gamma102*invPsi02 + 2*Gamma202*invPsi02 + 
                2*Gamma302*invPsi02 + 2*Gamma003*invPsi03 + 2*Gamma103*invPsi03 + 2*Gamma203*invPsi03 + 2*Gamma303*invPsi03 + 
                Gamma011*invPsi11 + Gamma111*invPsi11 + Gamma211*invPsi11 + Gamma311*invPsi11 + 2*Gamma012*invPsi12 + 2*Gamma112*invPsi12 + 
                2*Gamma212*invPsi12 + 2*Gamma312*invPsi12 + 2*Gamma013*invPsi13 + 2*Gamma113*invPsi13 + 2*Gamma213*invPsi13 + 
                2*Gamma313*invPsi13 + Gamma022*invPsi22 + Gamma122*invPsi22 + Gamma222*invPsi22 + Gamma322*invPsi22 + 2*Gamma023*invPsi23 + 
                2*Gamma123*invPsi23 + 2*Gamma223*invPsi23 + 2*Gamma323*invPsi23 + Gamma033*invPsi33 + Gamma133*invPsi33 + 
                Gamma233*invPsi33 + Gamma333*invPsi33;
    vecGamma2 = Gamma000*invPsi00 + Gamma100*invPsi00 + Gamma200*invPsi00 + Gamma300*invPsi00 + 
                 2*Gamma001*invPsi01 + 2*Gamma101*invPsi01 + 2*Gamma201*invPsi01 + 2*Gamma301*invPsi01 + 2*Gamma002*invPsi02 + 
                 2*Gamma102*invPsi02 + 2*Gamma202*invPsi02 + 2*Gamma302*invPsi02 + 2*Gamma003*invPsi03 + 2*Gamma103*invPsi03 + 
                 2*Gamma203*invPsi03 + 2*Gamma303*invPsi03 + Gamma011*invPsi11 + Gamma111*invPsi11 + Gamma211*invPsi11 + Gamma311*invPsi11 + 
                 2*Gamma012*invPsi12 + 2*Gamma112*invPsi12 + 2*Gamma212*invPsi12 + 2*Gamma312*invPsi12 + 2*Gamma013*invPsi13 + 
                 2*Gamma113*invPsi13 + 2*Gamma213*invPsi13 + 2*Gamma313*invPsi13 + Gamma022*invPsi22 + Gamma122*invPsi22 + 
                 Gamma222*invPsi22 + Gamma322*invPsi22 + 2*Gamma023*invPsi23 + 2*Gamma123*invPsi23 + 2*Gamma223*invPsi23 + 
                 2*Gamma323*invPsi23 + Gamma033*invPsi33 + Gamma133*invPsi33 + Gamma233*invPsi33 + Gamma333*invPsi33;
    vecGamma3 = Gamma000*invPsi00 + Gamma100*invPsi00 + Gamma200*invPsi00 + Gamma300*invPsi00 + 2*Gamma001*invPsi01 + 2*Gamma101*invPsi01 + 
                 2*Gamma201*invPsi01 + 2*Gamma301*invPsi01 + 2*Gamma002*invPsi02 + 2*Gamma102*invPsi02 + 2*Gamma202*invPsi02 + 
                 2*Gamma302*invPsi02 + 2*Gamma003*invPsi03 + 2*Gamma103*invPsi03 + 2*Gamma203*invPsi03 + 2*Gamma303*invPsi03 + 
                 Gamma011*invPsi11 + Gamma111*invPsi11 + Gamma211*invPsi11 + Gamma311*invPsi11 + 2*Gamma012*invPsi12 + 2*Gamma112*invPsi12 + 
                 2*Gamma212*invPsi12 + 2*Gamma312*invPsi12 + 2*Gamma013*invPsi13 + 2*Gamma113*invPsi13 + 2*Gamma213*invPsi13 + 
                 2*Gamma313*invPsi13 + Gamma022*invPsi22 + Gamma122*invPsi22 + Gamma222*invPsi22 + Gamma322*invPsi22 + 2*Gamma023*invPsi23 + 
                 2*Gamma123*invPsi23 + 2*Gamma223*invPsi23 + 2*Gamma323*invPsi23 + Gamma033*invPsi33 + Gamma133*invPsi33 + 
                 Gamma233*invPsi33 + Gamma333*invPsi33;

    
    //RHS of Psi
    srcPsi00 = (-gamma1)*(N1*Phi100 + N2*Phi200 + N3*Phi300) - N*Pi00; 
    srcPsi01 = (-gamma1)*(N1*Phi101 + N2*Phi201 + N3*Phi301) - N*Pi01; 
    srcPsi02 = (-gamma1)*(N1*Phi102 + N2*Phi202 + N3*Phi302) - N*Pi02; 
    srcPsi03 = (-gamma1)*(N1*Phi103 + N2*Phi203 + N3*Phi303) - N*Pi03; 
    srcPsi11 = (-gamma1)*(N1*Phi111 + N2*Phi211 + N3*Phi311) - N*Pi11; 
    srcPsi12 = (-gamma1)*(N1*Phi112 + N2*Phi212 + N3*Phi312) - N*Pi12;
    srcPsi13 = (-gamma1)*(N1*Phi113 + N2*Phi213 + N3*Phi313) - N*Pi13; 
    srcPsi22 = (-gamma1)*(N1*Phi122 + N2*Phi222 + N3*Phi322) - N*Pi22; 
    srcPsi23 = (-gamma1)*(N1*Phi123 + N2*Phi223 + N3*Phi323) - N*Pi23; 
    srcPsi33 = (-gamma1)*(N1*Phi133 + N2*Phi233 + N3*Phi333) - N*Pi33;
    
     
    //RHS of Pi
    
    //TERM 1
    t1srcPi00 = 2*N*(invPsi00*(-(Power(Gamma000,2)*invPsi00) - 2*Gamma000*Gamma001*invPsi01 - 2*Gamma000*Gamma002*invPsi02 - 
            Power(Gamma001,2)*invPsi11 - 2*Gamma001*Gamma002*invPsi12 - Power(Gamma002,2)*invPsi22 + g11*Power(Phi100,2) + 
            2*g12*Phi100*Phi200 + g22*Power(Phi200,2) + 2*g13*Phi100*Phi300 + 2*g23*Phi200*Phi300 + g33*Power(Phi300,2) - 
            Power(Pi00,2)) + 2*invPsi01*(-(Gamma000*Gamma001*invPsi00) - Power(Gamma001,2)*invPsi01 - Gamma000*Gamma011*invPsi01 - 
            Gamma001*Gamma002*invPsi02 - Gamma000*Gamma012*invPsi02 - Gamma001*Gamma011*invPsi11 - Gamma002*Gamma011*invPsi12 - 
            Gamma001*Gamma012*invPsi12 - Gamma002*Gamma012*invPsi22 + g11*Phi100*Phi101 + g12*Phi101*Phi200 + g12*Phi100*Phi201 + 
            g22*Phi200*Phi201 + g13*Phi101*Phi300 + g23*Phi201*Phi300 + g13*Phi100*Phi301 + g23*Phi200*Phi301 + g33*Phi300*Phi301 - 
            Pi00*Pi01) + invPsi11*(-(Power(Gamma001,2)*invPsi00) - 2*Gamma001*Gamma011*invPsi01 - 2*Gamma001*Gamma012*invPsi02 - 
            Power(Gamma011,2)*invPsi11 - 2*Gamma011*Gamma012*invPsi12 - Power(Gamma012,2)*invPsi22 + g11*Power(Phi101,2) + 
            2*g12*Phi101*Phi201 + g22*Power(Phi201,2) + 2*g13*Phi101*Phi301 + 2*g23*Phi201*Phi301 + g33*Power(Phi301,2) - 
            Power(Pi01,2)) + 2*invPsi02*(-(Gamma000*Gamma002*invPsi00) - Gamma001*Gamma002*invPsi01 - Gamma000*Gamma012*invPsi01 - 
            Power(Gamma002,2)*invPsi02 - Gamma000*Gamma022*invPsi02 - Gamma001*Gamma012*invPsi11 - Gamma002*Gamma012*invPsi12 - 
            Gamma001*Gamma022*invPsi12 - Gamma002*Gamma022*invPsi22 + g11*Phi100*Phi102 + g12*Phi102*Phi200 + g12*Phi100*Phi202 + 
            g22*Phi200*Phi202 + g13*Phi102*Phi300 + g23*Phi202*Phi300 + g13*Phi100*Phi302 + g23*Phi200*Phi302 + g33*Phi300*Phi302 - 
            Pi00*Pi02) + 2*invPsi12*(-(Gamma001*Gamma002*invPsi00) - Gamma002*Gamma011*invPsi01 - Gamma001*Gamma012*invPsi01 - 
            Gamma002*Gamma012*invPsi02 - Gamma001*Gamma022*invPsi02 - Gamma011*Gamma012*invPsi11 - Power(Gamma012,2)*invPsi12 - 
            Gamma011*Gamma022*invPsi12 - Gamma012*Gamma022*invPsi22 + g11*Phi101*Phi102 + g12*Phi102*Phi201 + g12*Phi101*Phi202 + 
            g22*Phi201*Phi202 + g13*Phi102*Phi301 + g23*Phi202*Phi301 + g13*Phi101*Phi302 + g23*Phi201*Phi302 + g33*Phi301*Phi302 - 
            Pi01*Pi02) + invPsi22*(-(Power(Gamma002,2)*invPsi00) - 2*Gamma002*Gamma012*invPsi01 - 2*Gamma002*Gamma022*invPsi02 - 
            Power(Gamma012,2)*invPsi11 - 2*Gamma012*Gamma022*invPsi12 - Power(Gamma022,2)*invPsi22 + g11*Power(Phi102,2) + 
            2*g12*Phi102*Phi202 + g22*Power(Phi202,2) + 2*g13*Phi102*Phi302 + 2*g23*Phi202*Phi302 + g33*Power(Phi302,2) - 
            Power(Pi02,2)) + 2*invPsi03*(-(Gamma000*Gamma003*invPsi00) - Gamma001*Gamma003*invPsi01 - Gamma000*Gamma013*invPsi01 - 
            Gamma002*Gamma003*invPsi02 - Gamma000*Gamma023*invPsi02 - Gamma001*Gamma013*invPsi11 - Gamma002*Gamma013*invPsi12 - 
            Gamma001*Gamma023*invPsi12 - Gamma002*Gamma023*invPsi22 + g11*Phi100*Phi103 + g12*Phi103*Phi200 + g12*Phi100*Phi203 + 
            g22*Phi200*Phi203 + g13*Phi103*Phi300 + g23*Phi203*Phi300 + g13*Phi100*Phi303 + g23*Phi200*Phi303 + g33*Phi300*Phi303 - 
            Pi00*Pi03) + 2*invPsi13*(-(Gamma001*Gamma003*invPsi00) - Gamma003*Gamma011*invPsi01 - Gamma001*Gamma013*invPsi01 - 
            Gamma003*Gamma012*invPsi02 - Gamma001*Gamma023*invPsi02 - Gamma011*Gamma013*invPsi11 - Gamma012*Gamma013*invPsi12 - 
            Gamma011*Gamma023*invPsi12 - Gamma012*Gamma023*invPsi22 + g11*Phi101*Phi103 + g12*Phi103*Phi201 + g12*Phi101*Phi203 + 
            g22*Phi201*Phi203 + g13*Phi103*Phi301 + g23*Phi203*Phi301 + g13*Phi101*Phi303 + g23*Phi201*Phi303 + g33*Phi301*Phi303 - 
            Pi01*Pi03) + 2*invPsi23*(-(Gamma002*Gamma003*invPsi00) - Gamma003*Gamma012*invPsi01 - Gamma002*Gamma013*invPsi01 - 
            Gamma003*Gamma022*invPsi02 - Gamma002*Gamma023*invPsi02 - Gamma012*Gamma013*invPsi11 - Gamma013*Gamma022*invPsi12 - 
            Gamma012*Gamma023*invPsi12 - Gamma022*Gamma023*invPsi22 + g11*Phi102*Phi103 + g12*Phi103*Phi202 + g12*Phi102*Phi203 + 
            g22*Phi202*Phi203 + g13*Phi103*Phi302 + g23*Phi203*Phi302 + g13*Phi102*Phi303 + g23*Phi202*Phi303 + g33*Phi302*Phi303 - 
            Pi02*Pi03) + invPsi33*(-(Power(Gamma003,2)*invPsi00) - 2*Gamma003*Gamma013*invPsi01 - 2*Gamma003*Gamma023*invPsi02 - 
            Power(Gamma013,2)*invPsi11 - 2*Gamma013*Gamma023*invPsi12 - Power(Gamma023,2)*invPsi22 + g11*Power(Phi103,2) + 
            2*g12*Phi103*Phi203 + g22*Power(Phi203,2) + 2*g13*Phi103*Phi303 + 2*g23*Phi203*Phi303 + g33*Power(Phi303,2) - 
            Power(Pi03,2)));
    
    t1srcPi01 = 2*N*(invPsi00*(-(Gamma000*Gamma100*invPsi00) - Gamma001*Gamma100*invPsi01 - Gamma000*Gamma101*invPsi01 - 
            Gamma002*Gamma100*invPsi02 - Gamma000*Gamma102*invPsi02 - Gamma001*Gamma101*invPsi11 - Gamma002*Gamma101*invPsi12 - 
            Gamma001*Gamma102*invPsi12 - Gamma002*Gamma102*invPsi22 + g11*Phi100*Phi101 + g12*Phi101*Phi200 + g12*Phi100*Phi201 + 
            g22*Phi200*Phi201 + g13*Phi101*Phi300 + g23*Phi201*Phi300 + g13*Phi100*Phi301 + g23*Phi200*Phi301 + g33*Phi300*Phi301 - 
            Pi00*Pi01) + invPsi01*(-(Gamma001*Gamma100*invPsi00) - Gamma011*Gamma100*invPsi01 - Gamma001*Gamma101*invPsi01 - 
            Gamma012*Gamma100*invPsi02 - Gamma001*Gamma102*invPsi02 - Gamma011*Gamma101*invPsi11 - Gamma012*Gamma101*invPsi12 - 
            Gamma011*Gamma102*invPsi12 - Gamma012*Gamma102*invPsi22 + g11*Power(Phi101,2) + 2*g12*Phi101*Phi201 + 
            g22*Power(Phi201,2) + 2*g13*Phi101*Phi301 + 2*g23*Phi201*Phi301 + g33*Power(Phi301,2) - Power(Pi01,2)) + 
         invPsi02*(-(Gamma002*Gamma100*invPsi00) - Gamma012*Gamma100*invPsi01 - Gamma002*Gamma101*invPsi01 - 
            Gamma022*Gamma100*invPsi02 - Gamma002*Gamma102*invPsi02 - Gamma012*Gamma101*invPsi11 - Gamma022*Gamma101*invPsi12 - 
            Gamma012*Gamma102*invPsi12 - Gamma022*Gamma102*invPsi22 + g11*Phi101*Phi102 + g12*Phi102*Phi201 + g12*Phi101*Phi202 + 
            g22*Phi201*Phi202 + g13*Phi102*Phi301 + g23*Phi202*Phi301 + g13*Phi101*Phi302 + g23*Phi201*Phi302 + g33*Phi301*Phi302 - 
            Pi01*Pi02) + invPsi03*(-(Gamma003*Gamma100*invPsi00) - Gamma013*Gamma100*invPsi01 - Gamma003*Gamma101*invPsi01 - 
            Gamma023*Gamma100*invPsi02 - Gamma003*Gamma102*invPsi02 - Gamma013*Gamma101*invPsi11 - Gamma023*Gamma101*invPsi12 - 
            Gamma013*Gamma102*invPsi12 - Gamma023*Gamma102*invPsi22 + g11*Phi101*Phi103 + g12*Phi103*Phi201 + g12*Phi101*Phi203 + 
            g22*Phi201*Phi203 + g13*Phi103*Phi301 + g23*Phi203*Phi301 + g13*Phi101*Phi303 + g23*Phi201*Phi303 + g33*Phi301*Phi303 - 
            Pi01*Pi03) + invPsi01*(-(Gamma000*Gamma101*invPsi00) - Gamma001*Gamma101*invPsi01 - Gamma000*Gamma111*invPsi01 - 
            Gamma002*Gamma101*invPsi02 - Gamma000*Gamma112*invPsi02 - Gamma001*Gamma111*invPsi11 - Gamma002*Gamma111*invPsi12 - 
            Gamma001*Gamma112*invPsi12 - Gamma002*Gamma112*invPsi22 + g11*Phi100*Phi111 + g12*Phi111*Phi200 + g12*Phi100*Phi211 + 
            g22*Phi200*Phi211 + g13*Phi111*Phi300 + g23*Phi211*Phi300 + g13*Phi100*Phi311 + g23*Phi200*Phi311 + g33*Phi300*Phi311 - 
            Pi00*Pi11) + invPsi11*(-(Gamma001*Gamma101*invPsi00) - Gamma011*Gamma101*invPsi01 - Gamma001*Gamma111*invPsi01 - 
            Gamma012*Gamma101*invPsi02 - Gamma001*Gamma112*invPsi02 - Gamma011*Gamma111*invPsi11 - Gamma012*Gamma111*invPsi12 - 
            Gamma011*Gamma112*invPsi12 - Gamma012*Gamma112*invPsi22 + g11*Phi101*Phi111 + g12*Phi111*Phi201 + g12*Phi101*Phi211 + 
            g22*Phi201*Phi211 + g13*Phi111*Phi301 + g23*Phi211*Phi301 + g13*Phi101*Phi311 + g23*Phi201*Phi311 + g33*Phi301*Phi311 - 
            Pi01*Pi11) + invPsi12*(-(Gamma002*Gamma101*invPsi00) - Gamma012*Gamma101*invPsi01 - Gamma002*Gamma111*invPsi01 - 
            Gamma022*Gamma101*invPsi02 - Gamma002*Gamma112*invPsi02 - Gamma012*Gamma111*invPsi11 - Gamma022*Gamma111*invPsi12 - 
            Gamma012*Gamma112*invPsi12 - Gamma022*Gamma112*invPsi22 + g11*Phi102*Phi111 + g12*Phi111*Phi202 + g12*Phi102*Phi211 + 
            g22*Phi202*Phi211 + g13*Phi111*Phi302 + g23*Phi211*Phi302 + g13*Phi102*Phi311 + g23*Phi202*Phi311 + g33*Phi302*Phi311 - 
            Pi02*Pi11) + invPsi13*(-(Gamma003*Gamma101*invPsi00) - Gamma013*Gamma101*invPsi01 - Gamma003*Gamma111*invPsi01 - 
            Gamma023*Gamma101*invPsi02 - Gamma003*Gamma112*invPsi02 - Gamma013*Gamma111*invPsi11 - Gamma023*Gamma111*invPsi12 - 
            Gamma013*Gamma112*invPsi12 - Gamma023*Gamma112*invPsi22 + g11*Phi103*Phi111 + g12*Phi111*Phi203 + g12*Phi103*Phi211 + 
            g22*Phi203*Phi211 + g13*Phi111*Phi303 + g23*Phi211*Phi303 + g13*Phi103*Phi311 + g23*Phi203*Phi311 + g33*Phi303*Phi311 - 
            Pi03*Pi11) + invPsi02*(-(Gamma000*Gamma102*invPsi00) - Gamma001*Gamma102*invPsi01 - Gamma000*Gamma112*invPsi01 - 
            Gamma002*Gamma102*invPsi02 - Gamma000*Gamma122*invPsi02 - Gamma001*Gamma112*invPsi11 - Gamma002*Gamma112*invPsi12 - 
            Gamma001*Gamma122*invPsi12 - Gamma002*Gamma122*invPsi22 + g11*Phi100*Phi112 + g12*Phi112*Phi200 + g12*Phi100*Phi212 + 
            g22*Phi200*Phi212 + g13*Phi112*Phi300 + g23*Phi212*Phi300 + g13*Phi100*Phi312 + g23*Phi200*Phi312 + g33*Phi300*Phi312 - 
            Pi00*Pi12) + invPsi12*(-(Gamma001*Gamma102*invPsi00) - Gamma011*Gamma102*invPsi01 - Gamma001*Gamma112*invPsi01 - 
            Gamma012*Gamma102*invPsi02 - Gamma001*Gamma122*invPsi02 - Gamma011*Gamma112*invPsi11 - Gamma012*Gamma112*invPsi12 - 
            Gamma011*Gamma122*invPsi12 - Gamma012*Gamma122*invPsi22 + g11*Phi101*Phi112 + g12*Phi112*Phi201 + g12*Phi101*Phi212 + 
            g22*Phi201*Phi212 + g13*Phi112*Phi301 + g23*Phi212*Phi301 + g13*Phi101*Phi312 + g23*Phi201*Phi312 + g33*Phi301*Phi312 - 
            Pi01*Pi12) + invPsi22*(-(Gamma002*Gamma102*invPsi00) - Gamma012*Gamma102*invPsi01 - Gamma002*Gamma112*invPsi01 - 
            Gamma022*Gamma102*invPsi02 - Gamma002*Gamma122*invPsi02 - Gamma012*Gamma112*invPsi11 - Gamma022*Gamma112*invPsi12 - 
            Gamma012*Gamma122*invPsi12 - Gamma022*Gamma122*invPsi22 + g11*Phi102*Phi112 + g12*Phi112*Phi202 + g12*Phi102*Phi212 + 
            g22*Phi202*Phi212 + g13*Phi112*Phi302 + g23*Phi212*Phi302 + g13*Phi102*Phi312 + g23*Phi202*Phi312 + g33*Phi302*Phi312 - 
            Pi02*Pi12) + invPsi23*(-(Gamma003*Gamma102*invPsi00) - Gamma013*Gamma102*invPsi01 - Gamma003*Gamma112*invPsi01 - 
            Gamma023*Gamma102*invPsi02 - Gamma003*Gamma122*invPsi02 - Gamma013*Gamma112*invPsi11 - Gamma023*Gamma112*invPsi12 - 
            Gamma013*Gamma122*invPsi12 - Gamma023*Gamma122*invPsi22 + g11*Phi103*Phi112 + g12*Phi112*Phi203 + g12*Phi103*Phi212 + 
            g22*Phi203*Phi212 + g13*Phi112*Phi303 + g23*Phi212*Phi303 + g13*Phi103*Phi312 + g23*Phi203*Phi312 + g33*Phi303*Phi312 - 
            Pi03*Pi12) + invPsi03*(-(Gamma000*Gamma103*invPsi00) - Gamma001*Gamma103*invPsi01 - Gamma000*Gamma113*invPsi01 - 
            Gamma002*Gamma103*invPsi02 - Gamma000*Gamma123*invPsi02 - Gamma001*Gamma113*invPsi11 - Gamma002*Gamma113*invPsi12 - 
            Gamma001*Gamma123*invPsi12 - Gamma002*Gamma123*invPsi22 + g11*Phi100*Phi113 + g12*Phi113*Phi200 + g12*Phi100*Phi213 + 
            g22*Phi200*Phi213 + g13*Phi113*Phi300 + g23*Phi213*Phi300 + g13*Phi100*Phi313 + g23*Phi200*Phi313 + g33*Phi300*Phi313 - 
            Pi00*Pi13) + invPsi13*(-(Gamma001*Gamma103*invPsi00) - Gamma011*Gamma103*invPsi01 - Gamma001*Gamma113*invPsi01 - 
            Gamma012*Gamma103*invPsi02 - Gamma001*Gamma123*invPsi02 - Gamma011*Gamma113*invPsi11 - Gamma012*Gamma113*invPsi12 - 
            Gamma011*Gamma123*invPsi12 - Gamma012*Gamma123*invPsi22 + g11*Phi101*Phi113 + g12*Phi113*Phi201 + g12*Phi101*Phi213 + 
            g22*Phi201*Phi213 + g13*Phi113*Phi301 + g23*Phi213*Phi301 + g13*Phi101*Phi313 + g23*Phi201*Phi313 + g33*Phi301*Phi313 - 
            Pi01*Pi13) + invPsi23*(-(Gamma002*Gamma103*invPsi00) - Gamma012*Gamma103*invPsi01 - Gamma002*Gamma113*invPsi01 - 
            Gamma022*Gamma103*invPsi02 - Gamma002*Gamma123*invPsi02 - Gamma012*Gamma113*invPsi11 - Gamma022*Gamma113*invPsi12 - 
            Gamma012*Gamma123*invPsi12 - Gamma022*Gamma123*invPsi22 + g11*Phi102*Phi113 + g12*Phi113*Phi202 + g12*Phi102*Phi213 + 
            g22*Phi202*Phi213 + g13*Phi113*Phi302 + g23*Phi213*Phi302 + g13*Phi102*Phi313 + g23*Phi202*Phi313 + g33*Phi302*Phi313 - 
            Pi02*Pi13) + invPsi33*(-(Gamma003*Gamma103*invPsi00) - Gamma013*Gamma103*invPsi01 - Gamma003*Gamma113*invPsi01 - 
            Gamma023*Gamma103*invPsi02 - Gamma003*Gamma123*invPsi02 - Gamma013*Gamma113*invPsi11 - Gamma023*Gamma113*invPsi12 - 
            Gamma013*Gamma123*invPsi12 - Gamma023*Gamma123*invPsi22 + g11*Phi103*Phi113 + g12*Phi113*Phi203 + g12*Phi103*Phi213 + 
            g22*Phi203*Phi213 + g13*Phi113*Phi303 + g23*Phi213*Phi303 + g13*Phi103*Phi313 + g23*Phi203*Phi313 + g33*Phi303*Phi313 - 
            Pi03*Pi13));
    
    t1srcPi02 = 2*N*(invPsi00*(-(Gamma000*Gamma200*invPsi00) - Gamma001*Gamma200*invPsi01 - Gamma000*Gamma201*invPsi01 - 
            Gamma002*Gamma200*invPsi02 - Gamma000*Gamma202*invPsi02 - Gamma001*Gamma201*invPsi11 - Gamma002*Gamma201*invPsi12 - 
            Gamma001*Gamma202*invPsi12 - Gamma002*Gamma202*invPsi22 + g11*Phi100*Phi102 + g12*Phi102*Phi200 + g12*Phi100*Phi202 + 
            g22*Phi200*Phi202 + g13*Phi102*Phi300 + g23*Phi202*Phi300 + g13*Phi100*Phi302 + g23*Phi200*Phi302 + g33*Phi300*Phi302 - 
            Pi00*Pi02) + invPsi01*(-(Gamma001*Gamma200*invPsi00) - Gamma011*Gamma200*invPsi01 - Gamma001*Gamma201*invPsi01 - 
            Gamma012*Gamma200*invPsi02 - Gamma001*Gamma202*invPsi02 - Gamma011*Gamma201*invPsi11 - Gamma012*Gamma201*invPsi12 - 
            Gamma011*Gamma202*invPsi12 - Gamma012*Gamma202*invPsi22 + g11*Phi101*Phi102 + g12*Phi102*Phi201 + g12*Phi101*Phi202 + 
            g22*Phi201*Phi202 + g13*Phi102*Phi301 + g23*Phi202*Phi301 + g13*Phi101*Phi302 + g23*Phi201*Phi302 + g33*Phi301*Phi302 - 
            Pi01*Pi02) + invPsi02*(-(Gamma002*Gamma200*invPsi00) - Gamma012*Gamma200*invPsi01 - Gamma002*Gamma201*invPsi01 - 
            Gamma022*Gamma200*invPsi02 - Gamma002*Gamma202*invPsi02 - Gamma012*Gamma201*invPsi11 - Gamma022*Gamma201*invPsi12 - 
            Gamma012*Gamma202*invPsi12 - Gamma022*Gamma202*invPsi22 + g11*Power(Phi102,2) + 2*g12*Phi102*Phi202 + 
            g22*Power(Phi202,2) + 2*g13*Phi102*Phi302 + 2*g23*Phi202*Phi302 + g33*Power(Phi302,2) - Power(Pi02,2)) + 
         invPsi03*(-(Gamma003*Gamma200*invPsi00) - Gamma013*Gamma200*invPsi01 - Gamma003*Gamma201*invPsi01 - 
            Gamma023*Gamma200*invPsi02 - Gamma003*Gamma202*invPsi02 - Gamma013*Gamma201*invPsi11 - Gamma023*Gamma201*invPsi12 - 
            Gamma013*Gamma202*invPsi12 - Gamma023*Gamma202*invPsi22 + g11*Phi102*Phi103 + g12*Phi103*Phi202 + g12*Phi102*Phi203 + 
            g22*Phi202*Phi203 + g13*Phi103*Phi302 + g23*Phi203*Phi302 + g13*Phi102*Phi303 + g23*Phi202*Phi303 + g33*Phi302*Phi303 - 
            Pi02*Pi03) + invPsi01*(-(Gamma000*Gamma201*invPsi00) - Gamma001*Gamma201*invPsi01 - Gamma000*Gamma211*invPsi01 - 
            Gamma002*Gamma201*invPsi02 - Gamma000*Gamma212*invPsi02 - Gamma001*Gamma211*invPsi11 - Gamma002*Gamma211*invPsi12 - 
            Gamma001*Gamma212*invPsi12 - Gamma002*Gamma212*invPsi22 + g11*Phi100*Phi112 + g12*Phi112*Phi200 + g12*Phi100*Phi212 + 
            g22*Phi200*Phi212 + g13*Phi112*Phi300 + g23*Phi212*Phi300 + g13*Phi100*Phi312 + g23*Phi200*Phi312 + g33*Phi300*Phi312 - 
            Pi00*Pi12) + invPsi11*(-(Gamma001*Gamma201*invPsi00) - Gamma011*Gamma201*invPsi01 - Gamma001*Gamma211*invPsi01 - 
            Gamma012*Gamma201*invPsi02 - Gamma001*Gamma212*invPsi02 - Gamma011*Gamma211*invPsi11 - Gamma012*Gamma211*invPsi12 - 
            Gamma011*Gamma212*invPsi12 - Gamma012*Gamma212*invPsi22 + g11*Phi101*Phi112 + g12*Phi112*Phi201 + g12*Phi101*Phi212 + 
            g22*Phi201*Phi212 + g13*Phi112*Phi301 + g23*Phi212*Phi301 + g13*Phi101*Phi312 + g23*Phi201*Phi312 + g33*Phi301*Phi312 - 
            Pi01*Pi12) + invPsi12*(-(Gamma002*Gamma201*invPsi00) - Gamma012*Gamma201*invPsi01 - Gamma002*Gamma211*invPsi01 - 
            Gamma022*Gamma201*invPsi02 - Gamma002*Gamma212*invPsi02 - Gamma012*Gamma211*invPsi11 - Gamma022*Gamma211*invPsi12 - 
            Gamma012*Gamma212*invPsi12 - Gamma022*Gamma212*invPsi22 + g11*Phi102*Phi112 + g12*Phi112*Phi202 + g12*Phi102*Phi212 + 
            g22*Phi202*Phi212 + g13*Phi112*Phi302 + g23*Phi212*Phi302 + g13*Phi102*Phi312 + g23*Phi202*Phi312 + g33*Phi302*Phi312 - 
            Pi02*Pi12) + invPsi13*(-(Gamma003*Gamma201*invPsi00) - Gamma013*Gamma201*invPsi01 - Gamma003*Gamma211*invPsi01 - 
            Gamma023*Gamma201*invPsi02 - Gamma003*Gamma212*invPsi02 - Gamma013*Gamma211*invPsi11 - Gamma023*Gamma211*invPsi12 - 
            Gamma013*Gamma212*invPsi12 - Gamma023*Gamma212*invPsi22 + g11*Phi103*Phi112 + g12*Phi112*Phi203 + g12*Phi103*Phi212 + 
            g22*Phi203*Phi212 + g13*Phi112*Phi303 + g23*Phi212*Phi303 + g13*Phi103*Phi312 + g23*Phi203*Phi312 + g33*Phi303*Phi312 - 
            Pi03*Pi12) + invPsi02*(-(Gamma000*Gamma202*invPsi00) - Gamma001*Gamma202*invPsi01 - Gamma000*Gamma212*invPsi01 - 
            Gamma002*Gamma202*invPsi02 - Gamma000*Gamma222*invPsi02 - Gamma001*Gamma212*invPsi11 - Gamma002*Gamma212*invPsi12 - 
            Gamma001*Gamma222*invPsi12 - Gamma002*Gamma222*invPsi22 + g11*Phi100*Phi122 + g12*Phi122*Phi200 + g12*Phi100*Phi222 + 
            g22*Phi200*Phi222 + g13*Phi122*Phi300 + g23*Phi222*Phi300 + g13*Phi100*Phi322 + g23*Phi200*Phi322 + g33*Phi300*Phi322 - 
            Pi00*Pi22) + invPsi12*(-(Gamma001*Gamma202*invPsi00) - Gamma011*Gamma202*invPsi01 - Gamma001*Gamma212*invPsi01 - 
            Gamma012*Gamma202*invPsi02 - Gamma001*Gamma222*invPsi02 - Gamma011*Gamma212*invPsi11 - Gamma012*Gamma212*invPsi12 - 
            Gamma011*Gamma222*invPsi12 - Gamma012*Gamma222*invPsi22 + g11*Phi101*Phi122 + g12*Phi122*Phi201 + g12*Phi101*Phi222 + 
            g22*Phi201*Phi222 + g13*Phi122*Phi301 + g23*Phi222*Phi301 + g13*Phi101*Phi322 + g23*Phi201*Phi322 + g33*Phi301*Phi322 - 
            Pi01*Pi22) + invPsi22*(-(Gamma002*Gamma202*invPsi00) - Gamma012*Gamma202*invPsi01 - Gamma002*Gamma212*invPsi01 - 
            Gamma022*Gamma202*invPsi02 - Gamma002*Gamma222*invPsi02 - Gamma012*Gamma212*invPsi11 - Gamma022*Gamma212*invPsi12 - 
            Gamma012*Gamma222*invPsi12 - Gamma022*Gamma222*invPsi22 + g11*Phi102*Phi122 + g12*Phi122*Phi202 + g12*Phi102*Phi222 + 
            g22*Phi202*Phi222 + g13*Phi122*Phi302 + g23*Phi222*Phi302 + g13*Phi102*Phi322 + g23*Phi202*Phi322 + g33*Phi302*Phi322 - 
            Pi02*Pi22) + invPsi23*(-(Gamma003*Gamma202*invPsi00) - Gamma013*Gamma202*invPsi01 - Gamma003*Gamma212*invPsi01 - 
            Gamma023*Gamma202*invPsi02 - Gamma003*Gamma222*invPsi02 - Gamma013*Gamma212*invPsi11 - Gamma023*Gamma212*invPsi12 - 
            Gamma013*Gamma222*invPsi12 - Gamma023*Gamma222*invPsi22 + g11*Phi103*Phi122 + g12*Phi122*Phi203 + g12*Phi103*Phi222 + 
            g22*Phi203*Phi222 + g13*Phi122*Phi303 + g23*Phi222*Phi303 + g13*Phi103*Phi322 + g23*Phi203*Phi322 + g33*Phi303*Phi322 - 
            Pi03*Pi22) + invPsi03*(-(Gamma000*Gamma203*invPsi00) - Gamma001*Gamma203*invPsi01 - Gamma000*Gamma213*invPsi01 - 
            Gamma002*Gamma203*invPsi02 - Gamma000*Gamma223*invPsi02 - Gamma001*Gamma213*invPsi11 - Gamma002*Gamma213*invPsi12 - 
            Gamma001*Gamma223*invPsi12 - Gamma002*Gamma223*invPsi22 + g11*Phi100*Phi123 + g12*Phi123*Phi200 + g12*Phi100*Phi223 + 
            g22*Phi200*Phi223 + g13*Phi123*Phi300 + g23*Phi223*Phi300 + g13*Phi100*Phi323 + g23*Phi200*Phi323 + g33*Phi300*Phi323 - 
            Pi00*Pi23) + invPsi13*(-(Gamma001*Gamma203*invPsi00) - Gamma011*Gamma203*invPsi01 - Gamma001*Gamma213*invPsi01 - 
            Gamma012*Gamma203*invPsi02 - Gamma001*Gamma223*invPsi02 - Gamma011*Gamma213*invPsi11 - Gamma012*Gamma213*invPsi12 - 
            Gamma011*Gamma223*invPsi12 - Gamma012*Gamma223*invPsi22 + g11*Phi101*Phi123 + g12*Phi123*Phi201 + g12*Phi101*Phi223 + 
            g22*Phi201*Phi223 + g13*Phi123*Phi301 + g23*Phi223*Phi301 + g13*Phi101*Phi323 + g23*Phi201*Phi323 + g33*Phi301*Phi323 - 
            Pi01*Pi23) + invPsi23*(-(Gamma002*Gamma203*invPsi00) - Gamma012*Gamma203*invPsi01 - Gamma002*Gamma213*invPsi01 - 
            Gamma022*Gamma203*invPsi02 - Gamma002*Gamma223*invPsi02 - Gamma012*Gamma213*invPsi11 - Gamma022*Gamma213*invPsi12 - 
            Gamma012*Gamma223*invPsi12 - Gamma022*Gamma223*invPsi22 + g11*Phi102*Phi123 + g12*Phi123*Phi202 + g12*Phi102*Phi223 + 
            g22*Phi202*Phi223 + g13*Phi123*Phi302 + g23*Phi223*Phi302 + g13*Phi102*Phi323 + g23*Phi202*Phi323 + g33*Phi302*Phi323 - 
            Pi02*Pi23) + invPsi33*(-(Gamma003*Gamma203*invPsi00) - Gamma013*Gamma203*invPsi01 - Gamma003*Gamma213*invPsi01 - 
            Gamma023*Gamma203*invPsi02 - Gamma003*Gamma223*invPsi02 - Gamma013*Gamma213*invPsi11 - Gamma023*Gamma213*invPsi12 - 
            Gamma013*Gamma223*invPsi12 - Gamma023*Gamma223*invPsi22 + g11*Phi103*Phi123 + g12*Phi123*Phi203 + g12*Phi103*Phi223 + 
            g22*Phi203*Phi223 + g13*Phi123*Phi303 + g23*Phi223*Phi303 + g13*Phi103*Phi323 + g23*Phi203*Phi323 + g33*Phi303*Phi323 - 
            Pi03*Pi23));
    
    t1srcPi03 = 2*N*(invPsi00*(-(Gamma000*Gamma300*invPsi00) - Gamma001*Gamma300*invPsi01 - Gamma000*Gamma301*invPsi01 - 
            Gamma002*Gamma300*invPsi02 - Gamma000*Gamma302*invPsi02 - Gamma001*Gamma301*invPsi11 - Gamma002*Gamma301*invPsi12 - 
            Gamma001*Gamma302*invPsi12 - Gamma002*Gamma302*invPsi22 + g11*Phi100*Phi103 + g12*Phi103*Phi200 + g12*Phi100*Phi203 + 
            g22*Phi200*Phi203 + g13*Phi103*Phi300 + g23*Phi203*Phi300 + g13*Phi100*Phi303 + g23*Phi200*Phi303 + g33*Phi300*Phi303 - 
            Pi00*Pi03) + invPsi01*(-(Gamma001*Gamma300*invPsi00) - Gamma011*Gamma300*invPsi01 - Gamma001*Gamma301*invPsi01 - 
            Gamma012*Gamma300*invPsi02 - Gamma001*Gamma302*invPsi02 - Gamma011*Gamma301*invPsi11 - Gamma012*Gamma301*invPsi12 - 
            Gamma011*Gamma302*invPsi12 - Gamma012*Gamma302*invPsi22 + g11*Phi101*Phi103 + g12*Phi103*Phi201 + g12*Phi101*Phi203 + 
            g22*Phi201*Phi203 + g13*Phi103*Phi301 + g23*Phi203*Phi301 + g13*Phi101*Phi303 + g23*Phi201*Phi303 + g33*Phi301*Phi303 - 
            Pi01*Pi03) + invPsi02*(-(Gamma002*Gamma300*invPsi00) - Gamma012*Gamma300*invPsi01 - Gamma002*Gamma301*invPsi01 - 
            Gamma022*Gamma300*invPsi02 - Gamma002*Gamma302*invPsi02 - Gamma012*Gamma301*invPsi11 - Gamma022*Gamma301*invPsi12 - 
            Gamma012*Gamma302*invPsi12 - Gamma022*Gamma302*invPsi22 + g11*Phi102*Phi103 + g12*Phi103*Phi202 + g12*Phi102*Phi203 + 
            g22*Phi202*Phi203 + g13*Phi103*Phi302 + g23*Phi203*Phi302 + g13*Phi102*Phi303 + g23*Phi202*Phi303 + g33*Phi302*Phi303 - 
            Pi02*Pi03) + invPsi03*(-(Gamma003*Gamma300*invPsi00) - Gamma013*Gamma300*invPsi01 - Gamma003*Gamma301*invPsi01 - 
            Gamma023*Gamma300*invPsi02 - Gamma003*Gamma302*invPsi02 - Gamma013*Gamma301*invPsi11 - Gamma023*Gamma301*invPsi12 - 
            Gamma013*Gamma302*invPsi12 - Gamma023*Gamma302*invPsi22 + g11*Power(Phi103,2) + 2*g12*Phi103*Phi203 + 
            g22*Power(Phi203,2) + 2*g13*Phi103*Phi303 + 2*g23*Phi203*Phi303 + g33*Power(Phi303,2) - Power(Pi03,2)) + 
         invPsi01*(-(Gamma000*Gamma301*invPsi00) - Gamma001*Gamma301*invPsi01 - Gamma000*Gamma311*invPsi01 - 
            Gamma002*Gamma301*invPsi02 - Gamma000*Gamma312*invPsi02 - Gamma001*Gamma311*invPsi11 - Gamma002*Gamma311*invPsi12 - 
            Gamma001*Gamma312*invPsi12 - Gamma002*Gamma312*invPsi22 + g11*Phi100*Phi113 + g12*Phi113*Phi200 + g12*Phi100*Phi213 + 
            g22*Phi200*Phi213 + g13*Phi113*Phi300 + g23*Phi213*Phi300 + g13*Phi100*Phi313 + g23*Phi200*Phi313 + g33*Phi300*Phi313 - 
            Pi00*Pi13) + invPsi11*(-(Gamma001*Gamma301*invPsi00) - Gamma011*Gamma301*invPsi01 - Gamma001*Gamma311*invPsi01 - 
            Gamma012*Gamma301*invPsi02 - Gamma001*Gamma312*invPsi02 - Gamma011*Gamma311*invPsi11 - Gamma012*Gamma311*invPsi12 - 
            Gamma011*Gamma312*invPsi12 - Gamma012*Gamma312*invPsi22 + g11*Phi101*Phi113 + g12*Phi113*Phi201 + g12*Phi101*Phi213 + 
            g22*Phi201*Phi213 + g13*Phi113*Phi301 + g23*Phi213*Phi301 + g13*Phi101*Phi313 + g23*Phi201*Phi313 + g33*Phi301*Phi313 - 
            Pi01*Pi13) + invPsi12*(-(Gamma002*Gamma301*invPsi00) - Gamma012*Gamma301*invPsi01 - Gamma002*Gamma311*invPsi01 - 
            Gamma022*Gamma301*invPsi02 - Gamma002*Gamma312*invPsi02 - Gamma012*Gamma311*invPsi11 - Gamma022*Gamma311*invPsi12 - 
            Gamma012*Gamma312*invPsi12 - Gamma022*Gamma312*invPsi22 + g11*Phi102*Phi113 + g12*Phi113*Phi202 + g12*Phi102*Phi213 + 
            g22*Phi202*Phi213 + g13*Phi113*Phi302 + g23*Phi213*Phi302 + g13*Phi102*Phi313 + g23*Phi202*Phi313 + g33*Phi302*Phi313 - 
            Pi02*Pi13) + invPsi13*(-(Gamma003*Gamma301*invPsi00) - Gamma013*Gamma301*invPsi01 - Gamma003*Gamma311*invPsi01 - 
            Gamma023*Gamma301*invPsi02 - Gamma003*Gamma312*invPsi02 - Gamma013*Gamma311*invPsi11 - Gamma023*Gamma311*invPsi12 - 
            Gamma013*Gamma312*invPsi12 - Gamma023*Gamma312*invPsi22 + g11*Phi103*Phi113 + g12*Phi113*Phi203 + g12*Phi103*Phi213 + 
            g22*Phi203*Phi213 + g13*Phi113*Phi303 + g23*Phi213*Phi303 + g13*Phi103*Phi313 + g23*Phi203*Phi313 + g33*Phi303*Phi313 - 
            Pi03*Pi13) + invPsi02*(-(Gamma000*Gamma302*invPsi00) - Gamma001*Gamma302*invPsi01 - Gamma000*Gamma312*invPsi01 - 
            Gamma002*Gamma302*invPsi02 - Gamma000*Gamma322*invPsi02 - Gamma001*Gamma312*invPsi11 - Gamma002*Gamma312*invPsi12 - 
            Gamma001*Gamma322*invPsi12 - Gamma002*Gamma322*invPsi22 + g11*Phi100*Phi123 + g12*Phi123*Phi200 + g12*Phi100*Phi223 + 
            g22*Phi200*Phi223 + g13*Phi123*Phi300 + g23*Phi223*Phi300 + g13*Phi100*Phi323 + g23*Phi200*Phi323 + g33*Phi300*Phi323 - 
            Pi00*Pi23) + invPsi12*(-(Gamma001*Gamma302*invPsi00) - Gamma011*Gamma302*invPsi01 - Gamma001*Gamma312*invPsi01 - 
            Gamma012*Gamma302*invPsi02 - Gamma001*Gamma322*invPsi02 - Gamma011*Gamma312*invPsi11 - Gamma012*Gamma312*invPsi12 - 
            Gamma011*Gamma322*invPsi12 - Gamma012*Gamma322*invPsi22 + g11*Phi101*Phi123 + g12*Phi123*Phi201 + g12*Phi101*Phi223 + 
            g22*Phi201*Phi223 + g13*Phi123*Phi301 + g23*Phi223*Phi301 + g13*Phi101*Phi323 + g23*Phi201*Phi323 + g33*Phi301*Phi323 - 
            Pi01*Pi23) + invPsi22*(-(Gamma002*Gamma302*invPsi00) - Gamma012*Gamma302*invPsi01 - Gamma002*Gamma312*invPsi01 - 
            Gamma022*Gamma302*invPsi02 - Gamma002*Gamma322*invPsi02 - Gamma012*Gamma312*invPsi11 - Gamma022*Gamma312*invPsi12 - 
            Gamma012*Gamma322*invPsi12 - Gamma022*Gamma322*invPsi22 + g11*Phi102*Phi123 + g12*Phi123*Phi202 + g12*Phi102*Phi223 + 
            g22*Phi202*Phi223 + g13*Phi123*Phi302 + g23*Phi223*Phi302 + g13*Phi102*Phi323 + g23*Phi202*Phi323 + g33*Phi302*Phi323 - 
            Pi02*Pi23) + invPsi23*(-(Gamma003*Gamma302*invPsi00) - Gamma013*Gamma302*invPsi01 - Gamma003*Gamma312*invPsi01 - 
            Gamma023*Gamma302*invPsi02 - Gamma003*Gamma322*invPsi02 - Gamma013*Gamma312*invPsi11 - Gamma023*Gamma312*invPsi12 - 
            Gamma013*Gamma322*invPsi12 - Gamma023*Gamma322*invPsi22 + g11*Phi103*Phi123 + g12*Phi123*Phi203 + g12*Phi103*Phi223 + 
            g22*Phi203*Phi223 + g13*Phi123*Phi303 + g23*Phi223*Phi303 + g13*Phi103*Phi323 + g23*Phi203*Phi323 + g33*Phi303*Phi323 - 
            Pi03*Pi23) + invPsi03*(-(Gamma000*Gamma303*invPsi00) - Gamma001*Gamma303*invPsi01 - Gamma000*Gamma313*invPsi01 - 
            Gamma002*Gamma303*invPsi02 - Gamma000*Gamma323*invPsi02 - Gamma001*Gamma313*invPsi11 - Gamma002*Gamma313*invPsi12 - 
            Gamma001*Gamma323*invPsi12 - Gamma002*Gamma323*invPsi22 + g11*Phi100*Phi133 + g12*Phi133*Phi200 + g12*Phi100*Phi233 + 
            g22*Phi200*Phi233 + g13*Phi133*Phi300 + g23*Phi233*Phi300 + g13*Phi100*Phi333 + g23*Phi200*Phi333 + g33*Phi300*Phi333 - 
            Pi00*Pi33) + invPsi13*(-(Gamma001*Gamma303*invPsi00) - Gamma011*Gamma303*invPsi01 - Gamma001*Gamma313*invPsi01 - 
            Gamma012*Gamma303*invPsi02 - Gamma001*Gamma323*invPsi02 - Gamma011*Gamma313*invPsi11 - Gamma012*Gamma313*invPsi12 - 
            Gamma011*Gamma323*invPsi12 - Gamma012*Gamma323*invPsi22 + g11*Phi101*Phi133 + g12*Phi133*Phi201 + g12*Phi101*Phi233 + 
            g22*Phi201*Phi233 + g13*Phi133*Phi301 + g23*Phi233*Phi301 + g13*Phi101*Phi333 + g23*Phi201*Phi333 + g33*Phi301*Phi333 - 
            Pi01*Pi33) + invPsi23*(-(Gamma002*Gamma303*invPsi00) - Gamma012*Gamma303*invPsi01 - Gamma002*Gamma313*invPsi01 - 
            Gamma022*Gamma303*invPsi02 - Gamma002*Gamma323*invPsi02 - Gamma012*Gamma313*invPsi11 - Gamma022*Gamma313*invPsi12 - 
            Gamma012*Gamma323*invPsi12 - Gamma022*Gamma323*invPsi22 + g11*Phi102*Phi133 + g12*Phi133*Phi202 + g12*Phi102*Phi233 + 
            g22*Phi202*Phi233 + g13*Phi133*Phi302 + g23*Phi233*Phi302 + g13*Phi102*Phi333 + g23*Phi202*Phi333 + g33*Phi302*Phi333 - 
            Pi02*Pi33) + invPsi33*(-(Gamma003*Gamma303*invPsi00) - Gamma013*Gamma303*invPsi01 - Gamma003*Gamma313*invPsi01 - 
            Gamma023*Gamma303*invPsi02 - Gamma003*Gamma323*invPsi02 - Gamma013*Gamma313*invPsi11 - Gamma023*Gamma313*invPsi12 - 
            Gamma013*Gamma323*invPsi12 - Gamma023*Gamma323*invPsi22 + g11*Phi103*Phi133 + g12*Phi133*Phi203 + g12*Phi103*Phi233 + 
            g22*Phi203*Phi233 + g13*Phi133*Phi303 + g23*Phi233*Phi303 + g13*Phi103*Phi333 + g23*Phi203*Phi333 + g33*Phi303*Phi333 - 
            Pi03*Pi33));
    
    t1srcPi11 = 2*N*(invPsi00*(-(Power(Gamma100,2)*invPsi00) - 2*Gamma100*Gamma101*invPsi01 - 2*Gamma100*Gamma102*invPsi02 - 
            Power(Gamma101,2)*invPsi11 - 2*Gamma101*Gamma102*invPsi12 - Power(Gamma102,2)*invPsi22 + g11*Power(Phi101,2) + 
            2*g12*Phi101*Phi201 + g22*Power(Phi201,2) + 2*g13*Phi101*Phi301 + 2*g23*Phi201*Phi301 + g33*Power(Phi301,2) - 
            Power(Pi01,2)) + 2*invPsi01*(-(Gamma100*Gamma101*invPsi00) - Power(Gamma101,2)*invPsi01 - Gamma100*Gamma111*invPsi01 - 
            Gamma101*Gamma102*invPsi02 - Gamma100*Gamma112*invPsi02 - Gamma101*Gamma111*invPsi11 - Gamma102*Gamma111*invPsi12 - 
            Gamma101*Gamma112*invPsi12 - Gamma102*Gamma112*invPsi22 + g11*Phi101*Phi111 + g12*Phi111*Phi201 + g12*Phi101*Phi211 + 
            g22*Phi201*Phi211 + g13*Phi111*Phi301 + g23*Phi211*Phi301 + g13*Phi101*Phi311 + g23*Phi201*Phi311 + g33*Phi301*Phi311 - 
            Pi01*Pi11) + invPsi11*(-(Power(Gamma101,2)*invPsi00) - 2*Gamma101*Gamma111*invPsi01 - 2*Gamma101*Gamma112*invPsi02 - 
            Power(Gamma111,2)*invPsi11 - 2*Gamma111*Gamma112*invPsi12 - Power(Gamma112,2)*invPsi22 + g11*Power(Phi111,2) + 
            2*g12*Phi111*Phi211 + g22*Power(Phi211,2) + 2*g13*Phi111*Phi311 + 2*g23*Phi211*Phi311 + g33*Power(Phi311,2) - 
            Power(Pi11,2)) + 2*invPsi02*(-(Gamma100*Gamma102*invPsi00) - Gamma101*Gamma102*invPsi01 - Gamma100*Gamma112*invPsi01 - 
            Power(Gamma102,2)*invPsi02 - Gamma100*Gamma122*invPsi02 - Gamma101*Gamma112*invPsi11 - Gamma102*Gamma112*invPsi12 - 
            Gamma101*Gamma122*invPsi12 - Gamma102*Gamma122*invPsi22 + g11*Phi101*Phi112 + g12*Phi112*Phi201 + g12*Phi101*Phi212 + 
            g22*Phi201*Phi212 + g13*Phi112*Phi301 + g23*Phi212*Phi301 + g13*Phi101*Phi312 + g23*Phi201*Phi312 + g33*Phi301*Phi312 - 
            Pi01*Pi12) + 2*invPsi12*(-(Gamma101*Gamma102*invPsi00) - Gamma102*Gamma111*invPsi01 - Gamma101*Gamma112*invPsi01 - 
            Gamma102*Gamma112*invPsi02 - Gamma101*Gamma122*invPsi02 - Gamma111*Gamma112*invPsi11 - Power(Gamma112,2)*invPsi12 - 
            Gamma111*Gamma122*invPsi12 - Gamma112*Gamma122*invPsi22 + g11*Phi111*Phi112 + g12*Phi112*Phi211 + g12*Phi111*Phi212 + 
            g22*Phi211*Phi212 + g13*Phi112*Phi311 + g23*Phi212*Phi311 + g13*Phi111*Phi312 + g23*Phi211*Phi312 + g33*Phi311*Phi312 - 
            Pi11*Pi12) + invPsi22*(-(Power(Gamma102,2)*invPsi00) - 2*Gamma102*Gamma112*invPsi01 - 2*Gamma102*Gamma122*invPsi02 - 
            Power(Gamma112,2)*invPsi11 - 2*Gamma112*Gamma122*invPsi12 - Power(Gamma122,2)*invPsi22 + g11*Power(Phi112,2) + 
            2*g12*Phi112*Phi212 + g22*Power(Phi212,2) + 2*g13*Phi112*Phi312 + 2*g23*Phi212*Phi312 + g33*Power(Phi312,2) - 
            Power(Pi12,2)) + 2*invPsi03*(-(Gamma100*Gamma103*invPsi00) - Gamma101*Gamma103*invPsi01 - Gamma100*Gamma113*invPsi01 - 
            Gamma102*Gamma103*invPsi02 - Gamma100*Gamma123*invPsi02 - Gamma101*Gamma113*invPsi11 - Gamma102*Gamma113*invPsi12 - 
            Gamma101*Gamma123*invPsi12 - Gamma102*Gamma123*invPsi22 + g11*Phi101*Phi113 + g12*Phi113*Phi201 + g12*Phi101*Phi213 + 
            g22*Phi201*Phi213 + g13*Phi113*Phi301 + g23*Phi213*Phi301 + g13*Phi101*Phi313 + g23*Phi201*Phi313 + g33*Phi301*Phi313 - 
            Pi01*Pi13) + 2*invPsi13*(-(Gamma101*Gamma103*invPsi00) - Gamma103*Gamma111*invPsi01 - Gamma101*Gamma113*invPsi01 - 
            Gamma103*Gamma112*invPsi02 - Gamma101*Gamma123*invPsi02 - Gamma111*Gamma113*invPsi11 - Gamma112*Gamma113*invPsi12 - 
            Gamma111*Gamma123*invPsi12 - Gamma112*Gamma123*invPsi22 + g11*Phi111*Phi113 + g12*Phi113*Phi211 + g12*Phi111*Phi213 + 
            g22*Phi211*Phi213 + g13*Phi113*Phi311 + g23*Phi213*Phi311 + g13*Phi111*Phi313 + g23*Phi211*Phi313 + g33*Phi311*Phi313 - 
            Pi11*Pi13) + 2*invPsi23*(-(Gamma102*Gamma103*invPsi00) - Gamma103*Gamma112*invPsi01 - Gamma102*Gamma113*invPsi01 - 
            Gamma103*Gamma122*invPsi02 - Gamma102*Gamma123*invPsi02 - Gamma112*Gamma113*invPsi11 - Gamma113*Gamma122*invPsi12 - 
            Gamma112*Gamma123*invPsi12 - Gamma122*Gamma123*invPsi22 + g11*Phi112*Phi113 + g12*Phi113*Phi212 + g12*Phi112*Phi213 + 
            g22*Phi212*Phi213 + g13*Phi113*Phi312 + g23*Phi213*Phi312 + g13*Phi112*Phi313 + g23*Phi212*Phi313 + g33*Phi312*Phi313 - 
            Pi12*Pi13) + invPsi33*(-(Power(Gamma103,2)*invPsi00) - 2*Gamma103*Gamma113*invPsi01 - 2*Gamma103*Gamma123*invPsi02 - 
            Power(Gamma113,2)*invPsi11 - 2*Gamma113*Gamma123*invPsi12 - Power(Gamma123,2)*invPsi22 + g11*Power(Phi113,2) + 
            2*g12*Phi113*Phi213 + g22*Power(Phi213,2) + 2*g13*Phi113*Phi313 + 2*g23*Phi213*Phi313 + g33*Power(Phi313,2) - 
            Power(Pi13,2)));
    
    t1srcPi12 = 2*N*(invPsi00*(-(Gamma100*Gamma200*invPsi00) - Gamma101*Gamma200*invPsi01 - Gamma100*Gamma201*invPsi01 - 
            Gamma102*Gamma200*invPsi02 - Gamma100*Gamma202*invPsi02 - Gamma101*Gamma201*invPsi11 - Gamma102*Gamma201*invPsi12 - 
            Gamma101*Gamma202*invPsi12 - Gamma102*Gamma202*invPsi22 + g11*Phi101*Phi102 + g12*Phi102*Phi201 + g12*Phi101*Phi202 + 
            g22*Phi201*Phi202 + g13*Phi102*Phi301 + g23*Phi202*Phi301 + g13*Phi101*Phi302 + g23*Phi201*Phi302 + g33*Phi301*Phi302 - 
            Pi01*Pi02) + invPsi01*(-(Gamma101*Gamma200*invPsi00) - Gamma111*Gamma200*invPsi01 - Gamma101*Gamma201*invPsi01 - 
            Gamma112*Gamma200*invPsi02 - Gamma101*Gamma202*invPsi02 - Gamma111*Gamma201*invPsi11 - Gamma112*Gamma201*invPsi12 - 
            Gamma111*Gamma202*invPsi12 - Gamma112*Gamma202*invPsi22 + g11*Phi102*Phi111 + g12*Phi111*Phi202 + g12*Phi102*Phi211 + 
            g22*Phi202*Phi211 + g13*Phi111*Phi302 + g23*Phi211*Phi302 + g13*Phi102*Phi311 + g23*Phi202*Phi311 + g33*Phi302*Phi311 - 
            Pi02*Pi11) + invPsi01*(-(Gamma100*Gamma201*invPsi00) - Gamma101*Gamma201*invPsi01 - Gamma100*Gamma211*invPsi01 - 
            Gamma102*Gamma201*invPsi02 - Gamma100*Gamma212*invPsi02 - Gamma101*Gamma211*invPsi11 - Gamma102*Gamma211*invPsi12 - 
            Gamma101*Gamma212*invPsi12 - Gamma102*Gamma212*invPsi22 + g11*Phi101*Phi112 + g12*Phi112*Phi201 + g12*Phi101*Phi212 + 
            g22*Phi201*Phi212 + g13*Phi112*Phi301 + g23*Phi212*Phi301 + g13*Phi101*Phi312 + g23*Phi201*Phi312 + g33*Phi301*Phi312 - 
            Pi01*Pi12) + invPsi02*(-(Gamma102*Gamma200*invPsi00) - Gamma112*Gamma200*invPsi01 - Gamma102*Gamma201*invPsi01 - 
            Gamma122*Gamma200*invPsi02 - Gamma102*Gamma202*invPsi02 - Gamma112*Gamma201*invPsi11 - Gamma122*Gamma201*invPsi12 - 
            Gamma112*Gamma202*invPsi12 - Gamma122*Gamma202*invPsi22 + g11*Phi102*Phi112 + g12*Phi112*Phi202 + g12*Phi102*Phi212 + 
            g22*Phi202*Phi212 + g13*Phi112*Phi302 + g23*Phi212*Phi302 + g13*Phi102*Phi312 + g23*Phi202*Phi312 + g33*Phi302*Phi312 - 
            Pi02*Pi12) + invPsi11*(-(Gamma101*Gamma201*invPsi00) - Gamma111*Gamma201*invPsi01 - Gamma101*Gamma211*invPsi01 - 
            Gamma112*Gamma201*invPsi02 - Gamma101*Gamma212*invPsi02 - Gamma111*Gamma211*invPsi11 - Gamma112*Gamma211*invPsi12 - 
            Gamma111*Gamma212*invPsi12 - Gamma112*Gamma212*invPsi22 + g11*Phi111*Phi112 + g12*Phi112*Phi211 + g12*Phi111*Phi212 + 
            g22*Phi211*Phi212 + g13*Phi112*Phi311 + g23*Phi212*Phi311 + g13*Phi111*Phi312 + g23*Phi211*Phi312 + g33*Phi311*Phi312 - 
            Pi11*Pi12) + invPsi12*(-(Gamma102*Gamma201*invPsi00) - Gamma112*Gamma201*invPsi01 - Gamma102*Gamma211*invPsi01 - 
            Gamma122*Gamma201*invPsi02 - Gamma102*Gamma212*invPsi02 - Gamma112*Gamma211*invPsi11 - Gamma122*Gamma211*invPsi12 - 
            Gamma112*Gamma212*invPsi12 - Gamma122*Gamma212*invPsi22 + g11*Power(Phi112,2) + 2*g12*Phi112*Phi212 + 
            g22*Power(Phi212,2) + 2*g13*Phi112*Phi312 + 2*g23*Phi212*Phi312 + g33*Power(Phi312,2) - Power(Pi12,2)) + 
         invPsi03*(-(Gamma103*Gamma200*invPsi00) - Gamma113*Gamma200*invPsi01 - Gamma103*Gamma201*invPsi01 - 
            Gamma123*Gamma200*invPsi02 - Gamma103*Gamma202*invPsi02 - Gamma113*Gamma201*invPsi11 - Gamma123*Gamma201*invPsi12 - 
            Gamma113*Gamma202*invPsi12 - Gamma123*Gamma202*invPsi22 + g11*Phi102*Phi113 + g12*Phi113*Phi202 + g12*Phi102*Phi213 + 
            g22*Phi202*Phi213 + g13*Phi113*Phi302 + g23*Phi213*Phi302 + g13*Phi102*Phi313 + g23*Phi202*Phi313 + g33*Phi302*Phi313 - 
            Pi02*Pi13) + invPsi13*(-(Gamma103*Gamma201*invPsi00) - Gamma113*Gamma201*invPsi01 - Gamma103*Gamma211*invPsi01 - 
            Gamma123*Gamma201*invPsi02 - Gamma103*Gamma212*invPsi02 - Gamma113*Gamma211*invPsi11 - Gamma123*Gamma211*invPsi12 - 
            Gamma113*Gamma212*invPsi12 - Gamma123*Gamma212*invPsi22 + g11*Phi112*Phi113 + g12*Phi113*Phi212 + g12*Phi112*Phi213 + 
            g22*Phi212*Phi213 + g13*Phi113*Phi312 + g23*Phi213*Phi312 + g13*Phi112*Phi313 + g23*Phi212*Phi313 + g33*Phi312*Phi313 - 
            Pi12*Pi13) + invPsi02*(-(Gamma100*Gamma202*invPsi00) - Gamma101*Gamma202*invPsi01 - Gamma100*Gamma212*invPsi01 - 
            Gamma102*Gamma202*invPsi02 - Gamma100*Gamma222*invPsi02 - Gamma101*Gamma212*invPsi11 - Gamma102*Gamma212*invPsi12 - 
            Gamma101*Gamma222*invPsi12 - Gamma102*Gamma222*invPsi22 + g11*Phi101*Phi122 + g12*Phi122*Phi201 + g12*Phi101*Phi222 + 
            g22*Phi201*Phi222 + g13*Phi122*Phi301 + g23*Phi222*Phi301 + g13*Phi101*Phi322 + g23*Phi201*Phi322 + g33*Phi301*Phi322 - 
            Pi01*Pi22) + invPsi12*(-(Gamma101*Gamma202*invPsi00) - Gamma111*Gamma202*invPsi01 - Gamma101*Gamma212*invPsi01 - 
            Gamma112*Gamma202*invPsi02 - Gamma101*Gamma222*invPsi02 - Gamma111*Gamma212*invPsi11 - Gamma112*Gamma212*invPsi12 - 
            Gamma111*Gamma222*invPsi12 - Gamma112*Gamma222*invPsi22 + g11*Phi111*Phi122 + g12*Phi122*Phi211 + g12*Phi111*Phi222 + 
            g22*Phi211*Phi222 + g13*Phi122*Phi311 + g23*Phi222*Phi311 + g13*Phi111*Phi322 + g23*Phi211*Phi322 + g33*Phi311*Phi322 - 
            Pi11*Pi22) + invPsi22*(-(Gamma102*Gamma202*invPsi00) - Gamma112*Gamma202*invPsi01 - Gamma102*Gamma212*invPsi01 - 
            Gamma122*Gamma202*invPsi02 - Gamma102*Gamma222*invPsi02 - Gamma112*Gamma212*invPsi11 - Gamma122*Gamma212*invPsi12 - 
            Gamma112*Gamma222*invPsi12 - Gamma122*Gamma222*invPsi22 + g11*Phi112*Phi122 + g12*Phi122*Phi212 + g12*Phi112*Phi222 + 
            g22*Phi212*Phi222 + g13*Phi122*Phi312 + g23*Phi222*Phi312 + g13*Phi112*Phi322 + g23*Phi212*Phi322 + g33*Phi312*Phi322 - 
            Pi12*Pi22) + invPsi23*(-(Gamma103*Gamma202*invPsi00) - Gamma113*Gamma202*invPsi01 - Gamma103*Gamma212*invPsi01 - 
            Gamma123*Gamma202*invPsi02 - Gamma103*Gamma222*invPsi02 - Gamma113*Gamma212*invPsi11 - Gamma123*Gamma212*invPsi12 - 
            Gamma113*Gamma222*invPsi12 - Gamma123*Gamma222*invPsi22 + g11*Phi113*Phi122 + g12*Phi122*Phi213 + g12*Phi113*Phi222 + 
            g22*Phi213*Phi222 + g13*Phi122*Phi313 + g23*Phi222*Phi313 + g13*Phi113*Phi322 + g23*Phi213*Phi322 + g33*Phi313*Phi322 - 
            Pi13*Pi22) + invPsi03*(-(Gamma100*Gamma203*invPsi00) - Gamma101*Gamma203*invPsi01 - Gamma100*Gamma213*invPsi01 - 
            Gamma102*Gamma203*invPsi02 - Gamma100*Gamma223*invPsi02 - Gamma101*Gamma213*invPsi11 - Gamma102*Gamma213*invPsi12 - 
            Gamma101*Gamma223*invPsi12 - Gamma102*Gamma223*invPsi22 + g11*Phi101*Phi123 + g12*Phi123*Phi201 + g12*Phi101*Phi223 + 
            g22*Phi201*Phi223 + g13*Phi123*Phi301 + g23*Phi223*Phi301 + g13*Phi101*Phi323 + g23*Phi201*Phi323 + g33*Phi301*Phi323 - 
            Pi01*Pi23) + invPsi13*(-(Gamma101*Gamma203*invPsi00) - Gamma111*Gamma203*invPsi01 - Gamma101*Gamma213*invPsi01 - 
            Gamma112*Gamma203*invPsi02 - Gamma101*Gamma223*invPsi02 - Gamma111*Gamma213*invPsi11 - Gamma112*Gamma213*invPsi12 - 
            Gamma111*Gamma223*invPsi12 - Gamma112*Gamma223*invPsi22 + g11*Phi111*Phi123 + g12*Phi123*Phi211 + g12*Phi111*Phi223 + 
            g22*Phi211*Phi223 + g13*Phi123*Phi311 + g23*Phi223*Phi311 + g13*Phi111*Phi323 + g23*Phi211*Phi323 + g33*Phi311*Phi323 - 
            Pi11*Pi23) + invPsi23*(-(Gamma102*Gamma203*invPsi00) - Gamma112*Gamma203*invPsi01 - Gamma102*Gamma213*invPsi01 - 
            Gamma122*Gamma203*invPsi02 - Gamma102*Gamma223*invPsi02 - Gamma112*Gamma213*invPsi11 - Gamma122*Gamma213*invPsi12 - 
            Gamma112*Gamma223*invPsi12 - Gamma122*Gamma223*invPsi22 + g11*Phi112*Phi123 + g12*Phi123*Phi212 + g12*Phi112*Phi223 + 
            g22*Phi212*Phi223 + g13*Phi123*Phi312 + g23*Phi223*Phi312 + g13*Phi112*Phi323 + g23*Phi212*Phi323 + g33*Phi312*Phi323 - 
            Pi12*Pi23) + invPsi33*(-(Gamma103*Gamma203*invPsi00) - Gamma113*Gamma203*invPsi01 - Gamma103*Gamma213*invPsi01 - 
            Gamma123*Gamma203*invPsi02 - Gamma103*Gamma223*invPsi02 - Gamma113*Gamma213*invPsi11 - Gamma123*Gamma213*invPsi12 - 
            Gamma113*Gamma223*invPsi12 - Gamma123*Gamma223*invPsi22 + g11*Phi113*Phi123 + g12*Phi123*Phi213 + g12*Phi113*Phi223 + 
            g22*Phi213*Phi223 + g13*Phi123*Phi313 + g23*Phi223*Phi313 + g13*Phi113*Phi323 + g23*Phi213*Phi323 + g33*Phi313*Phi323 - 
            Pi13*Pi23));
    
    t1srcPi13 = 2*N*(invPsi00*(-(Gamma100*Gamma300*invPsi00) - Gamma101*Gamma300*invPsi01 - Gamma100*Gamma301*invPsi01 - 
            Gamma102*Gamma300*invPsi02 - Gamma100*Gamma302*invPsi02 - Gamma101*Gamma301*invPsi11 - Gamma102*Gamma301*invPsi12 - 
            Gamma101*Gamma302*invPsi12 - Gamma102*Gamma302*invPsi22 + g11*Phi101*Phi103 + g12*Phi103*Phi201 + g12*Phi101*Phi203 + 
            g22*Phi201*Phi203 + g13*Phi103*Phi301 + g23*Phi203*Phi301 + g13*Phi101*Phi303 + g23*Phi201*Phi303 + g33*Phi301*Phi303 - 
            Pi01*Pi03) + invPsi01*(-(Gamma101*Gamma300*invPsi00) - Gamma111*Gamma300*invPsi01 - Gamma101*Gamma301*invPsi01 - 
            Gamma112*Gamma300*invPsi02 - Gamma101*Gamma302*invPsi02 - Gamma111*Gamma301*invPsi11 - Gamma112*Gamma301*invPsi12 - 
            Gamma111*Gamma302*invPsi12 - Gamma112*Gamma302*invPsi22 + g11*Phi103*Phi111 + g12*Phi111*Phi203 + g12*Phi103*Phi211 + 
            g22*Phi203*Phi211 + g13*Phi111*Phi303 + g23*Phi211*Phi303 + g13*Phi103*Phi311 + g23*Phi203*Phi311 + g33*Phi303*Phi311 - 
            Pi03*Pi11) + invPsi02*(-(Gamma102*Gamma300*invPsi00) - Gamma112*Gamma300*invPsi01 - Gamma102*Gamma301*invPsi01 - 
            Gamma122*Gamma300*invPsi02 - Gamma102*Gamma302*invPsi02 - Gamma112*Gamma301*invPsi11 - Gamma122*Gamma301*invPsi12 - 
            Gamma112*Gamma302*invPsi12 - Gamma122*Gamma302*invPsi22 + g11*Phi103*Phi112 + g12*Phi112*Phi203 + g12*Phi103*Phi212 + 
            g22*Phi203*Phi212 + g13*Phi112*Phi303 + g23*Phi212*Phi303 + g13*Phi103*Phi312 + g23*Phi203*Phi312 + g33*Phi303*Phi312 - 
            Pi03*Pi12) + invPsi01*(-(Gamma100*Gamma301*invPsi00) - Gamma101*Gamma301*invPsi01 - Gamma100*Gamma311*invPsi01 - 
            Gamma102*Gamma301*invPsi02 - Gamma100*Gamma312*invPsi02 - Gamma101*Gamma311*invPsi11 - Gamma102*Gamma311*invPsi12 - 
            Gamma101*Gamma312*invPsi12 - Gamma102*Gamma312*invPsi22 + g11*Phi101*Phi113 + g12*Phi113*Phi201 + g12*Phi101*Phi213 + 
            g22*Phi201*Phi213 + g13*Phi113*Phi301 + g23*Phi213*Phi301 + g13*Phi101*Phi313 + g23*Phi201*Phi313 + g33*Phi301*Phi313 - 
            Pi01*Pi13) + invPsi03*(-(Gamma103*Gamma300*invPsi00) - Gamma113*Gamma300*invPsi01 - Gamma103*Gamma301*invPsi01 - 
            Gamma123*Gamma300*invPsi02 - Gamma103*Gamma302*invPsi02 - Gamma113*Gamma301*invPsi11 - Gamma123*Gamma301*invPsi12 - 
            Gamma113*Gamma302*invPsi12 - Gamma123*Gamma302*invPsi22 + g11*Phi103*Phi113 + g12*Phi113*Phi203 + g12*Phi103*Phi213 + 
            g22*Phi203*Phi213 + g13*Phi113*Phi303 + g23*Phi213*Phi303 + g13*Phi103*Phi313 + g23*Phi203*Phi313 + g33*Phi303*Phi313 - 
            Pi03*Pi13) + invPsi11*(-(Gamma101*Gamma301*invPsi00) - Gamma111*Gamma301*invPsi01 - Gamma101*Gamma311*invPsi01 - 
            Gamma112*Gamma301*invPsi02 - Gamma101*Gamma312*invPsi02 - Gamma111*Gamma311*invPsi11 - Gamma112*Gamma311*invPsi12 - 
            Gamma111*Gamma312*invPsi12 - Gamma112*Gamma312*invPsi22 + g11*Phi111*Phi113 + g12*Phi113*Phi211 + g12*Phi111*Phi213 + 
            g22*Phi211*Phi213 + g13*Phi113*Phi311 + g23*Phi213*Phi311 + g13*Phi111*Phi313 + g23*Phi211*Phi313 + g33*Phi311*Phi313 - 
            Pi11*Pi13) + invPsi12*(-(Gamma102*Gamma301*invPsi00) - Gamma112*Gamma301*invPsi01 - Gamma102*Gamma311*invPsi01 - 
            Gamma122*Gamma301*invPsi02 - Gamma102*Gamma312*invPsi02 - Gamma112*Gamma311*invPsi11 - Gamma122*Gamma311*invPsi12 - 
            Gamma112*Gamma312*invPsi12 - Gamma122*Gamma312*invPsi22 + g11*Phi112*Phi113 + g12*Phi113*Phi212 + g12*Phi112*Phi213 + 
            g22*Phi212*Phi213 + g13*Phi113*Phi312 + g23*Phi213*Phi312 + g13*Phi112*Phi313 + g23*Phi212*Phi313 + g33*Phi312*Phi313 - 
            Pi12*Pi13) + invPsi13*(-(Gamma103*Gamma301*invPsi00) - Gamma113*Gamma301*invPsi01 - Gamma103*Gamma311*invPsi01 - 
            Gamma123*Gamma301*invPsi02 - Gamma103*Gamma312*invPsi02 - Gamma113*Gamma311*invPsi11 - Gamma123*Gamma311*invPsi12 - 
            Gamma113*Gamma312*invPsi12 - Gamma123*Gamma312*invPsi22 + g11*Power(Phi113,2) + 2*g12*Phi113*Phi213 + 
            g22*Power(Phi213,2) + 2*g13*Phi113*Phi313 + 2*g23*Phi213*Phi313 + g33*Power(Phi313,2) - Power(Pi13,2)) + 
         invPsi02*(-(Gamma100*Gamma302*invPsi00) - Gamma101*Gamma302*invPsi01 - Gamma100*Gamma312*invPsi01 - 
            Gamma102*Gamma302*invPsi02 - Gamma100*Gamma322*invPsi02 - Gamma101*Gamma312*invPsi11 - Gamma102*Gamma312*invPsi12 - 
            Gamma101*Gamma322*invPsi12 - Gamma102*Gamma322*invPsi22 + g11*Phi101*Phi123 + g12*Phi123*Phi201 + g12*Phi101*Phi223 + 
            g22*Phi201*Phi223 + g13*Phi123*Phi301 + g23*Phi223*Phi301 + g13*Phi101*Phi323 + g23*Phi201*Phi323 + g33*Phi301*Phi323 - 
            Pi01*Pi23) + invPsi12*(-(Gamma101*Gamma302*invPsi00) - Gamma111*Gamma302*invPsi01 - Gamma101*Gamma312*invPsi01 - 
            Gamma112*Gamma302*invPsi02 - Gamma101*Gamma322*invPsi02 - Gamma111*Gamma312*invPsi11 - Gamma112*Gamma312*invPsi12 - 
            Gamma111*Gamma322*invPsi12 - Gamma112*Gamma322*invPsi22 + g11*Phi111*Phi123 + g12*Phi123*Phi211 + g12*Phi111*Phi223 + 
            g22*Phi211*Phi223 + g13*Phi123*Phi311 + g23*Phi223*Phi311 + g13*Phi111*Phi323 + g23*Phi211*Phi323 + g33*Phi311*Phi323 - 
            Pi11*Pi23) + invPsi22*(-(Gamma102*Gamma302*invPsi00) - Gamma112*Gamma302*invPsi01 - Gamma102*Gamma312*invPsi01 - 
            Gamma122*Gamma302*invPsi02 - Gamma102*Gamma322*invPsi02 - Gamma112*Gamma312*invPsi11 - Gamma122*Gamma312*invPsi12 - 
            Gamma112*Gamma322*invPsi12 - Gamma122*Gamma322*invPsi22 + g11*Phi112*Phi123 + g12*Phi123*Phi212 + g12*Phi112*Phi223 + 
            g22*Phi212*Phi223 + g13*Phi123*Phi312 + g23*Phi223*Phi312 + g13*Phi112*Phi323 + g23*Phi212*Phi323 + g33*Phi312*Phi323 - 
            Pi12*Pi23) + invPsi23*(-(Gamma103*Gamma302*invPsi00) - Gamma113*Gamma302*invPsi01 - Gamma103*Gamma312*invPsi01 - 
            Gamma123*Gamma302*invPsi02 - Gamma103*Gamma322*invPsi02 - Gamma113*Gamma312*invPsi11 - Gamma123*Gamma312*invPsi12 - 
            Gamma113*Gamma322*invPsi12 - Gamma123*Gamma322*invPsi22 + g11*Phi113*Phi123 + g12*Phi123*Phi213 + g12*Phi113*Phi223 + 
            g22*Phi213*Phi223 + g13*Phi123*Phi313 + g23*Phi223*Phi313 + g13*Phi113*Phi323 + g23*Phi213*Phi323 + g33*Phi313*Phi323 - 
            Pi13*Pi23) + invPsi03*(-(Gamma100*Gamma303*invPsi00) - Gamma101*Gamma303*invPsi01 - Gamma100*Gamma313*invPsi01 - 
            Gamma102*Gamma303*invPsi02 - Gamma100*Gamma323*invPsi02 - Gamma101*Gamma313*invPsi11 - Gamma102*Gamma313*invPsi12 - 
            Gamma101*Gamma323*invPsi12 - Gamma102*Gamma323*invPsi22 + g11*Phi101*Phi133 + g12*Phi133*Phi201 + g12*Phi101*Phi233 + 
            g22*Phi201*Phi233 + g13*Phi133*Phi301 + g23*Phi233*Phi301 + g13*Phi101*Phi333 + g23*Phi201*Phi333 + g33*Phi301*Phi333 - 
            Pi01*Pi33) + invPsi13*(-(Gamma101*Gamma303*invPsi00) - Gamma111*Gamma303*invPsi01 - Gamma101*Gamma313*invPsi01 - 
            Gamma112*Gamma303*invPsi02 - Gamma101*Gamma323*invPsi02 - Gamma111*Gamma313*invPsi11 - Gamma112*Gamma313*invPsi12 - 
            Gamma111*Gamma323*invPsi12 - Gamma112*Gamma323*invPsi22 + g11*Phi111*Phi133 + g12*Phi133*Phi211 + g12*Phi111*Phi233 + 
            g22*Phi211*Phi233 + g13*Phi133*Phi311 + g23*Phi233*Phi311 + g13*Phi111*Phi333 + g23*Phi211*Phi333 + g33*Phi311*Phi333 - 
            Pi11*Pi33) + invPsi23*(-(Gamma102*Gamma303*invPsi00) - Gamma112*Gamma303*invPsi01 - Gamma102*Gamma313*invPsi01 - 
            Gamma122*Gamma303*invPsi02 - Gamma102*Gamma323*invPsi02 - Gamma112*Gamma313*invPsi11 - Gamma122*Gamma313*invPsi12 - 
            Gamma112*Gamma323*invPsi12 - Gamma122*Gamma323*invPsi22 + g11*Phi112*Phi133 + g12*Phi133*Phi212 + g12*Phi112*Phi233 + 
            g22*Phi212*Phi233 + g13*Phi133*Phi312 + g23*Phi233*Phi312 + g13*Phi112*Phi333 + g23*Phi212*Phi333 + g33*Phi312*Phi333 - 
            Pi12*Pi33) + invPsi33*(-(Gamma103*Gamma303*invPsi00) - Gamma113*Gamma303*invPsi01 - Gamma103*Gamma313*invPsi01 - 
            Gamma123*Gamma303*invPsi02 - Gamma103*Gamma323*invPsi02 - Gamma113*Gamma313*invPsi11 - Gamma123*Gamma313*invPsi12 - 
            Gamma113*Gamma323*invPsi12 - Gamma123*Gamma323*invPsi22 + g11*Phi113*Phi133 + g12*Phi133*Phi213 + g12*Phi113*Phi233 + 
            g22*Phi213*Phi233 + g13*Phi133*Phi313 + g23*Phi233*Phi313 + g13*Phi113*Phi333 + g23*Phi213*Phi333 + g33*Phi313*Phi333 - 
            Pi13*Pi33));
    
    t1srcPi22 = 2*N*(invPsi00*(-(Power(Gamma200,2)*invPsi00) - 2*Gamma200*Gamma201*invPsi01 - 2*Gamma200*Gamma202*invPsi02 - 
            Power(Gamma201,2)*invPsi11 - 2*Gamma201*Gamma202*invPsi12 - Power(Gamma202,2)*invPsi22 + g11*Power(Phi102,2) + 
            2*g12*Phi102*Phi202 + g22*Power(Phi202,2) + 2*g13*Phi102*Phi302 + 2*g23*Phi202*Phi302 + g33*Power(Phi302,2) - 
            Power(Pi02,2)) + 2*invPsi01*(-(Gamma200*Gamma201*invPsi00) - Power(Gamma201,2)*invPsi01 - Gamma200*Gamma211*invPsi01 - 
            Gamma201*Gamma202*invPsi02 - Gamma200*Gamma212*invPsi02 - Gamma201*Gamma211*invPsi11 - Gamma202*Gamma211*invPsi12 - 
            Gamma201*Gamma212*invPsi12 - Gamma202*Gamma212*invPsi22 + g11*Phi102*Phi112 + g12*Phi112*Phi202 + g12*Phi102*Phi212 + 
            g22*Phi202*Phi212 + g13*Phi112*Phi302 + g23*Phi212*Phi302 + g13*Phi102*Phi312 + g23*Phi202*Phi312 + g33*Phi302*Phi312 - 
            Pi02*Pi12) + invPsi11*(-(Power(Gamma201,2)*invPsi00) - 2*Gamma201*Gamma211*invPsi01 - 2*Gamma201*Gamma212*invPsi02 - 
            Power(Gamma211,2)*invPsi11 - 2*Gamma211*Gamma212*invPsi12 - Power(Gamma212,2)*invPsi22 + g11*Power(Phi112,2) + 
            2*g12*Phi112*Phi212 + g22*Power(Phi212,2) + 2*g13*Phi112*Phi312 + 2*g23*Phi212*Phi312 + g33*Power(Phi312,2) - 
            Power(Pi12,2)) + 2*invPsi02*(-(Gamma200*Gamma202*invPsi00) - Gamma201*Gamma202*invPsi01 - Gamma200*Gamma212*invPsi01 - 
            Power(Gamma202,2)*invPsi02 - Gamma200*Gamma222*invPsi02 - Gamma201*Gamma212*invPsi11 - Gamma202*Gamma212*invPsi12 - 
            Gamma201*Gamma222*invPsi12 - Gamma202*Gamma222*invPsi22 + g11*Phi102*Phi122 + g12*Phi122*Phi202 + g12*Phi102*Phi222 + 
            g22*Phi202*Phi222 + g13*Phi122*Phi302 + g23*Phi222*Phi302 + g13*Phi102*Phi322 + g23*Phi202*Phi322 + g33*Phi302*Phi322 - 
            Pi02*Pi22) + 2*invPsi12*(-(Gamma201*Gamma202*invPsi00) - Gamma202*Gamma211*invPsi01 - Gamma201*Gamma212*invPsi01 - 
            Gamma202*Gamma212*invPsi02 - Gamma201*Gamma222*invPsi02 - Gamma211*Gamma212*invPsi11 - Power(Gamma212,2)*invPsi12 - 
            Gamma211*Gamma222*invPsi12 - Gamma212*Gamma222*invPsi22 + g11*Phi112*Phi122 + g12*Phi122*Phi212 + g12*Phi112*Phi222 + 
            g22*Phi212*Phi222 + g13*Phi122*Phi312 + g23*Phi222*Phi312 + g13*Phi112*Phi322 + g23*Phi212*Phi322 + g33*Phi312*Phi322 - 
            Pi12*Pi22) + invPsi22*(-(Power(Gamma202,2)*invPsi00) - 2*Gamma202*Gamma212*invPsi01 - 2*Gamma202*Gamma222*invPsi02 - 
            Power(Gamma212,2)*invPsi11 - 2*Gamma212*Gamma222*invPsi12 - Power(Gamma222,2)*invPsi22 + g11*Power(Phi122,2) + 
            2*g12*Phi122*Phi222 + g22*Power(Phi222,2) + 2*g13*Phi122*Phi322 + 2*g23*Phi222*Phi322 + g33*Power(Phi322,2) - 
            Power(Pi22,2)) + 2*invPsi03*(-(Gamma200*Gamma203*invPsi00) - Gamma201*Gamma203*invPsi01 - Gamma200*Gamma213*invPsi01 - 
            Gamma202*Gamma203*invPsi02 - Gamma200*Gamma223*invPsi02 - Gamma201*Gamma213*invPsi11 - Gamma202*Gamma213*invPsi12 - 
            Gamma201*Gamma223*invPsi12 - Gamma202*Gamma223*invPsi22 + g11*Phi102*Phi123 + g12*Phi123*Phi202 + g12*Phi102*Phi223 + 
            g22*Phi202*Phi223 + g13*Phi123*Phi302 + g23*Phi223*Phi302 + g13*Phi102*Phi323 + g23*Phi202*Phi323 + g33*Phi302*Phi323 - 
            Pi02*Pi23) + 2*invPsi13*(-(Gamma201*Gamma203*invPsi00) - Gamma203*Gamma211*invPsi01 - Gamma201*Gamma213*invPsi01 - 
            Gamma203*Gamma212*invPsi02 - Gamma201*Gamma223*invPsi02 - Gamma211*Gamma213*invPsi11 - Gamma212*Gamma213*invPsi12 - 
            Gamma211*Gamma223*invPsi12 - Gamma212*Gamma223*invPsi22 + g11*Phi112*Phi123 + g12*Phi123*Phi212 + g12*Phi112*Phi223 + 
            g22*Phi212*Phi223 + g13*Phi123*Phi312 + g23*Phi223*Phi312 + g13*Phi112*Phi323 + g23*Phi212*Phi323 + g33*Phi312*Phi323 - 
            Pi12*Pi23) + 2*invPsi23*(-(Gamma202*Gamma203*invPsi00) - Gamma203*Gamma212*invPsi01 - Gamma202*Gamma213*invPsi01 - 
            Gamma203*Gamma222*invPsi02 - Gamma202*Gamma223*invPsi02 - Gamma212*Gamma213*invPsi11 - Gamma213*Gamma222*invPsi12 - 
            Gamma212*Gamma223*invPsi12 - Gamma222*Gamma223*invPsi22 + g11*Phi122*Phi123 + g12*Phi123*Phi222 + g12*Phi122*Phi223 + 
            g22*Phi222*Phi223 + g13*Phi123*Phi322 + g23*Phi223*Phi322 + g13*Phi122*Phi323 + g23*Phi222*Phi323 + g33*Phi322*Phi323 - 
            Pi22*Pi23) + invPsi33*(-(Power(Gamma203,2)*invPsi00) - 2*Gamma203*Gamma213*invPsi01 - 2*Gamma203*Gamma223*invPsi02 - 
            Power(Gamma213,2)*invPsi11 - 2*Gamma213*Gamma223*invPsi12 - Power(Gamma223,2)*invPsi22 + g11*Power(Phi123,2) + 
            2*g12*Phi123*Phi223 + g22*Power(Phi223,2) + 2*g13*Phi123*Phi323 + 2*g23*Phi223*Phi323 + g33*Power(Phi323,2) - 
            Power(Pi23,2)));
    
    t1srcPi23 = 2*N*(invPsi00*(-(Gamma200*Gamma300*invPsi00) - Gamma201*Gamma300*invPsi01 - Gamma200*Gamma301*invPsi01 - 
            Gamma202*Gamma300*invPsi02 - Gamma200*Gamma302*invPsi02 - Gamma201*Gamma301*invPsi11 - Gamma202*Gamma301*invPsi12 - 
            Gamma201*Gamma302*invPsi12 - Gamma202*Gamma302*invPsi22 + g11*Phi102*Phi103 + g12*Phi103*Phi202 + g12*Phi102*Phi203 + 
            g22*Phi202*Phi203 + g13*Phi103*Phi302 + g23*Phi203*Phi302 + g13*Phi102*Phi303 + g23*Phi202*Phi303 + g33*Phi302*Phi303 - 
            Pi02*Pi03) + invPsi01*(-(Gamma201*Gamma300*invPsi00) - Gamma211*Gamma300*invPsi01 - Gamma201*Gamma301*invPsi01 - 
            Gamma212*Gamma300*invPsi02 - Gamma201*Gamma302*invPsi02 - Gamma211*Gamma301*invPsi11 - Gamma212*Gamma301*invPsi12 - 
            Gamma211*Gamma302*invPsi12 - Gamma212*Gamma302*invPsi22 + g11*Phi103*Phi112 + g12*Phi112*Phi203 + g12*Phi103*Phi212 + 
            g22*Phi203*Phi212 + g13*Phi112*Phi303 + g23*Phi212*Phi303 + g13*Phi103*Phi312 + g23*Phi203*Phi312 + g33*Phi303*Phi312 - 
            Pi03*Pi12) + invPsi01*(-(Gamma200*Gamma301*invPsi00) - Gamma201*Gamma301*invPsi01 - Gamma200*Gamma311*invPsi01 - 
            Gamma202*Gamma301*invPsi02 - Gamma200*Gamma312*invPsi02 - Gamma201*Gamma311*invPsi11 - Gamma202*Gamma311*invPsi12 - 
            Gamma201*Gamma312*invPsi12 - Gamma202*Gamma312*invPsi22 + g11*Phi102*Phi113 + g12*Phi113*Phi202 + g12*Phi102*Phi213 + 
            g22*Phi202*Phi213 + g13*Phi113*Phi302 + g23*Phi213*Phi302 + g13*Phi102*Phi313 + g23*Phi202*Phi313 + g33*Phi302*Phi313 - 
            Pi02*Pi13) + invPsi11*(-(Gamma201*Gamma301*invPsi00) - Gamma211*Gamma301*invPsi01 - Gamma201*Gamma311*invPsi01 - 
            Gamma212*Gamma301*invPsi02 - Gamma201*Gamma312*invPsi02 - Gamma211*Gamma311*invPsi11 - Gamma212*Gamma311*invPsi12 - 
            Gamma211*Gamma312*invPsi12 - Gamma212*Gamma312*invPsi22 + g11*Phi112*Phi113 + g12*Phi113*Phi212 + g12*Phi112*Phi213 + 
            g22*Phi212*Phi213 + g13*Phi113*Phi312 + g23*Phi213*Phi312 + g13*Phi112*Phi313 + g23*Phi212*Phi313 + g33*Phi312*Phi313 - 
            Pi12*Pi13) + invPsi02*(-(Gamma202*Gamma300*invPsi00) - Gamma212*Gamma300*invPsi01 - Gamma202*Gamma301*invPsi01 - 
            Gamma222*Gamma300*invPsi02 - Gamma202*Gamma302*invPsi02 - Gamma212*Gamma301*invPsi11 - Gamma222*Gamma301*invPsi12 - 
            Gamma212*Gamma302*invPsi12 - Gamma222*Gamma302*invPsi22 + g11*Phi103*Phi122 + g12*Phi122*Phi203 + g12*Phi103*Phi222 + 
            g22*Phi203*Phi222 + g13*Phi122*Phi303 + g23*Phi222*Phi303 + g13*Phi103*Phi322 + g23*Phi203*Phi322 + g33*Phi303*Phi322 - 
            Pi03*Pi22) + invPsi12*(-(Gamma202*Gamma301*invPsi00) - Gamma212*Gamma301*invPsi01 - Gamma202*Gamma311*invPsi01 - 
            Gamma222*Gamma301*invPsi02 - Gamma202*Gamma312*invPsi02 - Gamma212*Gamma311*invPsi11 - Gamma222*Gamma311*invPsi12 - 
            Gamma212*Gamma312*invPsi12 - Gamma222*Gamma312*invPsi22 + g11*Phi113*Phi122 + g12*Phi122*Phi213 + g12*Phi113*Phi222 + 
            g22*Phi213*Phi222 + g13*Phi122*Phi313 + g23*Phi222*Phi313 + g13*Phi113*Phi322 + g23*Phi213*Phi322 + g33*Phi313*Phi322 - 
            Pi13*Pi22) + invPsi02*(-(Gamma200*Gamma302*invPsi00) - Gamma201*Gamma302*invPsi01 - Gamma200*Gamma312*invPsi01 - 
            Gamma202*Gamma302*invPsi02 - Gamma200*Gamma322*invPsi02 - Gamma201*Gamma312*invPsi11 - Gamma202*Gamma312*invPsi12 - 
            Gamma201*Gamma322*invPsi12 - Gamma202*Gamma322*invPsi22 + g11*Phi102*Phi123 + g12*Phi123*Phi202 + g12*Phi102*Phi223 + 
            g22*Phi202*Phi223 + g13*Phi123*Phi302 + g23*Phi223*Phi302 + g13*Phi102*Phi323 + g23*Phi202*Phi323 + g33*Phi302*Phi323 - 
            Pi02*Pi23) + invPsi03*(-(Gamma203*Gamma300*invPsi00) - Gamma213*Gamma300*invPsi01 - Gamma203*Gamma301*invPsi01 - 
            Gamma223*Gamma300*invPsi02 - Gamma203*Gamma302*invPsi02 - Gamma213*Gamma301*invPsi11 - Gamma223*Gamma301*invPsi12 - 
            Gamma213*Gamma302*invPsi12 - Gamma223*Gamma302*invPsi22 + g11*Phi103*Phi123 + g12*Phi123*Phi203 + g12*Phi103*Phi223 + 
            g22*Phi203*Phi223 + g13*Phi123*Phi303 + g23*Phi223*Phi303 + g13*Phi103*Phi323 + g23*Phi203*Phi323 + g33*Phi303*Phi323 - 
            Pi03*Pi23) + invPsi12*(-(Gamma201*Gamma302*invPsi00) - Gamma211*Gamma302*invPsi01 - Gamma201*Gamma312*invPsi01 - 
            Gamma212*Gamma302*invPsi02 - Gamma201*Gamma322*invPsi02 - Gamma211*Gamma312*invPsi11 - Gamma212*Gamma312*invPsi12 - 
            Gamma211*Gamma322*invPsi12 - Gamma212*Gamma322*invPsi22 + g11*Phi112*Phi123 + g12*Phi123*Phi212 + g12*Phi112*Phi223 + 
            g22*Phi212*Phi223 + g13*Phi123*Phi312 + g23*Phi223*Phi312 + g13*Phi112*Phi323 + g23*Phi212*Phi323 + g33*Phi312*Phi323 - 
            Pi12*Pi23) + invPsi13*(-(Gamma203*Gamma301*invPsi00) - Gamma213*Gamma301*invPsi01 - Gamma203*Gamma311*invPsi01 - 
            Gamma223*Gamma301*invPsi02 - Gamma203*Gamma312*invPsi02 - Gamma213*Gamma311*invPsi11 - Gamma223*Gamma311*invPsi12 - 
            Gamma213*Gamma312*invPsi12 - Gamma223*Gamma312*invPsi22 + g11*Phi113*Phi123 + g12*Phi123*Phi213 + g12*Phi113*Phi223 + 
            g22*Phi213*Phi223 + g13*Phi123*Phi313 + g23*Phi223*Phi313 + g13*Phi113*Phi323 + g23*Phi213*Phi323 + g33*Phi313*Phi323 - 
            Pi13*Pi23) + invPsi22*(-(Gamma202*Gamma302*invPsi00) - Gamma212*Gamma302*invPsi01 - Gamma202*Gamma312*invPsi01 - 
            Gamma222*Gamma302*invPsi02 - Gamma202*Gamma322*invPsi02 - Gamma212*Gamma312*invPsi11 - Gamma222*Gamma312*invPsi12 - 
            Gamma212*Gamma322*invPsi12 - Gamma222*Gamma322*invPsi22 + g11*Phi122*Phi123 + g12*Phi123*Phi222 + g12*Phi122*Phi223 + 
            g22*Phi222*Phi223 + g13*Phi123*Phi322 + g23*Phi223*Phi322 + g13*Phi122*Phi323 + g23*Phi222*Phi323 + g33*Phi322*Phi323 - 
            Pi22*Pi23) + invPsi23*(-(Gamma203*Gamma302*invPsi00) - Gamma213*Gamma302*invPsi01 - Gamma203*Gamma312*invPsi01 - 
            Gamma223*Gamma302*invPsi02 - Gamma203*Gamma322*invPsi02 - Gamma213*Gamma312*invPsi11 - Gamma223*Gamma312*invPsi12 - 
            Gamma213*Gamma322*invPsi12 - Gamma223*Gamma322*invPsi22 + g11*Power(Phi123,2) + 2*g12*Phi123*Phi223 + 
            g22*Power(Phi223,2) + 2*g13*Phi123*Phi323 + 2*g23*Phi223*Phi323 + g33*Power(Phi323,2) - Power(Pi23,2)) + 
         invPsi03*(-(Gamma200*Gamma303*invPsi00) - Gamma201*Gamma303*invPsi01 - Gamma200*Gamma313*invPsi01 - 
            Gamma202*Gamma303*invPsi02 - Gamma200*Gamma323*invPsi02 - Gamma201*Gamma313*invPsi11 - Gamma202*Gamma313*invPsi12 - 
            Gamma201*Gamma323*invPsi12 - Gamma202*Gamma323*invPsi22 + g11*Phi102*Phi133 + g12*Phi133*Phi202 + g12*Phi102*Phi233 + 
            g22*Phi202*Phi233 + g13*Phi133*Phi302 + g23*Phi233*Phi302 + g13*Phi102*Phi333 + g23*Phi202*Phi333 + g33*Phi302*Phi333 - 
            Pi02*Pi33) + invPsi13*(-(Gamma201*Gamma303*invPsi00) - Gamma211*Gamma303*invPsi01 - Gamma201*Gamma313*invPsi01 - 
            Gamma212*Gamma303*invPsi02 - Gamma201*Gamma323*invPsi02 - Gamma211*Gamma313*invPsi11 - Gamma212*Gamma313*invPsi12 - 
            Gamma211*Gamma323*invPsi12 - Gamma212*Gamma323*invPsi22 + g11*Phi112*Phi133 + g12*Phi133*Phi212 + g12*Phi112*Phi233 + 
            g22*Phi212*Phi233 + g13*Phi133*Phi312 + g23*Phi233*Phi312 + g13*Phi112*Phi333 + g23*Phi212*Phi333 + g33*Phi312*Phi333 - 
            Pi12*Pi33) + invPsi23*(-(Gamma202*Gamma303*invPsi00) - Gamma212*Gamma303*invPsi01 - Gamma202*Gamma313*invPsi01 - 
            Gamma222*Gamma303*invPsi02 - Gamma202*Gamma323*invPsi02 - Gamma212*Gamma313*invPsi11 - Gamma222*Gamma313*invPsi12 - 
            Gamma212*Gamma323*invPsi12 - Gamma222*Gamma323*invPsi22 + g11*Phi122*Phi133 + g12*Phi133*Phi222 + g12*Phi122*Phi233 + 
            g22*Phi222*Phi233 + g13*Phi133*Phi322 + g23*Phi233*Phi322 + g13*Phi122*Phi333 + g23*Phi222*Phi333 + g33*Phi322*Phi333 - 
            Pi22*Pi33) + invPsi33*(-(Gamma203*Gamma303*invPsi00) - Gamma213*Gamma303*invPsi01 - Gamma203*Gamma313*invPsi01 - 
            Gamma223*Gamma303*invPsi02 - Gamma203*Gamma323*invPsi02 - Gamma213*Gamma313*invPsi11 - Gamma223*Gamma313*invPsi12 - 
            Gamma213*Gamma323*invPsi12 - Gamma223*Gamma323*invPsi22 + g11*Phi123*Phi133 + g12*Phi133*Phi223 + g12*Phi123*Phi233 + 
            g22*Phi223*Phi233 + g13*Phi133*Phi323 + g23*Phi233*Phi323 + g13*Phi123*Phi333 + g23*Phi223*Phi333 + g33*Phi323*Phi333 - 
            Pi23*Pi33));
    
    t1srcPi33 = 2*N*(invPsi00*(-(Power(Gamma300,2)*invPsi00) - 2*Gamma300*Gamma301*invPsi01 - 2*Gamma300*Gamma302*invPsi02 - 
            Power(Gamma301,2)*invPsi11 - 2*Gamma301*Gamma302*invPsi12 - Power(Gamma302,2)*invPsi22 + g11*Power(Phi103,2) + 
            2*g12*Phi103*Phi203 + g22*Power(Phi203,2) + 2*g13*Phi103*Phi303 + 2*g23*Phi203*Phi303 + g33*Power(Phi303,2) - 
            Power(Pi03,2)) + 2*invPsi01*(-(Gamma300*Gamma301*invPsi00) - Power(Gamma301,2)*invPsi01 - Gamma300*Gamma311*invPsi01 - 
            Gamma301*Gamma302*invPsi02 - Gamma300*Gamma312*invPsi02 - Gamma301*Gamma311*invPsi11 - Gamma302*Gamma311*invPsi12 - 
            Gamma301*Gamma312*invPsi12 - Gamma302*Gamma312*invPsi22 + g11*Phi103*Phi113 + g12*Phi113*Phi203 + g12*Phi103*Phi213 + 
            g22*Phi203*Phi213 + g13*Phi113*Phi303 + g23*Phi213*Phi303 + g13*Phi103*Phi313 + g23*Phi203*Phi313 + g33*Phi303*Phi313 - 
            Pi03*Pi13) + invPsi11*(-(Power(Gamma301,2)*invPsi00) - 2*Gamma301*Gamma311*invPsi01 - 2*Gamma301*Gamma312*invPsi02 - 
            Power(Gamma311,2)*invPsi11 - 2*Gamma311*Gamma312*invPsi12 - Power(Gamma312,2)*invPsi22 + g11*Power(Phi113,2) + 
            2*g12*Phi113*Phi213 + g22*Power(Phi213,2) + 2*g13*Phi113*Phi313 + 2*g23*Phi213*Phi313 + g33*Power(Phi313,2) - 
            Power(Pi13,2)) + 2*invPsi02*(-(Gamma300*Gamma302*invPsi00) - Gamma301*Gamma302*invPsi01 - Gamma300*Gamma312*invPsi01 - 
            Power(Gamma302,2)*invPsi02 - Gamma300*Gamma322*invPsi02 - Gamma301*Gamma312*invPsi11 - Gamma302*Gamma312*invPsi12 - 
            Gamma301*Gamma322*invPsi12 - Gamma302*Gamma322*invPsi22 + g11*Phi103*Phi123 + g12*Phi123*Phi203 + g12*Phi103*Phi223 + 
            g22*Phi203*Phi223 + g13*Phi123*Phi303 + g23*Phi223*Phi303 + g13*Phi103*Phi323 + g23*Phi203*Phi323 + g33*Phi303*Phi323 - 
            Pi03*Pi23) + 2*invPsi12*(-(Gamma301*Gamma302*invPsi00) - Gamma302*Gamma311*invPsi01 - Gamma301*Gamma312*invPsi01 - 
            Gamma302*Gamma312*invPsi02 - Gamma301*Gamma322*invPsi02 - Gamma311*Gamma312*invPsi11 - Power(Gamma312,2)*invPsi12 - 
            Gamma311*Gamma322*invPsi12 - Gamma312*Gamma322*invPsi22 + g11*Phi113*Phi123 + g12*Phi123*Phi213 + g12*Phi113*Phi223 + 
            g22*Phi213*Phi223 + g13*Phi123*Phi313 + g23*Phi223*Phi313 + g13*Phi113*Phi323 + g23*Phi213*Phi323 + g33*Phi313*Phi323 - 
            Pi13*Pi23) + invPsi22*(-(Power(Gamma302,2)*invPsi00) - 2*Gamma302*Gamma312*invPsi01 - 2*Gamma302*Gamma322*invPsi02 - 
            Power(Gamma312,2)*invPsi11 - 2*Gamma312*Gamma322*invPsi12 - Power(Gamma322,2)*invPsi22 + g11*Power(Phi123,2) + 
            2*g12*Phi123*Phi223 + g22*Power(Phi223,2) + 2*g13*Phi123*Phi323 + 2*g23*Phi223*Phi323 + g33*Power(Phi323,2) - 
            Power(Pi23,2)) + 2*invPsi03*(-(Gamma300*Gamma303*invPsi00) - Gamma301*Gamma303*invPsi01 - Gamma300*Gamma313*invPsi01 - 
            Gamma302*Gamma303*invPsi02 - Gamma300*Gamma323*invPsi02 - Gamma301*Gamma313*invPsi11 - Gamma302*Gamma313*invPsi12 - 
            Gamma301*Gamma323*invPsi12 - Gamma302*Gamma323*invPsi22 + g11*Phi103*Phi133 + g12*Phi133*Phi203 + g12*Phi103*Phi233 + 
            g22*Phi203*Phi233 + g13*Phi133*Phi303 + g23*Phi233*Phi303 + g13*Phi103*Phi333 + g23*Phi203*Phi333 + g33*Phi303*Phi333 - 
            Pi03*Pi33) + 2*invPsi13*(-(Gamma301*Gamma303*invPsi00) - Gamma303*Gamma311*invPsi01 - Gamma301*Gamma313*invPsi01 - 
            Gamma303*Gamma312*invPsi02 - Gamma301*Gamma323*invPsi02 - Gamma311*Gamma313*invPsi11 - Gamma312*Gamma313*invPsi12 - 
            Gamma311*Gamma323*invPsi12 - Gamma312*Gamma323*invPsi22 + g11*Phi113*Phi133 + g12*Phi133*Phi213 + g12*Phi113*Phi233 + 
            g22*Phi213*Phi233 + g13*Phi133*Phi313 + g23*Phi233*Phi313 + g13*Phi113*Phi333 + g23*Phi213*Phi333 + g33*Phi313*Phi333 - 
            Pi13*Pi33) + 2*invPsi23*(-(Gamma302*Gamma303*invPsi00) - Gamma303*Gamma312*invPsi01 - Gamma302*Gamma313*invPsi01 - 
            Gamma303*Gamma322*invPsi02 - Gamma302*Gamma323*invPsi02 - Gamma312*Gamma313*invPsi11 - Gamma313*Gamma322*invPsi12 - 
            Gamma312*Gamma323*invPsi12 - Gamma322*Gamma323*invPsi22 + g11*Phi123*Phi133 + g12*Phi133*Phi223 + g12*Phi123*Phi233 + 
            g22*Phi223*Phi233 + g13*Phi133*Phi323 + g23*Phi233*Phi323 + g13*Phi123*Phi333 + g23*Phi223*Phi333 + g33*Phi323*Phi333 - 
            Pi23*Pi33) + invPsi33*(-(Power(Gamma303,2)*invPsi00) - 2*Gamma303*Gamma313*invPsi01 - 2*Gamma303*Gamma323*invPsi02 - 
            Power(Gamma313,2)*invPsi11 - 2*Gamma313*Gamma323*invPsi12 - Power(Gamma323,2)*invPsi22 + g11*Power(Phi133,2) + 
            2*g12*Phi133*Phi233 + g22*Power(Phi233,2) + 2*g13*Phi133*Phi333 + 2*g23*Phi233*Phi333 + g33*Power(Phi333,2) - 
            Power(Pi33,2)));
    
    //TERM 2
    t2srcPi00 = -((2*dH00 - 2*Gamma000*H0*invPsi00 - 2*Gamma100*H0*invPsi01 - 2*Gamma000*H1*invPsi01 - 2*Gamma200*H0*invPsi02 - 
           2*Gamma000*H2*invPsi02 - 2*Gamma300*H0*invPsi03 - 2*Gamma000*H3*invPsi03 - 2*Gamma100*H1*invPsi11 - 
           2*Gamma200*H1*invPsi12 - 2*Gamma100*H2*invPsi12 - 2*Gamma300*H1*invPsi13 - 2*Gamma100*H3*invPsi13 - 
           2*Gamma200*H2*invPsi22 - 2*Gamma300*H2*invPsi23 - 2*Gamma200*H3*invPsi23 - 2*Gamma300*H3*invPsi33)*N);
    
    t2srcPi01 = -((dH01 + dH10 - 2*Gamma001*H0*invPsi00 - 2*Gamma101*H0*invPsi01 - 2*Gamma001*H1*invPsi01 - 2*Gamma201*H0*invPsi02 - 
           2*Gamma001*H2*invPsi02 - 2*Gamma301*H0*invPsi03 - 2*Gamma001*H3*invPsi03 - 2*Gamma101*H1*invPsi11 - 
           2*Gamma201*H1*invPsi12 - 2*Gamma101*H2*invPsi12 - 2*Gamma301*H1*invPsi13 - 2*Gamma101*H3*invPsi13 - 
           2*Gamma201*H2*invPsi22 - 2*Gamma301*H2*invPsi23 - 2*Gamma201*H3*invPsi23 - 2*Gamma301*H3*invPsi33)*N);
    
    t2srcPi02 = -((dH02 + dH20 - 2*Gamma002*H0*invPsi00 - 2*Gamma102*H0*invPsi01 - 2*Gamma002*H1*invPsi01 - 2*Gamma202*H0*invPsi02 - 
           2*Gamma002*H2*invPsi02 - 2*Gamma302*H0*invPsi03 - 2*Gamma002*H3*invPsi03 - 2*Gamma102*H1*invPsi11 - 
           2*Gamma202*H1*invPsi12 - 2*Gamma102*H2*invPsi12 - 2*Gamma302*H1*invPsi13 - 2*Gamma102*H3*invPsi13 - 
           2*Gamma202*H2*invPsi22 - 2*Gamma302*H2*invPsi23 - 2*Gamma202*H3*invPsi23 - 2*Gamma302*H3*invPsi33)*N);
    
    t2srcPi03 = -((dH03 + dH30 - 2*Gamma003*H0*invPsi00 - 2*Gamma103*H0*invPsi01 - 2*Gamma003*H1*invPsi01 - 2*Gamma203*H0*invPsi02 - 
           2*Gamma003*H2*invPsi02 - 2*Gamma303*H0*invPsi03 - 2*Gamma003*H3*invPsi03 - 2*Gamma103*H1*invPsi11 - 
           2*Gamma203*H1*invPsi12 - 2*Gamma103*H2*invPsi12 - 2*Gamma303*H1*invPsi13 - 2*Gamma103*H3*invPsi13 - 
           2*Gamma203*H2*invPsi22 - 2*Gamma303*H2*invPsi23 - 2*Gamma203*H3*invPsi23 - 2*Gamma303*H3*invPsi33)*N);
    
    t2srcPi11 = -((2*dH11 - 2*Gamma011*H0*invPsi00 - 2*Gamma111*H0*invPsi01 - 2*Gamma011*H1*invPsi01 - 2*Gamma211*H0*invPsi02 - 
           2*Gamma011*H2*invPsi02 - 2*Gamma311*H0*invPsi03 - 2*Gamma011*H3*invPsi03 - 2*Gamma111*H1*invPsi11 - 
           2*Gamma211*H1*invPsi12 - 2*Gamma111*H2*invPsi12 - 2*Gamma311*H1*invPsi13 - 2*Gamma111*H3*invPsi13 - 
           2*Gamma211*H2*invPsi22 - 2*Gamma311*H2*invPsi23 - 2*Gamma211*H3*invPsi23 - 2*Gamma311*H3*invPsi33)*N);
    
    t2srcPi12 = -((dH12 + dH21 - 2*Gamma012*H0*invPsi00 - 2*Gamma112*H0*invPsi01 - 2*Gamma012*H1*invPsi01 - 2*Gamma212*H0*invPsi02 - 
           2*Gamma012*H2*invPsi02 - 2*Gamma312*H0*invPsi03 - 2*Gamma012*H3*invPsi03 - 2*Gamma112*H1*invPsi11 - 
           2*Gamma212*H1*invPsi12 - 2*Gamma112*H2*invPsi12 - 2*Gamma312*H1*invPsi13 - 2*Gamma112*H3*invPsi13 - 
           2*Gamma212*H2*invPsi22 - 2*Gamma312*H2*invPsi23 - 2*Gamma212*H3*invPsi23 - 2*Gamma312*H3*invPsi33)*N);
    
    t2srcPi13 = -((dH13 + dH31 - 2*Gamma013*H0*invPsi00 - 2*Gamma113*H0*invPsi01 - 2*Gamma013*H1*invPsi01 - 2*Gamma213*H0*invPsi02 - 
           2*Gamma013*H2*invPsi02 - 2*Gamma313*H0*invPsi03 - 2*Gamma013*H3*invPsi03 - 2*Gamma113*H1*invPsi11 - 
           2*Gamma213*H1*invPsi12 - 2*Gamma113*H2*invPsi12 - 2*Gamma313*H1*invPsi13 - 2*Gamma113*H3*invPsi13 - 
           2*Gamma213*H2*invPsi22 - 2*Gamma313*H2*invPsi23 - 2*Gamma213*H3*invPsi23 - 2*Gamma313*H3*invPsi33)*N);
    
    t2srcPi22 = -((2*dH22 - 2*Gamma022*H0*invPsi00 - 2*Gamma122*H0*invPsi01 - 2*Gamma022*H1*invPsi01 - 2*Gamma222*H0*invPsi02 - 
           2*Gamma022*H2*invPsi02 - 2*Gamma322*H0*invPsi03 - 2*Gamma022*H3*invPsi03 - 2*Gamma122*H1*invPsi11 - 
           2*Gamma222*H1*invPsi12 - 2*Gamma122*H2*invPsi12 - 2*Gamma322*H1*invPsi13 - 2*Gamma122*H3*invPsi13 - 
           2*Gamma222*H2*invPsi22 - 2*Gamma322*H2*invPsi23 - 2*Gamma222*H3*invPsi23 - 2*Gamma322*H3*invPsi33)*N);
    
    t2srcPi23 = -((dH23 + dH32 - 2*Gamma023*H0*invPsi00 - 2*Gamma123*H0*invPsi01 - 2*Gamma023*H1*invPsi01 - 2*Gamma223*H0*invPsi02 - 
           2*Gamma023*H2*invPsi02 - 2*Gamma323*H0*invPsi03 - 2*Gamma023*H3*invPsi03 - 2*Gamma123*H1*invPsi11 - 
           2*Gamma223*H1*invPsi12 - 2*Gamma123*H2*invPsi12 - 2*Gamma323*H1*invPsi13 - 2*Gamma123*H3*invPsi13 - 
           2*Gamma223*H2*invPsi22 - 2*Gamma323*H2*invPsi23 - 2*Gamma223*H3*invPsi23 - 2*Gamma323*H3*invPsi33)*N);
    
    t2srcPi33 = -((2*dH33 - 2*Gamma033*H0*invPsi00 - 2*Gamma133*H0*invPsi01 - 2*Gamma033*H1*invPsi01 - 2*Gamma233*H0*invPsi02 - 
           2*Gamma033*H2*invPsi02 - 2*Gamma333*H0*invPsi03 - 2*Gamma033*H3*invPsi03 - 2*Gamma133*H1*invPsi11 - 
           2*Gamma233*H1*invPsi12 - 2*Gamma133*H2*invPsi12 - 2*Gamma333*H1*invPsi13 - 2*Gamma133*H3*invPsi13 - 
           2*Gamma233*H2*invPsi22 - 2*Gamma333*H2*invPsi23 - 2*Gamma233*H3*invPsi23 - 2*Gamma333*H3*invPsi33)*N);
    
    //TERM 3
    
    t3srcPi00 = -0.5*N*(Power(Pi00,2)*Power(t0,2) + 2*Pi00*Pi01*t0*t1 + Pi00*Pi11*Power(t1,2) + 2*Pi00*Pi02*t0*t2 + 
           2*Pi00*Pi12*t1*t2 + Pi00*Pi22*Power(t2,2) + 2*Pi00*Pi03*t0*t3 + 2*Pi00*Pi13*t1*t3 + 2*Pi00*Pi23*t2*t3 + 
           Pi00*Pi33*Power(t3,2)); 
    
    t3srcPi01 = -0.5*N*(Pi00*Pi01*Power(t0,2) + 2*Power(Pi01,2)*t0*t1 + Pi01*Pi11*Power(t1,2) + 
           2*Pi01*Pi02*t0*t2 + 2*Pi01*Pi12*t1*t2 + Pi01*Pi22*Power(t2,2) + 2*Pi01*Pi03*t0*t3 + 2*Pi01*Pi13*t1*t3 + 
           2*Pi01*Pi23*t2*t3 + Pi01*Pi33*Power(t3,2));
    
    t3srcPi02 = -0.5*N*(Pi00*Pi02*Power(t0,2) + 2*Pi01*Pi02*t0*t1 + Pi02*Pi11*Power(t1,2) + 2*Power(Pi02,2)*t0*t2 + 2*Pi02*Pi12*t1*t2 + 
           Pi02*Pi22*Power(t2,2) + 2*Pi02*Pi03*t0*t3 + 2*Pi02*Pi13*t1*t3 + 2*Pi02*Pi23*t2*t3 + Pi02*Pi33*Power(t3,2));
    
    t3srcPi03 = -0.5*N*(Pi00*Pi03*Power(t0,2) + 2*Pi01*Pi03*t0*t1 + Pi03*Pi11*Power(t1,2) + 2*Pi02*Pi03*t0*t2 + 2*Pi03*Pi12*t1*t2 + 
           Pi03*Pi22*Power(t2,2) + 2*Power(Pi03,2)*t0*t3 + 2*Pi03*Pi13*t1*t3 + 2*Pi03*Pi23*t2*t3 + Pi03*Pi33*Power(t3,2));
    
    t3srcPi11 = -0.5*N*(Pi00*Pi11*Power(t0,2) + 2*Pi01*Pi11*t0*t1 + Power(Pi11,2)*Power(t1,2) + 2*Pi02*Pi11*t0*t2 + 2*Pi11*Pi12*t1*t2 + 
           Pi11*Pi22*Power(t2,2) + 2*Pi03*Pi11*t0*t3 + 2*Pi11*Pi13*t1*t3 + 2*Pi11*Pi23*t2*t3 + Pi11*Pi33*Power(t3,2));
    
    t3srcPi12 = -0.5*N*(Pi00*Pi12*Power(t0,2) + 2*Pi01*Pi12*t0*t1 + Pi11*Pi12*Power(t1,2) + 2*Pi02*Pi12*t0*t2 + 2*Power(Pi12,2)*t1*t2 + 
           Pi12*Pi22*Power(t2,2) + 2*Pi03*Pi12*t0*t3 + 2*Pi12*Pi13*t1*t3 + 2*Pi12*Pi23*t2*t3 + Pi12*Pi33*Power(t3,2));
    
    t3srcPi13 = -0.5*N*(Pi00*Pi13*Power(t0,2) + 2*Pi01*Pi13*t0*t1 + Pi11*Pi13*Power(t1,2) + 2*Pi02*Pi13*t0*t2 + 2*Pi12*Pi13*t1*t2 + 
           Pi13*Pi22*Power(t2,2) + 2*Pi03*Pi13*t0*t3 + 2*Power(Pi13,2)*t1*t3 + 2*Pi13*Pi23*t2*t3 + Pi13*Pi33*Power(t3,2));
    
    t3srcPi22 = -0.5*N*(Pi00*Pi22*Power(t0,2) + 2*Pi01*Pi22*t0*t1 + Pi11*Pi22*Power(t1,2) + 2*Pi02*Pi22*t0*t2 + 2*Pi12*Pi22*t1*t2 + 
           Power(Pi22,2)*Power(t2,2) + 2*Pi03*Pi22*t0*t3 + 2*Pi13*Pi22*t1*t3 + 2*Pi22*Pi23*t2*t3 + Pi22*Pi33*Power(t3,2));
    
    t3srcPi23 = -0.5*N*(Pi00*Pi23*Power(t0,2) + 2*Pi01*Pi23*t0*t1 + Pi11*Pi23*Power(t1,2) + 2*Pi02*Pi23*t0*t2 + 2*Pi12*Pi23*t1*t2 + 
           Pi22*Pi23*Power(t2,2) + 2*Pi03*Pi23*t0*t3 + 2*Pi13*Pi23*t1*t3 + 2*Power(Pi23,2)*t2*t3 + Pi23*Pi33*Power(t3,2));
    
    t3srcPi33 = -0.5*N*(Pi00*Pi33*Power(t0,2) + 2*Pi01*Pi33*t0*t1 + Pi11*Pi33*Power(t1,2) + 2*Pi02*Pi33*t0*t2 + 2*Pi12*Pi33*t1*t2 + 
           Pi22*Pi33*Power(t2,2) + 2*Pi03*Pi33*t0*t3 + 2*Pi13*Pi33*t1*t3 + 2*Pi23*Pi33*t2*t3 + Power(Pi33,2)*Power(t3,2));
    
    //TERM 4
    t4srcPi00 = -(N*(g11*Phi100*Pi00*t0 + g12*Phi200*Pi00*t0 + g13*Phi300*Pi00*t0 + g12*Phi100*Pi01*t0 + g22*Phi200*Pi01*t0 + 
             g23*Phi300*Pi01*t0 + g13*Phi100*Pi02*t0 + g23*Phi200*Pi02*t0 + g33*Phi300*Pi02*t0 + g11*Phi100*Pi01*t1 + 
             g12*Phi200*Pi01*t1 + g13*Phi300*Pi01*t1 + g12*Phi100*Pi11*t1 + g22*Phi200*Pi11*t1 + g23*Phi300*Pi11*t1 + 
             g13*Phi100*Pi12*t1 + g23*Phi200*Pi12*t1 + g33*Phi300*Pi12*t1 + g11*Phi100*Pi02*t2 + g12*Phi200*Pi02*t2 + 
             g13*Phi300*Pi02*t2 + g12*Phi100*Pi12*t2 + g22*Phi200*Pi12*t2 + g23*Phi300*Pi12*t2 + g13*Phi100*Pi22*t2 + 
             g23*Phi200*Pi22*t2 + g33*Phi300*Pi22*t2 + g11*Phi100*Pi03*t3 + g12*Phi200*Pi03*t3 + g13*Phi300*Pi03*t3 + 
             g12*Phi100*Pi13*t3 + g22*Phi200*Pi13*t3 + g23*Phi300*Pi13*t3 + g13*Phi100*Pi23*t3 + g23*Phi200*Pi23*t3 + 
             g33*Phi300*Pi23*t3));
    
    t4srcPi01 = -(N*(g11*Phi101*Pi00*t0 + g12*Phi201*Pi00*t0 + g13*Phi301*Pi00*t0 + g12*Phi101*Pi01*t0 + 
             g22*Phi201*Pi01*t0 + g23*Phi301*Pi01*t0 + g13*Phi101*Pi02*t0 + g23*Phi201*Pi02*t0 + g33*Phi301*Pi02*t0 + 
             g11*Phi101*Pi01*t1 + g12*Phi201*Pi01*t1 + g13*Phi301*Pi01*t1 + g12*Phi101*Pi11*t1 + g22*Phi201*Pi11*t1 + 
             g23*Phi301*Pi11*t1 + g13*Phi101*Pi12*t1 + g23*Phi201*Pi12*t1 + g33*Phi301*Pi12*t1 + g11*Phi101*Pi02*t2 + 
             g12*Phi201*Pi02*t2 + g13*Phi301*Pi02*t2 + g12*Phi101*Pi12*t2 + g22*Phi201*Pi12*t2 + g23*Phi301*Pi12*t2 + 
             g13*Phi101*Pi22*t2 + g23*Phi201*Pi22*t2 + g33*Phi301*Pi22*t2 + g11*Phi101*Pi03*t3 + g12*Phi201*Pi03*t3 + 
             g13*Phi301*Pi03*t3 + g12*Phi101*Pi13*t3 + g22*Phi201*Pi13*t3 + g23*Phi301*Pi13*t3 + g13*Phi101*Pi23*t3 + 
             g23*Phi201*Pi23*t3 + g33*Phi301*Pi23*t3));
                
    t4srcPi02 = -(N*(g11*Phi102*Pi00*t0 + g12*Phi202*Pi00*t0 + g13*Phi302*Pi00*t0 + g12*Phi102*Pi01*t0 + g22*Phi202*Pi01*t0 + 
             g23*Phi302*Pi01*t0 + g13*Phi102*Pi02*t0 + g23*Phi202*Pi02*t0 + g33*Phi302*Pi02*t0 + g11*Phi102*Pi01*t1 + 
             g12*Phi202*Pi01*t1 + g13*Phi302*Pi01*t1 + g12*Phi102*Pi11*t1 + g22*Phi202*Pi11*t1 + g23*Phi302*Pi11*t1 + 
             g13*Phi102*Pi12*t1 + g23*Phi202*Pi12*t1 + g33*Phi302*Pi12*t1 + g11*Phi102*Pi02*t2 + g12*Phi202*Pi02*t2 + 
             g13*Phi302*Pi02*t2 + g12*Phi102*Pi12*t2 + g22*Phi202*Pi12*t2 + g23*Phi302*Pi12*t2 + g13*Phi102*Pi22*t2 + 
             g23*Phi202*Pi22*t2 + g33*Phi302*Pi22*t2 + g11*Phi102*Pi03*t3 + g12*Phi202*Pi03*t3 + g13*Phi302*Pi03*t3 + 
             g12*Phi102*Pi13*t3 + g22*Phi202*Pi13*t3 + g23*Phi302*Pi13*t3 + g13*Phi102*Pi23*t3 + g23*Phi202*Pi23*t3 + 
             g33*Phi302*Pi23*t3));
                
    t4srcPi03 = -(N*(g11*Phi103*Pi00*t0 + g12*Phi203*Pi00*t0 + g13*Phi303*Pi00*t0 + g12*Phi103*Pi01*t0 + 
             g22*Phi203*Pi01*t0 + g23*Phi303*Pi01*t0 + g13*Phi103*Pi02*t0 + g23*Phi203*Pi02*t0 + g33*Phi303*Pi02*t0 + 
             g11*Phi103*Pi01*t1 + g12*Phi203*Pi01*t1 + g13*Phi303*Pi01*t1 + g12*Phi103*Pi11*t1 + g22*Phi203*Pi11*t1 + 
             g23*Phi303*Pi11*t1 + g13*Phi103*Pi12*t1 + g23*Phi203*Pi12*t1 + g33*Phi303*Pi12*t1 + g11*Phi103*Pi02*t2 + 
             g12*Phi203*Pi02*t2 + g13*Phi303*Pi02*t2 + g12*Phi103*Pi12*t2 + g22*Phi203*Pi12*t2 + g23*Phi303*Pi12*t2 + 
             g13*Phi103*Pi22*t2 + g23*Phi203*Pi22*t2 + g33*Phi303*Pi22*t2 + g11*Phi103*Pi03*t3 + g12*Phi203*Pi03*t3 + 
             g13*Phi303*Pi03*t3 + g12*Phi103*Pi13*t3 + g22*Phi203*Pi13*t3 + g23*Phi303*Pi13*t3 + g13*Phi103*Pi23*t3 + 
             g23*Phi203*Pi23*t3 + g33*Phi303*Pi23*t3));
    
    t4srcPi11 = -(N*(g11*Phi111*Pi00*t0 + g12*Phi211*Pi00*t0 + g13*Phi311*Pi00*t0 + g12*Phi111*Pi01*t0 + 
             g22*Phi211*Pi01*t0 + g23*Phi311*Pi01*t0 + g13*Phi111*Pi02*t0 + g23*Phi211*Pi02*t0 + g33*Phi311*Pi02*t0 + 
             g11*Phi111*Pi01*t1 + g12*Phi211*Pi01*t1 + g13*Phi311*Pi01*t1 + g12*Phi111*Pi11*t1 + g22*Phi211*Pi11*t1 + 
             g23*Phi311*Pi11*t1 + g13*Phi111*Pi12*t1 + g23*Phi211*Pi12*t1 + g33*Phi311*Pi12*t1 + g11*Phi111*Pi02*t2 + 
             g12*Phi211*Pi02*t2 + g13*Phi311*Pi02*t2 + g12*Phi111*Pi12*t2 + g22*Phi211*Pi12*t2 + g23*Phi311*Pi12*t2 + 
             g13*Phi111*Pi22*t2 + g23*Phi211*Pi22*t2 + g33*Phi311*Pi22*t2 + g11*Phi111*Pi03*t3 + g12*Phi211*Pi03*t3 + 
             g13*Phi311*Pi03*t3 + g12*Phi111*Pi13*t3 + g22*Phi211*Pi13*t3 + g23*Phi311*Pi13*t3 + g13*Phi111*Pi23*t3 + 
             g23*Phi211*Pi23*t3 + g33*Phi311*Pi23*t3));
    
    t4srcPi12 = -(N*(g11*Phi112*Pi00*t0 + g12*Phi212*Pi00*t0 + g13*Phi312*Pi00*t0 + g12*Phi112*Pi01*t0 + g22*Phi212*Pi01*t0 + 
             g23*Phi312*Pi01*t0 + g13*Phi112*Pi02*t0 + g23*Phi212*Pi02*t0 + g33*Phi312*Pi02*t0 + g11*Phi112*Pi01*t1 + 
             g12*Phi212*Pi01*t1 + g13*Phi312*Pi01*t1 + g12*Phi112*Pi11*t1 + g22*Phi212*Pi11*t1 + g23*Phi312*Pi11*t1 + 
             g13*Phi112*Pi12*t1 + g23*Phi212*Pi12*t1 + g33*Phi312*Pi12*t1 + g11*Phi112*Pi02*t2 + g12*Phi212*Pi02*t2 + 
             g13*Phi312*Pi02*t2 + g12*Phi112*Pi12*t2 + g22*Phi212*Pi12*t2 + g23*Phi312*Pi12*t2 + g13*Phi112*Pi22*t2 + 
             g23*Phi212*Pi22*t2 + g33*Phi312*Pi22*t2 + g11*Phi112*Pi03*t3 + g12*Phi212*Pi03*t3 + g13*Phi312*Pi03*t3 + 
             g12*Phi112*Pi13*t3 + g22*Phi212*Pi13*t3 + g23*Phi312*Pi13*t3 + g13*Phi112*Pi23*t3 + g23*Phi212*Pi23*t3 + 
             g33*Phi312*Pi23*t3)); 
    
    t4srcPi13 = -(N*(g11*Phi113*Pi00*t0 + g12*Phi213*Pi00*t0 + g13*Phi313*Pi00*t0 + g12*Phi113*Pi01*t0 + 
             g22*Phi213*Pi01*t0 + g23*Phi313*Pi01*t0 + g13*Phi113*Pi02*t0 + g23*Phi213*Pi02*t0 + g33*Phi313*Pi02*t0 + 
             g11*Phi113*Pi01*t1 + g12*Phi213*Pi01*t1 + g13*Phi313*Pi01*t1 + g12*Phi113*Pi11*t1 + g22*Phi213*Pi11*t1 + 
             g23*Phi313*Pi11*t1 + g13*Phi113*Pi12*t1 + g23*Phi213*Pi12*t1 + g33*Phi313*Pi12*t1 + g11*Phi113*Pi02*t2 + 
             g12*Phi213*Pi02*t2 + g13*Phi313*Pi02*t2 + g12*Phi113*Pi12*t2 + g22*Phi213*Pi12*t2 + g23*Phi313*Pi12*t2 + 
             g13*Phi113*Pi22*t2 + g23*Phi213*Pi22*t2 + g33*Phi313*Pi22*t2 + g11*Phi113*Pi03*t3 + g12*Phi213*Pi03*t3 + 
             g13*Phi313*Pi03*t3 + g12*Phi113*Pi13*t3 + g22*Phi213*Pi13*t3 + g23*Phi313*Pi13*t3 + g13*Phi113*Pi23*t3 + 
             g23*Phi213*Pi23*t3 + g33*Phi313*Pi23*t3));
    


    t4srcPi22 = -(N*(g11*Phi122*Pi00*t0 + g12*Phi222*Pi00*t0 + g13*Phi322*Pi00*t0 + g12*Phi122*Pi01*t0 + g22*Phi222*Pi01*t0 + 
             g23*Phi322*Pi01*t0 + g13*Phi122*Pi02*t0 + g23*Phi222*Pi02*t0 + g33*Phi322*Pi02*t0 + g11*Phi122*Pi01*t1 + 
             g12*Phi222*Pi01*t1 + g13*Phi322*Pi01*t1 + g12*Phi122*Pi11*t1 + g22*Phi222*Pi11*t1 + g23*Phi322*Pi11*t1 + 
             g13*Phi122*Pi12*t1 + g23*Phi222*Pi12*t1 + g33*Phi322*Pi12*t1 + g11*Phi122*Pi02*t2 + g12*Phi222*Pi02*t2 + 
             g13*Phi322*Pi02*t2 + g12*Phi122*Pi12*t2 + g22*Phi222*Pi12*t2 + g23*Phi322*Pi12*t2 + g13*Phi122*Pi22*t2 + 
             g23*Phi222*Pi22*t2 + g33*Phi322*Pi22*t2 + g11*Phi122*Pi03*t3 + g12*Phi222*Pi03*t3 + g13*Phi322*Pi03*t3 + 
             g12*Phi122*Pi13*t3 + g22*Phi222*Pi13*t3 + g23*Phi322*Pi13*t3 + g13*Phi122*Pi23*t3 + g23*Phi222*Pi23*t3 + 
             g33*Phi322*Pi23*t3));
    
    t4srcPi23 = -(N*(g11*Phi123*Pi00*t0 + g12*Phi223*Pi00*t0 + g13*Phi323*Pi00*t0 + g12*Phi123*Pi01*t0 + 
             g22*Phi223*Pi01*t0 + g23*Phi323*Pi01*t0 + g13*Phi123*Pi02*t0 + g23*Phi223*Pi02*t0 + g33*Phi323*Pi02*t0 + 
             g11*Phi123*Pi01*t1 + g12*Phi223*Pi01*t1 + g13*Phi323*Pi01*t1 + g12*Phi123*Pi11*t1 + g22*Phi223*Pi11*t1 + 
             g23*Phi323*Pi11*t1 + g13*Phi123*Pi12*t1 + g23*Phi223*Pi12*t1 + g33*Phi323*Pi12*t1 + g11*Phi123*Pi02*t2 + 
             g12*Phi223*Pi02*t2 + g13*Phi323*Pi02*t2 + g12*Phi123*Pi12*t2 + g22*Phi223*Pi12*t2 + g23*Phi323*Pi12*t2 + 
             g13*Phi123*Pi22*t2 + g23*Phi223*Pi22*t2 + g33*Phi323*Pi22*t2 + g11*Phi123*Pi03*t3 + g12*Phi223*Pi03*t3 + 
             g13*Phi323*Pi03*t3 + g12*Phi123*Pi13*t3 + g22*Phi223*Pi13*t3 + g23*Phi323*Pi13*t3 + g13*Phi123*Pi23*t3 + 
             g23*Phi223*Pi23*t3 + g33*Phi323*Pi23*t3));
    
    t4srcPi33 = -(N*(g11*Phi133*Pi00*t0 + g12*Phi233*Pi00*t0 + g13*Phi333*Pi00*t0 + g12*Phi133*Pi01*t0 + 
             g22*Phi233*Pi01*t0 + g23*Phi333*Pi01*t0 + g13*Phi133*Pi02*t0 + g23*Phi233*Pi02*t0 + g33*Phi333*Pi02*t0 + 
             g11*Phi133*Pi01*t1 + g12*Phi233*Pi01*t1 + g13*Phi333*Pi01*t1 + g12*Phi133*Pi11*t1 + g22*Phi233*Pi11*t1 + 
             g23*Phi333*Pi11*t1 + g13*Phi133*Pi12*t1 + g23*Phi233*Pi12*t1 + g33*Phi333*Pi12*t1 + g11*Phi133*Pi02*t2 + 
             g12*Phi233*Pi02*t2 + g13*Phi333*Pi02*t2 + g12*Phi133*Pi12*t2 + g22*Phi233*Pi12*t2 + g23*Phi333*Pi12*t2 + 
             g13*Phi133*Pi22*t2 + g23*Phi233*Pi22*t2 + g33*Phi333*Pi22*t2 + g11*Phi133*Pi03*t3 + g12*Phi233*Pi03*t3 + 
             g13*Phi333*Pi03*t3 + g12*Phi133*Pi13*t3 + g22*Phi233*Pi13*t3 + g23*Phi333*Pi13*t3 + g13*Phi133*Pi23*t3 + 
             g23*Phi233*Pi23*t3 + g33*Phi333*Pi23*t3));
    
    //TERM 5
    
    t5srcPi00 = gamma0*N*((2*t0 - Psi00*t0)*(H0 + vecGamma0) - Psi00*t1*(H1 + vecGamma1) - Psi00*t2*(H2 + vecGamma2) - 
           Psi00*t3*(H3 + vecGamma3));
    
    t5srcPi01 = gamma0*N*((-(Psi01*t0) + t1)*(H0 + vecGamma0) + (t0 - Psi01*t1)*(H1 + vecGamma1) - 
           Psi01*t2*(H2 + vecGamma2) - Psi01*t3*(H3 + vecGamma3));
    
    t5srcPi02 = gamma0*N*((-(Psi02*t0) + t2)*(H0 + vecGamma0) - Psi02*t1*(H1 + vecGamma1) + (t0 - Psi02*t2)*(H2 + vecGamma2) - 
           Psi02*t3*(H3 + vecGamma3));
    
    t5srcPi03 = gamma0*N*((-(Psi03*t0) + t3)*(H0 + vecGamma0) - Psi03*t1*(H1 + vecGamma1) - 
           Psi03*t2*(H2 + vecGamma2) + (t0 - Psi03*t3)*(H3 + vecGamma3));
    
    
    t5srcPi11 = gamma0*N*(-(Psi11*t0*(H0 + vecGamma0)) + (2*t1 - Psi11*t1)*(H1 + vecGamma1) - 
           Psi11*t2*(H2 + vecGamma2) - Psi11*t3*(H3 + vecGamma3));
     
    
    t5srcPi12 = gamma0*N*(-(Psi12*t0*(H0 + vecGamma0)) + (-(Psi12*t1) + t2)*(H1 + vecGamma1) + (t1 - Psi12*t2)*(H2 + vecGamma2) - 
           Psi12*t3*(H3 + vecGamma3));
    
    
    t5srcPi13 = gamma0*N*(-(Psi13*t0*(H0 + vecGamma0)) + (-(Psi13*t1) + t3)*(H1 + vecGamma1) - 
           Psi13*t2*(H2 + vecGamma2) + (t1 - Psi13*t3)*(H3 + vecGamma3));
    
    t5srcPi22 = gamma0*N*(-(Psi22*t0*(H0 + vecGamma0)) - Psi22*t1*(H1 + vecGamma1) + (2*t2 - Psi22*t2)*(H2 + vecGamma2) - 
           Psi22*t3*(H3 + vecGamma3));
    
    t5srcPi23 = gamma0*N*(-(Psi23*t0*(H0 + vecGamma0)) - Psi23*t1*(H1 + vecGamma1) + 
           (-(Psi23*t2) + t3)*(H2 + vecGamma2) + (t2 - Psi23*t3)*(H3 + vecGamma3));
     
    t5srcPi33 = gamma0*N*(-(Psi33*t0*(H0 + vecGamma0)) - Psi33*t1*(H1 + vecGamma1) - 
           Psi33*t2*(H2 + vecGamma2) + (2*t3 - Psi33*t3)*(H3 + vecGamma3));
    
    //TERM 6
    t6srcPi00 = -(gamma1*gamma2*(N1*Phi100 + N2*Phi200 + N3*Phi300));
    
    t6srcPi01 = -(gamma1*gamma2*(N1*Phi101 + N2*Phi201 + N3*Phi301));
    
    t6srcPi02 = -(gamma1*gamma2*(N1*Phi102 + N2*Phi202 + N3*Phi302));
    
    t6srcPi03 = -(gamma1*gamma2*(N1*Phi103 + N2*Phi203 + N3*Phi303));
    
    t6srcPi11 = -(gamma1*gamma2*(N1*Phi111 + N2*Phi211 + N3*Phi311));
    
    t6srcPi12 = -(gamma1*gamma2*(N1*Phi112 + N2*Phi212 + N3*Phi312));
    
    t6srcPi13 = -(gamma1*gamma2*(N1*Phi113 + N2*Phi213 + N3*Phi313));
    
    t6srcPi22 = -(gamma1*gamma2*(N1*Phi122 + N2*Phi222 + N3*Phi322));
    
    t6srcPi23 = -(gamma1*gamma2*(N1*Phi123 + N2*Phi223 + N3*Phi323));
    
    t6srcPi33 = -(gamma1*gamma2*(N1*Phi133 + N2*Phi233 + N3*Phi333));
    
    //Sum all terms of RHS Pi
    srcPi00 = t1srcPi00 + t2srcPi00 + t3srcPi00 + t4srcPi00 + t5srcPi00 + t6srcPi00;
    srcPi01 = t1srcPi01 + t2srcPi01 + t3srcPi01 + t4srcPi01 + t5srcPi01 + t6srcPi01;
    srcPi02 = t1srcPi02 + t2srcPi02 + t3srcPi02 + t4srcPi02 + t5srcPi02 + t6srcPi02;
    srcPi03 = t1srcPi03 + t2srcPi03 + t3srcPi03 + t4srcPi03 + t5srcPi03 + t6srcPi03;
    srcPi11 = t1srcPi11 + t2srcPi11 + t3srcPi11 + t4srcPi11 + t5srcPi11 + t6srcPi11;
    srcPi12 = t1srcPi12 + t2srcPi12 + t3srcPi12 + t4srcPi12 + t5srcPi12 + t6srcPi12;
    srcPi13 = t1srcPi13 + t2srcPi13 + t3srcPi13 + t4srcPi13 + t5srcPi13 + t6srcPi13;
    srcPi22 = t1srcPi22 + t2srcPi22 + t3srcPi22 + t4srcPi22 + t5srcPi22 + t6srcPi22;
    srcPi23 = t1srcPi23 + t2srcPi23 + t3srcPi23 + t4srcPi23 + t5srcPi23 + t6srcPi23;
    srcPi33 = t1srcPi33 + t2srcPi33 + t3srcPi33 + t4srcPi33 + t5srcPi33 + t6srcPi33;
    
    //RHS of Phi
    
    //TEMR 1 of Phi
    
    t1srcPhi100 = 0.5*N*Pi00*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi101 =  0.5*N*Pi01*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi102 = 0.5*N*Pi02*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi103 = 0.5*N*Pi03*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi111 = 0.5*N*Pi11*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi112 = 0.5*N*Pi12*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
        
    t1srcPhi113 = 0.5*N*Pi13*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi122 = 0.5*N*Pi22*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi123 = 0.5*N*Pi23*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi133 = 0.5*N*Pi33*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
            Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
    t1srcPhi200 = 0.5*N*Pi00*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi201 = 0.5*N*Pi01*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi202 = 0.5*N*Pi02*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi203 = 0.5*N*Pi03*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi211 = 0.5*N*Pi11*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi212 = 0.5*N*Pi12*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi213 = 0.5*N*Pi13*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi222 = 0.5*N*Pi22*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi223 = 0.5*N*Pi23*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi233 = 0.5*N*Pi33*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
            Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));
    
    t1srcPhi300 = 0.5*N*Pi00*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
            Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi301 = 0.5*N*Pi01*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
            Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi302 = 0.5*N*Pi02*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
            Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi303 = 0.5*N*Pi03*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
            Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi311 = 0.5*N*Pi11*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
             Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi312 = 0.5*N*Pi12*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                    Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi313 = 0.5*N*Pi13*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                    Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi322 = 0.5*N*Pi22*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                    Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi323 = 0.5*N*Pi23*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                    Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    t1srcPhi333 = 0.5*N*Pi33*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                    Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));
    
    
    //TERM 2 of Phi
    t2srcPhi100 =  N*(g11*Phi100*Phi101*t0 + 
                     g12*Phi100*Phi102*t0 + g13*Phi100*Phi103*t0 + 
                     g12*Phi101*Phi200*t0 + g22*Phi102*Phi200*t0 + 
                     g23*Phi103*Phi200*t0 + g13*Phi101*Phi300*t0 + 
                     g23*Phi102*Phi300*t0 + g33*Phi103*Phi300*t0 + 
                     g11*Phi100*Phi111*t1 + g12*Phi100*Phi112*t1 + 
                     g13*Phi100*Phi113*t1 + g12*Phi111*Phi200*t1 + 
                     g22*Phi112*Phi200*t1 + g23*Phi113*Phi200*t1 + 
                     g13*Phi111*Phi300*t1 + g23*Phi112*Phi300*t1 + 
                     g33*Phi113*Phi300*t1 + g11*Phi100*Phi112*t2 + 
                     g12*Phi100*Phi122*t2 + g13*Phi100*Phi123*t2 + 
                     g12*Phi112*Phi200*t2 + g22*Phi122*Phi200*t2 + 
                     g23*Phi123*Phi200*t2 + g13*Phi112*Phi300*t2 + 
                     g23*Phi122*Phi300*t2 + g33*Phi123*Phi300*t2 + 
                     g11*Phi100*Phi113*t3 + g12*Phi100*Phi123*t3 + 
                     g13*Phi100*Phi133*t3 + g12*Phi113*Phi200*t3 + 
                     g22*Phi123*Phi200*t3 + g23*Phi133*Phi200*t3 + 
                     g13*Phi113*Phi300*t3 + g23*Phi123*Phi300*t3 + 
                     g33*Phi133*Phi300*t3);   

    t2srcPhi101 = N*(g11*Power(Phi101,2)*t0 + g12*Phi101*Phi102*t0 + 
                     g13*Phi101*Phi103*t0 + g12*Phi101*Phi201*t0 + 
                     g22*Phi102*Phi201*t0 + g23*Phi103*Phi201*t0 + 
                     g13*Phi101*Phi301*t0 + g23*Phi102*Phi301*t0 + 
                     g33*Phi103*Phi301*t0 + g11*Phi101*Phi111*t1 + 
                     g12*Phi101*Phi112*t1 + g13*Phi101*Phi113*t1 + 
                     g12*Phi111*Phi201*t1 + g22*Phi112*Phi201*t1 + 
                     g23*Phi113*Phi201*t1 + g13*Phi111*Phi301*t1 + 
                     g23*Phi112*Phi301*t1 + g33*Phi113*Phi301*t1 + 
                     g11*Phi101*Phi112*t2 + g12*Phi101*Phi122*t2 + 
                     g13*Phi101*Phi123*t2 + g12*Phi112*Phi201*t2 + 
                     g22*Phi122*Phi201*t2 + g23*Phi123*Phi201*t2 + 
                     g13*Phi112*Phi301*t2 + g23*Phi122*Phi301*t2 + 
                     g33*Phi123*Phi301*t2 + g11*Phi101*Phi113*t3 + 
                     g12*Phi101*Phi123*t3 + g13*Phi101*Phi133*t3 + 
                     g12*Phi113*Phi201*t3 + g22*Phi123*Phi201*t3 + 
                     g23*Phi133*Phi201*t3 + g13*Phi113*Phi301*t3 + 
                     g23*Phi123*Phi301*t3 + g33*Phi133*Phi301*t3); 
    
    t2srcPhi102 = N*(g11*Phi101*Phi102*t0 + g12*Power(Phi102,2)*t0 + 
                     g13*Phi102*Phi103*t0 + g12*Phi101*Phi202*t0 + 
                     g22*Phi102*Phi202*t0 + g23*Phi103*Phi202*t0 + 
                     g13*Phi101*Phi302*t0 + g23*Phi102*Phi302*t0 + 
                     g33*Phi103*Phi302*t0 + g11*Phi102*Phi111*t1 + 
                     g12*Phi102*Phi112*t1 + g13*Phi102*Phi113*t1 + 
                     g12*Phi111*Phi202*t1 + g22*Phi112*Phi202*t1 + 
                     g23*Phi113*Phi202*t1 + g13*Phi111*Phi302*t1 + 
                     g23*Phi112*Phi302*t1 + g33*Phi113*Phi302*t1 + 
                     g11*Phi102*Phi112*t2 + g12*Phi102*Phi122*t2 + 
                     g13*Phi102*Phi123*t2 + g12*Phi112*Phi202*t2 + 
                     g22*Phi122*Phi202*t2 + g23*Phi123*Phi202*t2 + 
                     g13*Phi112*Phi302*t2 + g23*Phi122*Phi302*t2 + 
                     g33*Phi123*Phi302*t2 + g11*Phi102*Phi113*t3 + 
                     g12*Phi102*Phi123*t3 + g13*Phi102*Phi133*t3 + 
                     g12*Phi113*Phi202*t3 + g22*Phi123*Phi202*t3 + 
                     g23*Phi133*Phi202*t3 + g13*Phi113*Phi302*t3 + 
                     g23*Phi123*Phi302*t3 + g33*Phi133*Phi302*t3);    

    t2srcPhi103 = N*(g11*Phi101*Phi103*t0 + g12*Phi102*Phi103*t0 + 
                     g13*Power(Phi103,2)*t0 + g12*Phi101*Phi203*t0 + 
                     g22*Phi102*Phi203*t0 + g23*Phi103*Phi203*t0 + 
                     g13*Phi101*Phi303*t0 + g23*Phi102*Phi303*t0 + 
                     g33*Phi103*Phi303*t0 + g11*Phi103*Phi111*t1 + 
                     g12*Phi103*Phi112*t1 + g13*Phi103*Phi113*t1 + 
                     g12*Phi111*Phi203*t1 + g22*Phi112*Phi203*t1 + 
                     g23*Phi113*Phi203*t1 + g13*Phi111*Phi303*t1 + 
                     g23*Phi112*Phi303*t1 + g33*Phi113*Phi303*t1 + 
                     g11*Phi103*Phi112*t2 + g12*Phi103*Phi122*t2 + 
                     g13*Phi103*Phi123*t2 + g12*Phi112*Phi203*t2 + 
                     g22*Phi122*Phi203*t2 + g23*Phi123*Phi203*t2 + 
                     g13*Phi112*Phi303*t2 + g23*Phi122*Phi303*t2 + 
                     g33*Phi123*Phi303*t2 + g11*Phi103*Phi113*t3 + 
                     g12*Phi103*Phi123*t3 + g13*Phi103*Phi133*t3 + 
                     g12*Phi113*Phi203*t3 + g22*Phi123*Phi203*t3 + 
                     g23*Phi133*Phi203*t3 + g13*Phi113*Phi303*t3 + 
                     g23*Phi123*Phi303*t3 + g33*Phi133*Phi303*t3);
         
    t2srcPhi111 = N*(g11*Phi101*Phi111*t0 + g12*Phi102*Phi111*t0 + 
     g13*Phi103*Phi111*t0 + g12*Phi101*Phi211*t0 + 
     g22*Phi102*Phi211*t0 + g23*Phi103*Phi211*t0 + 
     g13*Phi101*Phi311*t0 + g23*Phi102*Phi311*t0 + 
     g33*Phi103*Phi311*t0 + g11*Power(Phi111,2)*t1 + 
     g12*Phi111*Phi112*t1 + g13*Phi111*Phi113*t1 + 
     g12*Phi111*Phi211*t1 + g22*Phi112*Phi211*t1 + 
     g23*Phi113*Phi211*t1 + g13*Phi111*Phi311*t1 + 
     g23*Phi112*Phi311*t1 + g33*Phi113*Phi311*t1 + 
     g11*Phi111*Phi112*t2 + g12*Phi111*Phi122*t2 + 
     g13*Phi111*Phi123*t2 + g12*Phi112*Phi211*t2 + 
     g22*Phi122*Phi211*t2 + g23*Phi123*Phi211*t2 + 
     g13*Phi112*Phi311*t2 + g23*Phi122*Phi311*t2 + 
     g33*Phi123*Phi311*t2 + g11*Phi111*Phi113*t3 + 
     g12*Phi111*Phi123*t3 + g13*Phi111*Phi133*t3 + 
     g12*Phi113*Phi211*t3 + g22*Phi123*Phi211*t3 + 
     g23*Phi133*Phi211*t3 + g13*Phi113*Phi311*t3 + 
     g23*Phi123*Phi311*t3 + g33*Phi133*Phi311*t3);  

    t2srcPhi112 = N*(g11*Phi101*Phi112*t0 + g12*Phi102*Phi112*t0 + 
     g13*Phi103*Phi112*t0 + g12*Phi101*Phi212*t0 + 
     g22*Phi102*Phi212*t0 + g23*Phi103*Phi212*t0 + 
     g13*Phi101*Phi312*t0 + g23*Phi102*Phi312*t0 + 
     g33*Phi103*Phi312*t0 + g11*Phi111*Phi112*t1 + 
     g12*Power(Phi112,2)*t1 + g13*Phi112*Phi113*t1 + 
     g12*Phi111*Phi212*t1 + g22*Phi112*Phi212*t1 + 
     g23*Phi113*Phi212*t1 + g13*Phi111*Phi312*t1 + 
     g23*Phi112*Phi312*t1 + g33*Phi113*Phi312*t1 + 
     g11*Power(Phi112,2)*t2 + g12*Phi112*Phi122*t2 + 
     g13*Phi112*Phi123*t2 + g12*Phi112*Phi212*t2 + 
     g22*Phi122*Phi212*t2 + g23*Phi123*Phi212*t2 + 
     g13*Phi112*Phi312*t2 + g23*Phi122*Phi312*t2 + 
     g33*Phi123*Phi312*t2 + g11*Phi112*Phi113*t3 + 
     g12*Phi112*Phi123*t3 + g13*Phi112*Phi133*t3 + 
     g12*Phi113*Phi212*t3 + g22*Phi123*Phi212*t3 + 
     g23*Phi133*Phi212*t3 + g13*Phi113*Phi312*t3 + 
     g23*Phi123*Phi312*t3 + g33*Phi133*Phi312*t3);
    
    t2srcPhi113 = N*(g11*Phi101*Phi113*t0 + g12*Phi102*Phi113*t0 + 
     g13*Phi103*Phi113*t0 + g12*Phi101*Phi213*t0 + 
     g22*Phi102*Phi213*t0 + g23*Phi103*Phi213*t0 + 
     g13*Phi101*Phi313*t0 + g23*Phi102*Phi313*t0 + 
     g33*Phi103*Phi313*t0 + g11*Phi111*Phi113*t1 + 
     g12*Phi112*Phi113*t1 + g13*Power(Phi113,2)*t1 + 
     g12*Phi111*Phi213*t1 + g22*Phi112*Phi213*t1 + 
     g23*Phi113*Phi213*t1 + g13*Phi111*Phi313*t1 + 
     g23*Phi112*Phi313*t1 + g33*Phi113*Phi313*t1 + 
     g11*Phi112*Phi113*t2 + g12*Phi113*Phi122*t2 + 
     g13*Phi113*Phi123*t2 + g12*Phi112*Phi213*t2 + 
     g22*Phi122*Phi213*t2 + g23*Phi123*Phi213*t2 + 
     g13*Phi112*Phi313*t2 + g23*Phi122*Phi313*t2 + 
     g33*Phi123*Phi313*t2 + g11*Power(Phi113,2)*t3 + 
     g12*Phi113*Phi123*t3 + g13*Phi113*Phi133*t3 + 
     g12*Phi113*Phi213*t3 + g22*Phi123*Phi213*t3 + 
     g23*Phi133*Phi213*t3 + g13*Phi113*Phi313*t3 + 
     g23*Phi123*Phi313*t3 + g33*Phi133*Phi313*t3); 
    
    t2srcPhi122 = N*(g11*Phi101*Phi122*t0 + g12*Phi102*Phi122*t0 + 
     g13*Phi103*Phi122*t0 + g12*Phi101*Phi222*t0 + 
     g22*Phi102*Phi222*t0 + g23*Phi103*Phi222*t0 + 
     g13*Phi101*Phi322*t0 + g23*Phi102*Phi322*t0 + 
     g33*Phi103*Phi322*t0 + g11*Phi111*Phi122*t1 + 
     g12*Phi112*Phi122*t1 + g13*Phi113*Phi122*t1 + 
     g12*Phi111*Phi222*t1 + g22*Phi112*Phi222*t1 + 
     g23*Phi113*Phi222*t1 + g13*Phi111*Phi322*t1 + 
     g23*Phi112*Phi322*t1 + g33*Phi113*Phi322*t1 + 
     g11*Phi112*Phi122*t2 + g12*Power(Phi122,2)*t2 + 
     g13*Phi122*Phi123*t2 + g12*Phi112*Phi222*t2 + 
     g22*Phi122*Phi222*t2 + g23*Phi123*Phi222*t2 + 
     g13*Phi112*Phi322*t2 + g23*Phi122*Phi322*t2 + 
     g33*Phi123*Phi322*t2 + g11*Phi113*Phi122*t3 + 
     g12*Phi122*Phi123*t3 + g13*Phi122*Phi133*t3 + 
     g12*Phi113*Phi222*t3 + g22*Phi123*Phi222*t3 + 
     g23*Phi133*Phi222*t3 + g13*Phi113*Phi322*t3 + 
     g23*Phi123*Phi322*t3 + g33*Phi133*Phi322*t3); 
    
    t2srcPhi123 = N*(g11*Phi101*Phi123*t0 + g12*Phi102*Phi123*t0 + 
     g13*Phi103*Phi123*t0 + g12*Phi101*Phi223*t0 + 
     g22*Phi102*Phi223*t0 + g23*Phi103*Phi223*t0 + 
     g13*Phi101*Phi323*t0 + g23*Phi102*Phi323*t0 + 
     g33*Phi103*Phi323*t0 + g11*Phi111*Phi123*t1 + 
     g12*Phi112*Phi123*t1 + g13*Phi113*Phi123*t1 + 
     g12*Phi111*Phi223*t1 + g22*Phi112*Phi223*t1 + 
     g23*Phi113*Phi223*t1 + g13*Phi111*Phi323*t1 + 
     g23*Phi112*Phi323*t1 + g33*Phi113*Phi323*t1 + 
     g11*Phi112*Phi123*t2 + g12*Phi122*Phi123*t2 + 
     g13*Power(Phi123,2)*t2 + g12*Phi112*Phi223*t2 + 
     g22*Phi122*Phi223*t2 + g23*Phi123*Phi223*t2 + 
     g13*Phi112*Phi323*t2 + g23*Phi122*Phi323*t2 + 
     g33*Phi123*Phi323*t2 + g11*Phi113*Phi123*t3 + 
     g12*Power(Phi123,2)*t3 + g13*Phi123*Phi133*t3 + 
     g12*Phi113*Phi223*t3 + g22*Phi123*Phi223*t3 + 
     g23*Phi133*Phi223*t3 + g13*Phi113*Phi323*t3 + 
     g23*Phi123*Phi323*t3 + g33*Phi133*Phi323*t3); 
    
    t2srcPhi133 = N*(g11*Phi101*Phi133*t0 + g12*Phi102*Phi133*t0 + 
     g13*Phi103*Phi133*t0 + g12*Phi101*Phi233*t0 + 
     g22*Phi102*Phi233*t0 + g23*Phi103*Phi233*t0 + 
     g13*Phi101*Phi333*t0 + g23*Phi102*Phi333*t0 + 
     g33*Phi103*Phi333*t0 + g11*Phi111*Phi133*t1 + 
     g12*Phi112*Phi133*t1 + g13*Phi113*Phi133*t1 + 
     g12*Phi111*Phi233*t1 + g22*Phi112*Phi233*t1 + 
     g23*Phi113*Phi233*t1 + g13*Phi111*Phi333*t1 + 
     g23*Phi112*Phi333*t1 + g33*Phi113*Phi333*t1 + 
     g11*Phi112*Phi133*t2 + g12*Phi122*Phi133*t2 + 
     g13*Phi123*Phi133*t2 + g12*Phi112*Phi233*t2 + 
     g22*Phi122*Phi233*t2 + g23*Phi123*Phi233*t2 + 
     g13*Phi112*Phi333*t2 + g23*Phi122*Phi333*t2 + 
     g33*Phi123*Phi333*t2 + g11*Phi113*Phi133*t3 + 
     g12*Phi123*Phi133*t3 + g13*Power(Phi133,2)*t3 + 
     g12*Phi113*Phi233*t3 + g22*Phi123*Phi233*t3 + 
     g23*Phi133*Phi233*t3 + g13*Phi113*Phi333*t3 + 
     g23*Phi123*Phi333*t3 + g33*Phi133*Phi333*t3); 
    
    t2srcPhi200 = N*(g11*Phi100*Phi201*t0 + g12*Phi200*Phi201*t0 + 
     g12*Phi100*Phi202*t0 + g22*Phi200*Phi202*t0 + 
     g13*Phi100*Phi203*t0 + g23*Phi200*Phi203*t0 + 
     g13*Phi201*Phi300*t0 + g23*Phi202*Phi300*t0 + 
     g33*Phi203*Phi300*t0 + g11*Phi100*Phi211*t1 + 
     g12*Phi200*Phi211*t1 + g12*Phi100*Phi212*t1 + 
     g22*Phi200*Phi212*t1 + g13*Phi100*Phi213*t1 + 
     g23*Phi200*Phi213*t1 + g13*Phi211*Phi300*t1 + 
     g23*Phi212*Phi300*t1 + g33*Phi213*Phi300*t1 + 
     g11*Phi100*Phi212*t2 + g12*Phi200*Phi212*t2 + 
     g12*Phi100*Phi222*t2 + g22*Phi200*Phi222*t2 + 
     g13*Phi100*Phi223*t2 + g23*Phi200*Phi223*t2 + 
     g13*Phi212*Phi300*t2 + g23*Phi222*Phi300*t2 + 
     g33*Phi223*Phi300*t2 + g11*Phi100*Phi213*t3 + 
     g12*Phi200*Phi213*t3 + g12*Phi100*Phi223*t3 + 
     g22*Phi200*Phi223*t3 + g13*Phi100*Phi233*t3 + 
     g23*Phi200*Phi233*t3 + g13*Phi213*Phi300*t3 + 
     g23*Phi223*Phi300*t3 + g33*Phi233*Phi300*t3); 
    
    t2srcPhi201 = N*(g11*Phi101*Phi201*t0 + g12*Power(Phi201,2)*t0 + 
     g12*Phi101*Phi202*t0 + g22*Phi201*Phi202*t0 + 
     g13*Phi101*Phi203*t0 + g23*Phi201*Phi203*t0 + 
     g13*Phi201*Phi301*t0 + g23*Phi202*Phi301*t0 + 
     g33*Phi203*Phi301*t0 + g11*Phi101*Phi211*t1 + 
     g12*Phi201*Phi211*t1 + g12*Phi101*Phi212*t1 + 
     g22*Phi201*Phi212*t1 + g13*Phi101*Phi213*t1 + 
     g23*Phi201*Phi213*t1 + g13*Phi211*Phi301*t1 + 
     g23*Phi212*Phi301*t1 + g33*Phi213*Phi301*t1 + 
     g11*Phi101*Phi212*t2 + g12*Phi201*Phi212*t2 + 
     g12*Phi101*Phi222*t2 + g22*Phi201*Phi222*t2 + 
     g13*Phi101*Phi223*t2 + g23*Phi201*Phi223*t2 + 
     g13*Phi212*Phi301*t2 + g23*Phi222*Phi301*t2 + 
     g33*Phi223*Phi301*t2 + g11*Phi101*Phi213*t3 + 
     g12*Phi201*Phi213*t3 + g12*Phi101*Phi223*t3 + 
     g22*Phi201*Phi223*t3 + g13*Phi101*Phi233*t3 + 
     g23*Phi201*Phi233*t3 + g13*Phi213*Phi301*t3 + 
     g23*Phi223*Phi301*t3 + g33*Phi233*Phi301*t3); 
    
    t2srcPhi202 = N*(g11*Phi102*Phi201*t0 + g12*Phi102*Phi202*t0 + 
     g12*Phi201*Phi202*t0 + g22*Power(Phi202,2)*t0 + 
     g13*Phi102*Phi203*t0 + g23*Phi202*Phi203*t0 + 
     g13*Phi201*Phi302*t0 + g23*Phi202*Phi302*t0 + 
     g33*Phi203*Phi302*t0 + g11*Phi102*Phi211*t1 + 
     g12*Phi202*Phi211*t1 + g12*Phi102*Phi212*t1 + 
     g22*Phi202*Phi212*t1 + g13*Phi102*Phi213*t1 + 
     g23*Phi202*Phi213*t1 + g13*Phi211*Phi302*t1 + 
     g23*Phi212*Phi302*t1 + g33*Phi213*Phi302*t1 + 
     g11*Phi102*Phi212*t2 + g12*Phi202*Phi212*t2 + 
     g12*Phi102*Phi222*t2 + g22*Phi202*Phi222*t2 + 
     g13*Phi102*Phi223*t2 + g23*Phi202*Phi223*t2 + 
     g13*Phi212*Phi302*t2 + g23*Phi222*Phi302*t2 + 
     g33*Phi223*Phi302*t2 + g11*Phi102*Phi213*t3 + 
     g12*Phi202*Phi213*t3 + g12*Phi102*Phi223*t3 + 
     g22*Phi202*Phi223*t3 + g13*Phi102*Phi233*t3 + 
     g23*Phi202*Phi233*t3 + g13*Phi213*Phi302*t3 + 
     g23*Phi223*Phi302*t3 + g33*Phi233*Phi302*t3); 
    
    t2srcPhi203 = N*(g11*Phi103*Phi201*t0 + g12*Phi103*Phi202*t0 + 
     g13*Phi103*Phi203*t0 + g12*Phi201*Phi203*t0 + 
     g22*Phi202*Phi203*t0 + g23*Power(Phi203,2)*t0 + 
     g13*Phi201*Phi303*t0 + g23*Phi202*Phi303*t0 + 
     g33*Phi203*Phi303*t0 + g11*Phi103*Phi211*t1 + 
     g12*Phi203*Phi211*t1 + g12*Phi103*Phi212*t1 + 
     g22*Phi203*Phi212*t1 + g13*Phi103*Phi213*t1 + 
     g23*Phi203*Phi213*t1 + g13*Phi211*Phi303*t1 + 
     g23*Phi212*Phi303*t1 + g33*Phi213*Phi303*t1 + 
     g11*Phi103*Phi212*t2 + g12*Phi203*Phi212*t2 + 
     g12*Phi103*Phi222*t2 + g22*Phi203*Phi222*t2 + 
     g13*Phi103*Phi223*t2 + g23*Phi203*Phi223*t2 + 
     g13*Phi212*Phi303*t2 + g23*Phi222*Phi303*t2 + 
     g33*Phi223*Phi303*t2 + g11*Phi103*Phi213*t3 + 
     g12*Phi203*Phi213*t3 + g12*Phi103*Phi223*t3 + 
     g22*Phi203*Phi223*t3 + g13*Phi103*Phi233*t3 + 
     g23*Phi203*Phi233*t3 + g13*Phi213*Phi303*t3 + 
     g23*Phi223*Phi303*t3 + g33*Phi233*Phi303*t3); 
    
    t2srcPhi211 = N*(g11*Phi111*Phi201*t0 + g12*Phi111*Phi202*t0 + 
     g13*Phi111*Phi203*t0 + g12*Phi201*Phi211*t0 + 
     g22*Phi202*Phi211*t0 + g23*Phi203*Phi211*t0 + 
     g13*Phi201*Phi311*t0 + g23*Phi202*Phi311*t0 + 
     g33*Phi203*Phi311*t0 + g11*Phi111*Phi211*t1 + 
     g12*Power(Phi211,2)*t1 + g12*Phi111*Phi212*t1 + 
     g22*Phi211*Phi212*t1 + g13*Phi111*Phi213*t1 + 
     g23*Phi211*Phi213*t1 + g13*Phi211*Phi311*t1 + 
     g23*Phi212*Phi311*t1 + g33*Phi213*Phi311*t1 + 
     g11*Phi111*Phi212*t2 + g12*Phi211*Phi212*t2 + 
     g12*Phi111*Phi222*t2 + g22*Phi211*Phi222*t2 + 
     g13*Phi111*Phi223*t2 + g23*Phi211*Phi223*t2 + 
     g13*Phi212*Phi311*t2 + g23*Phi222*Phi311*t2 + 
     g33*Phi223*Phi311*t2 + g11*Phi111*Phi213*t3 + 
     g12*Phi211*Phi213*t3 + g12*Phi111*Phi223*t3 + 
     g22*Phi211*Phi223*t3 + g13*Phi111*Phi233*t3 + 
     g23*Phi211*Phi233*t3 + g13*Phi213*Phi311*t3 + 
     g23*Phi223*Phi311*t3 + g33*Phi233*Phi311*t3); 

    t2srcPhi212 = N*(g11*Phi112*Phi201*t0 + g12*Phi112*Phi202*t0 + 
     g13*Phi112*Phi203*t0 + g12*Phi201*Phi212*t0 + 
     g22*Phi202*Phi212*t0 + g23*Phi203*Phi212*t0 + 
     g13*Phi201*Phi312*t0 + g23*Phi202*Phi312*t0 + 
     g33*Phi203*Phi312*t0 + g11*Phi112*Phi211*t1 + 
     g12*Phi112*Phi212*t1 + g12*Phi211*Phi212*t1 + 
     g22*Power(Phi212,2)*t1 + g13*Phi112*Phi213*t1 + 
     g23*Phi212*Phi213*t1 + g13*Phi211*Phi312*t1 + 
     g23*Phi212*Phi312*t1 + g33*Phi213*Phi312*t1 + 
     g11*Phi112*Phi212*t2 + g12*Power(Phi212,2)*t2 + 
     g12*Phi112*Phi222*t2 + g22*Phi212*Phi222*t2 + 
     g13*Phi112*Phi223*t2 + g23*Phi212*Phi223*t2 + 
     g13*Phi212*Phi312*t2 + g23*Phi222*Phi312*t2 + 
     g33*Phi223*Phi312*t2 + g11*Phi112*Phi213*t3 + 
     g12*Phi212*Phi213*t3 + g12*Phi112*Phi223*t3 + 
     g22*Phi212*Phi223*t3 + g13*Phi112*Phi233*t3 + 
     g23*Phi212*Phi233*t3 + g13*Phi213*Phi312*t3 + 
     g23*Phi223*Phi312*t3 + g33*Phi233*Phi312*t3); 
    
    t2srcPhi213 = N*(g11*Phi113*Phi201*t0 + g12*Phi113*Phi202*t0 + 
     g13*Phi113*Phi203*t0 + g12*Phi201*Phi213*t0 + 
     g22*Phi202*Phi213*t0 + g23*Phi203*Phi213*t0 + 
     g13*Phi201*Phi313*t0 + g23*Phi202*Phi313*t0 + 
     g33*Phi203*Phi313*t0 + g11*Phi113*Phi211*t1 + 
     g12*Phi113*Phi212*t1 + g13*Phi113*Phi213*t1 + 
     g12*Phi211*Phi213*t1 + g22*Phi212*Phi213*t1 + 
     g23*Power(Phi213,2)*t1 + g13*Phi211*Phi313*t1 + 
     g23*Phi212*Phi313*t1 + g33*Phi213*Phi313*t1 + 
     g11*Phi113*Phi212*t2 + g12*Phi212*Phi213*t2 + 
     g12*Phi113*Phi222*t2 + g22*Phi213*Phi222*t2 + 
     g13*Phi113*Phi223*t2 + g23*Phi213*Phi223*t2 + 
     g13*Phi212*Phi313*t2 + g23*Phi222*Phi313*t2 + 
     g33*Phi223*Phi313*t2 + g11*Phi113*Phi213*t3 + 
     g12*Power(Phi213,2)*t3 + g12*Phi113*Phi223*t3 + 
     g22*Phi213*Phi223*t3 + g13*Phi113*Phi233*t3 + 
     g23*Phi213*Phi233*t3 + g13*Phi213*Phi313*t3 + 
     g23*Phi223*Phi313*t3 + g33*Phi233*Phi313*t3);
    
    t2srcPhi222 = N*(g11*Phi122*Phi201*t0 + g12*Phi122*Phi202*t0 + 
     g13*Phi122*Phi203*t0 + g12*Phi201*Phi222*t0 + 
     g22*Phi202*Phi222*t0 + g23*Phi203*Phi222*t0 + 
     g13*Phi201*Phi322*t0 + g23*Phi202*Phi322*t0 + 
     g33*Phi203*Phi322*t0 + g11*Phi122*Phi211*t1 + 
     g12*Phi122*Phi212*t1 + g13*Phi122*Phi213*t1 + 
     g12*Phi211*Phi222*t1 + g22*Phi212*Phi222*t1 + 
     g23*Phi213*Phi222*t1 + g13*Phi211*Phi322*t1 + 
     g23*Phi212*Phi322*t1 + g33*Phi213*Phi322*t1 + 
     g11*Phi122*Phi212*t2 + g12*Phi122*Phi222*t2 + 
     g12*Phi212*Phi222*t2 + g22*Power(Phi222,2)*t2 + 
     g13*Phi122*Phi223*t2 + g23*Phi222*Phi223*t2 + 
     g13*Phi212*Phi322*t2 + g23*Phi222*Phi322*t2 + 
     g33*Phi223*Phi322*t2 + g11*Phi122*Phi213*t3 + 
     g12*Phi213*Phi222*t3 + g12*Phi122*Phi223*t3 + 
     g22*Phi222*Phi223*t3 + g13*Phi122*Phi233*t3 + 
     g23*Phi222*Phi233*t3 + g13*Phi213*Phi322*t3 + 
     g23*Phi223*Phi322*t3 + g33*Phi233*Phi322*t3); 
    
    t2srcPhi223 = N*(g11*Phi123*Phi201*t0 + g12*Phi123*Phi202*t0 + 
     g13*Phi123*Phi203*t0 + g12*Phi201*Phi223*t0 + 
     g22*Phi202*Phi223*t0 + g23*Phi203*Phi223*t0 + 
     g13*Phi201*Phi323*t0 + g23*Phi202*Phi323*t0 + 
     g33*Phi203*Phi323*t0 + g11*Phi123*Phi211*t1 + 
     g12*Phi123*Phi212*t1 + g13*Phi123*Phi213*t1 + 
     g12*Phi211*Phi223*t1 + g22*Phi212*Phi223*t1 + 
     g23*Phi213*Phi223*t1 + g13*Phi211*Phi323*t1 + 
     g23*Phi212*Phi323*t1 + g33*Phi213*Phi323*t1 + 
     g11*Phi123*Phi212*t2 + g12*Phi123*Phi222*t2 + 
     g13*Phi123*Phi223*t2 + g12*Phi212*Phi223*t2 + 
     g22*Phi222*Phi223*t2 + g23*Power(Phi223,2)*t2 + 
     g13*Phi212*Phi323*t2 + g23*Phi222*Phi323*t2 + 
     g33*Phi223*Phi323*t2 + g11*Phi123*Phi213*t3 + 
     g12*Phi123*Phi223*t3 + g12*Phi213*Phi223*t3 + 
     g22*Power(Phi223,2)*t3 + g13*Phi123*Phi233*t3 + 
     g23*Phi223*Phi233*t3 + g13*Phi213*Phi323*t3 + 
     g23*Phi223*Phi323*t3 + g33*Phi233*Phi323*t3);
    
    t2srcPhi233 = N*(g11*Phi133*Phi201*t0 + g12*Phi133*Phi202*t0 + 
     g13*Phi133*Phi203*t0 + g12*Phi201*Phi233*t0 + 
     g22*Phi202*Phi233*t0 + g23*Phi203*Phi233*t0 + 
     g13*Phi201*Phi333*t0 + g23*Phi202*Phi333*t0 + 
     g33*Phi203*Phi333*t0 + g11*Phi133*Phi211*t1 + 
     g12*Phi133*Phi212*t1 + g13*Phi133*Phi213*t1 + 
     g12*Phi211*Phi233*t1 + g22*Phi212*Phi233*t1 + 
     g23*Phi213*Phi233*t1 + g13*Phi211*Phi333*t1 + 
     g23*Phi212*Phi333*t1 + g33*Phi213*Phi333*t1 + 
     g11*Phi133*Phi212*t2 + g12*Phi133*Phi222*t2 + 
     g13*Phi133*Phi223*t2 + g12*Phi212*Phi233*t2 + 
     g22*Phi222*Phi233*t2 + g23*Phi223*Phi233*t2 + 
     g13*Phi212*Phi333*t2 + g23*Phi222*Phi333*t2 + 
     g33*Phi223*Phi333*t2 + g11*Phi133*Phi213*t3 + 
     g12*Phi133*Phi223*t3 + g13*Phi133*Phi233*t3 + 
     g12*Phi213*Phi233*t3 + g22*Phi223*Phi233*t3 + 
     g23*Power(Phi233,2)*t3 + g13*Phi213*Phi333*t3 + 
     g23*Phi223*Phi333*t3 + g33*Phi233*Phi333*t3); 
    
    t2srcPhi300 = N*(g11*Phi100*Phi301*t0 + g12*Phi200*Phi301*t0 + 
     g13*Phi300*Phi301*t0 + g12*Phi100*Phi302*t0 + 
     g22*Phi200*Phi302*t0 + g23*Phi300*Phi302*t0 + 
     g13*Phi100*Phi303*t0 + g23*Phi200*Phi303*t0 + 
     g33*Phi300*Phi303*t0 + g11*Phi100*Phi311*t1 + 
     g12*Phi200*Phi311*t1 + g13*Phi300*Phi311*t1 + 
     g12*Phi100*Phi312*t1 + g22*Phi200*Phi312*t1 + 
     g23*Phi300*Phi312*t1 + g13*Phi100*Phi313*t1 + 
     g23*Phi200*Phi313*t1 + g33*Phi300*Phi313*t1 + 
     g11*Phi100*Phi312*t2 + g12*Phi200*Phi312*t2 + 
     g13*Phi300*Phi312*t2 + g12*Phi100*Phi322*t2 + 
     g22*Phi200*Phi322*t2 + g23*Phi300*Phi322*t2 + 
     g13*Phi100*Phi323*t2 + g23*Phi200*Phi323*t2 + 
     g33*Phi300*Phi323*t2 + g11*Phi100*Phi313*t3 + 
     g12*Phi200*Phi313*t3 + g13*Phi300*Phi313*t3 + 
     g12*Phi100*Phi323*t3 + g22*Phi200*Phi323*t3 + 
     g23*Phi300*Phi323*t3 + g13*Phi100*Phi333*t3 + 
     g23*Phi200*Phi333*t3 + g33*Phi300*Phi333*t3); 
    
    t2srcPhi301 = N*(g11*Phi101*Phi301*t0 + g12*Phi201*Phi301*t0 + 
     g13*Power(Phi301,2)*t0 + g12*Phi101*Phi302*t0 + 
     g22*Phi201*Phi302*t0 + g23*Phi301*Phi302*t0 + 
     g13*Phi101*Phi303*t0 + g23*Phi201*Phi303*t0 + 
     g33*Phi301*Phi303*t0 + g11*Phi101*Phi311*t1 + 
     g12*Phi201*Phi311*t1 + g13*Phi301*Phi311*t1 + 
     g12*Phi101*Phi312*t1 + g22*Phi201*Phi312*t1 + 
     g23*Phi301*Phi312*t1 + g13*Phi101*Phi313*t1 + 
     g23*Phi201*Phi313*t1 + g33*Phi301*Phi313*t1 + 
     g11*Phi101*Phi312*t2 + g12*Phi201*Phi312*t2 + 
     g13*Phi301*Phi312*t2 + g12*Phi101*Phi322*t2 + 
     g22*Phi201*Phi322*t2 + g23*Phi301*Phi322*t2 + 
     g13*Phi101*Phi323*t2 + g23*Phi201*Phi323*t2 + 
     g33*Phi301*Phi323*t2 + g11*Phi101*Phi313*t3 + 
     g12*Phi201*Phi313*t3 + g13*Phi301*Phi313*t3 + 
     g12*Phi101*Phi323*t3 + g22*Phi201*Phi323*t3 + 
     g23*Phi301*Phi323*t3 + g13*Phi101*Phi333*t3 + 
     g23*Phi201*Phi333*t3 + g33*Phi301*Phi333*t3); 
    
    t2srcPhi302 = N*(g11*Phi102*Phi301*t0 + g12*Phi202*Phi301*t0 + 
     g12*Phi102*Phi302*t0 + g22*Phi202*Phi302*t0 + 
     g13*Phi301*Phi302*t0 + g23*Power(Phi302,2)*t0 + 
     g13*Phi102*Phi303*t0 + g23*Phi202*Phi303*t0 + 
     g33*Phi302*Phi303*t0 + g11*Phi102*Phi311*t1 + 
     g12*Phi202*Phi311*t1 + g13*Phi302*Phi311*t1 + 
     g12*Phi102*Phi312*t1 + g22*Phi202*Phi312*t1 + 
     g23*Phi302*Phi312*t1 + g13*Phi102*Phi313*t1 + 
     g23*Phi202*Phi313*t1 + g33*Phi302*Phi313*t1 + 
     g11*Phi102*Phi312*t2 + g12*Phi202*Phi312*t2 + 
     g13*Phi302*Phi312*t2 + g12*Phi102*Phi322*t2 + 
     g22*Phi202*Phi322*t2 + g23*Phi302*Phi322*t2 + 
     g13*Phi102*Phi323*t2 + g23*Phi202*Phi323*t2 + 
     g33*Phi302*Phi323*t2 + g11*Phi102*Phi313*t3 + 
     g12*Phi202*Phi313*t3 + g13*Phi302*Phi313*t3 + 
     g12*Phi102*Phi323*t3 + g22*Phi202*Phi323*t3 + 
     g23*Phi302*Phi323*t3 + g13*Phi102*Phi333*t3 + 
     g23*Phi202*Phi333*t3 + g33*Phi302*Phi333*t3); 
    
    t2srcPhi303 = N*(g11*Phi103*Phi301*t0 + g12*Phi203*Phi301*t0 + 
     g12*Phi103*Phi302*t0 + g22*Phi203*Phi302*t0 + 
     g13*Phi103*Phi303*t0 + g23*Phi203*Phi303*t0 + 
     g13*Phi301*Phi303*t0 + g23*Phi302*Phi303*t0 + 
     g33*Power(Phi303,2)*t0 + g11*Phi103*Phi311*t1 + 
     g12*Phi203*Phi311*t1 + g13*Phi303*Phi311*t1 + 
     g12*Phi103*Phi312*t1 + g22*Phi203*Phi312*t1 + 
     g23*Phi303*Phi312*t1 + g13*Phi103*Phi313*t1 + 
     g23*Phi203*Phi313*t1 + g33*Phi303*Phi313*t1 + 
     g11*Phi103*Phi312*t2 + g12*Phi203*Phi312*t2 + 
     g13*Phi303*Phi312*t2 + g12*Phi103*Phi322*t2 + 
     g22*Phi203*Phi322*t2 + g23*Phi303*Phi322*t2 + 
     g13*Phi103*Phi323*t2 + g23*Phi203*Phi323*t2 + 
     g33*Phi303*Phi323*t2 + g11*Phi103*Phi313*t3 + 
     g12*Phi203*Phi313*t3 + g13*Phi303*Phi313*t3 + 
     g12*Phi103*Phi323*t3 + g22*Phi203*Phi323*t3 + 
     g23*Phi303*Phi323*t3 + g13*Phi103*Phi333*t3 + 
     g23*Phi203*Phi333*t3 + g33*Phi303*Phi333*t3); 
    
    t2srcPhi311 = N*(g11*Phi111*Phi301*t0 + g12*Phi211*Phi301*t0 + 
     g12*Phi111*Phi302*t0 + g22*Phi211*Phi302*t0 + 
     g13*Phi111*Phi303*t0 + g23*Phi211*Phi303*t0 + 
     g13*Phi301*Phi311*t0 + g23*Phi302*Phi311*t0 + 
     g33*Phi303*Phi311*t0 + g11*Phi111*Phi311*t1 + 
     g12*Phi211*Phi311*t1 + g13*Power(Phi311,2)*t1 + 
     g12*Phi111*Phi312*t1 + g22*Phi211*Phi312*t1 + 
     g23*Phi311*Phi312*t1 + g13*Phi111*Phi313*t1 + 
     g23*Phi211*Phi313*t1 + g33*Phi311*Phi313*t1 + 
     g11*Phi111*Phi312*t2 + g12*Phi211*Phi312*t2 + 
     g13*Phi311*Phi312*t2 + g12*Phi111*Phi322*t2 + 
     g22*Phi211*Phi322*t2 + g23*Phi311*Phi322*t2 + 
     g13*Phi111*Phi323*t2 + g23*Phi211*Phi323*t2 + 
     g33*Phi311*Phi323*t2 + g11*Phi111*Phi313*t3 + 
     g12*Phi211*Phi313*t3 + g13*Phi311*Phi313*t3 + 
     g12*Phi111*Phi323*t3 + g22*Phi211*Phi323*t3 + 
     g23*Phi311*Phi323*t3 + g13*Phi111*Phi333*t3 + 
     g23*Phi211*Phi333*t3 + g33*Phi311*Phi333*t3); 
    
    t2srcPhi312 = N*(g11*Phi112*Phi301*t0 + g12*Phi212*Phi301*t0 + 
     g12*Phi112*Phi302*t0 + g22*Phi212*Phi302*t0 + 
     g13*Phi112*Phi303*t0 + g23*Phi212*Phi303*t0 + 
     g13*Phi301*Phi312*t0 + g23*Phi302*Phi312*t0 + 
     g33*Phi303*Phi312*t0 + g11*Phi112*Phi311*t1 + 
     g12*Phi212*Phi311*t1 + g12*Phi112*Phi312*t1 + 
     g22*Phi212*Phi312*t1 + g13*Phi311*Phi312*t1 + 
     g23*Power(Phi312,2)*t1 + g13*Phi112*Phi313*t1 + 
     g23*Phi212*Phi313*t1 + g33*Phi312*Phi313*t1 + 
     g11*Phi112*Phi312*t2 + g12*Phi212*Phi312*t2 + 
     g13*Power(Phi312,2)*t2 + g12*Phi112*Phi322*t2 + 
     g22*Phi212*Phi322*t2 + g23*Phi312*Phi322*t2 + 
     g13*Phi112*Phi323*t2 + g23*Phi212*Phi323*t2 + 
     g33*Phi312*Phi323*t2 + g11*Phi112*Phi313*t3 + 
     g12*Phi212*Phi313*t3 + g13*Phi312*Phi313*t3 + 
     g12*Phi112*Phi323*t3 + g22*Phi212*Phi323*t3 + 
     g23*Phi312*Phi323*t3 + g13*Phi112*Phi333*t3 + 
     g23*Phi212*Phi333*t3 + g33*Phi312*Phi333*t3); 
    
    t2srcPhi313 = N*(g11*Phi113*Phi301*t0 + g12*Phi213*Phi301*t0 + 
     g12*Phi113*Phi302*t0 + g22*Phi213*Phi302*t0 + 
     g13*Phi113*Phi303*t0 + g23*Phi213*Phi303*t0 + 
     g13*Phi301*Phi313*t0 + g23*Phi302*Phi313*t0 + 
     g33*Phi303*Phi313*t0 + g11*Phi113*Phi311*t1 + 
     g12*Phi213*Phi311*t1 + g12*Phi113*Phi312*t1 + 
     g22*Phi213*Phi312*t1 + g13*Phi113*Phi313*t1 + 
     g23*Phi213*Phi313*t1 + g13*Phi311*Phi313*t1 + 
     g23*Phi312*Phi313*t1 + g33*Power(Phi313,2)*t1 + 
     g11*Phi113*Phi312*t2 + g12*Phi213*Phi312*t2 + 
     g13*Phi312*Phi313*t2 + g12*Phi113*Phi322*t2 + 
     g22*Phi213*Phi322*t2 + g23*Phi313*Phi322*t2 + 
     g13*Phi113*Phi323*t2 + g23*Phi213*Phi323*t2 + 
     g33*Phi313*Phi323*t2 + g11*Phi113*Phi313*t3 + 
     g12*Phi213*Phi313*t3 + g13*Power(Phi313,2)*t3 + 
     g12*Phi113*Phi323*t3 + g22*Phi213*Phi323*t3 + 
     g23*Phi313*Phi323*t3 + g13*Phi113*Phi333*t3 + 
     g23*Phi213*Phi333*t3 + g33*Phi313*Phi333*t3); 
    
    t2srcPhi322 = N*(g11*Phi122*Phi301*t0 + g12*Phi222*Phi301*t0 + 
     g12*Phi122*Phi302*t0 + g22*Phi222*Phi302*t0 + 
     g13*Phi122*Phi303*t0 + g23*Phi222*Phi303*t0 + 
     g13*Phi301*Phi322*t0 + g23*Phi302*Phi322*t0 + 
     g33*Phi303*Phi322*t0 + g11*Phi122*Phi311*t1 + 
     g12*Phi222*Phi311*t1 + g12*Phi122*Phi312*t1 + 
     g22*Phi222*Phi312*t1 + g13*Phi122*Phi313*t1 + 
     g23*Phi222*Phi313*t1 + g13*Phi311*Phi322*t1 + 
     g23*Phi312*Phi322*t1 + g33*Phi313*Phi322*t1 + 
     g11*Phi122*Phi312*t2 + g12*Phi222*Phi312*t2 + 
     g12*Phi122*Phi322*t2 + g22*Phi222*Phi322*t2 + 
     g13*Phi312*Phi322*t2 + g23*Power(Phi322,2)*t2 + 
     g13*Phi122*Phi323*t2 + g23*Phi222*Phi323*t2 + 
     g33*Phi322*Phi323*t2 + g11*Phi122*Phi313*t3 + 
     g12*Phi222*Phi313*t3 + g13*Phi313*Phi322*t3 + 
     g12*Phi122*Phi323*t3 + g22*Phi222*Phi323*t3 + 
     g23*Phi322*Phi323*t3 + g13*Phi122*Phi333*t3 + 
     g23*Phi222*Phi333*t3 + g33*Phi322*Phi333*t3); 
    
    t2srcPhi323 = N*(g11*Phi123*Phi301*t0 + g12*Phi223*Phi301*t0 + 
     g12*Phi123*Phi302*t0 + g22*Phi223*Phi302*t0 + 
     g13*Phi123*Phi303*t0 + g23*Phi223*Phi303*t0 + 
     g13*Phi301*Phi323*t0 + g23*Phi302*Phi323*t0 + 
     g33*Phi303*Phi323*t0 + g11*Phi123*Phi311*t1 + 
     g12*Phi223*Phi311*t1 + g12*Phi123*Phi312*t1 + 
     g22*Phi223*Phi312*t1 + g13*Phi123*Phi313*t1 + 
     g23*Phi223*Phi313*t1 + g13*Phi311*Phi323*t1 + 
     g23*Phi312*Phi323*t1 + g33*Phi313*Phi323*t1 + 
     g11*Phi123*Phi312*t2 + g12*Phi223*Phi312*t2 + 
     g12*Phi123*Phi322*t2 + g22*Phi223*Phi322*t2 + 
     g13*Phi123*Phi323*t2 + g23*Phi223*Phi323*t2 + 
     g13*Phi312*Phi323*t2 + g23*Phi322*Phi323*t2 + 
     g33*Power(Phi323,2)*t2 + g11*Phi123*Phi313*t3 + 
     g12*Phi223*Phi313*t3 + g12*Phi123*Phi323*t3 + 
     g22*Phi223*Phi323*t3 + g13*Phi313*Phi323*t3 + 
     g23*Power(Phi323,2)*t3 + g13*Phi123*Phi333*t3 + 
     g23*Phi223*Phi333*t3 + g33*Phi323*Phi333*t3);
    
    t2srcPhi333 = N*(g11*Phi133*Phi301*t0 + g12*Phi233*Phi301*t0 + 
     g12*Phi133*Phi302*t0 + g22*Phi233*Phi302*t0 + 
     g13*Phi133*Phi303*t0 + g23*Phi233*Phi303*t0 + 
     g13*Phi301*Phi333*t0 + g23*Phi302*Phi333*t0 + 
     g33*Phi303*Phi333*t0 + g11*Phi133*Phi311*t1 + 
     g12*Phi233*Phi311*t1 + g12*Phi133*Phi312*t1 + 
     g22*Phi233*Phi312*t1 + g13*Phi133*Phi313*t1 + 
     g23*Phi233*Phi313*t1 + g13*Phi311*Phi333*t1 + 
     g23*Phi312*Phi333*t1 + g33*Phi313*Phi333*t1 + 
     g11*Phi133*Phi312*t2 + g12*Phi233*Phi312*t2 + 
     g12*Phi133*Phi322*t2 + g22*Phi233*Phi322*t2 + 
     g13*Phi133*Phi323*t2 + g23*Phi233*Phi323*t2 + 
     g13*Phi312*Phi333*t2 + g23*Phi322*Phi333*t2 + 
     g33*Phi323*Phi333*t2 + g11*Phi133*Phi313*t3 + 
     g12*Phi233*Phi313*t3 + g12*Phi133*Phi323*t3 + 
     g22*Phi233*Phi323*t3 + g13*Phi133*Phi333*t3 + 
     g23*Phi233*Phi333*t3 + g13*Phi313*Phi333*t3 + 
     g23*Phi323*Phi333*t3 + g33*Power(Phi333,2)*t3); 
    
    
    //TERM 3 of Phi
    t3srcPhi100 = -(gamma2*N*Phi100);
    t3srcPhi101 = -(gamma2*N*Phi101);
    t3srcPhi102 = -(gamma2*N*Phi102);
    t3srcPhi103 = -(gamma2*N*Phi103);
    t3srcPhi111 = -(gamma2*N*Phi111);
    t3srcPhi112 = -(gamma2*N*Phi112);
    t3srcPhi113 = -(gamma2*N*Phi113);
    t3srcPhi122 = -(gamma2*N*Phi122);
    t3srcPhi123 = -(gamma2*N*Phi123);
    t3srcPhi133 = -(gamma2*N*Phi133);
    
    t3srcPhi200 = -(gamma2*N*Phi200);
    t3srcPhi201 = -(gamma2*N*Phi201);
    t3srcPhi202 = -(gamma2*N*Phi202);
    t3srcPhi203 = -(gamma2*N*Phi203);
    t3srcPhi211 = -(gamma2*N*Phi211);
    t3srcPhi212 = -(gamma2*N*Phi212);
    t3srcPhi213 = -(gamma2*N*Phi213);
    t3srcPhi222 = -(gamma2*N*Phi222);
    t3srcPhi223 = -(gamma2*N*Phi223);
    t3srcPhi233 = -(gamma2*N*Phi233);
    
    t3srcPhi300 = -(gamma2*N*Phi300);
    t3srcPhi301 = -(gamma2*N*Phi301);
    t3srcPhi302 = -(gamma2*N*Phi302);
    t3srcPhi303 = -(gamma2*N*Phi303);
    t3srcPhi311 = -(gamma2*N*Phi311);
    t3srcPhi312 = -(gamma2*N*Phi312);
    t3srcPhi313 = -(gamma2*N*Phi313);
    t3srcPhi322 = -(gamma2*N*Phi322);
    t3srcPhi323 = -(gamma2*N*Phi323);
    t3srcPhi333 = -(gamma2*N*Phi333);
    
    //Sum all terms of RHS Phi
    srcPhi100 = t1srcPhi100 + t2srcPhi100 + t3srcPhi100;
    srcPhi101 = t1srcPhi101 + t2srcPhi101 + t3srcPhi101;
    srcPhi102 = t1srcPhi102 + t2srcPhi102 + t3srcPhi102;
    srcPhi103 = t1srcPhi103 + t2srcPhi103 + t3srcPhi103;
    srcPhi111 = t1srcPhi111 + t2srcPhi111 + t3srcPhi111;
    srcPhi112 = t1srcPhi112 + t2srcPhi112 + t3srcPhi112;
    srcPhi113 = t1srcPhi113 + t2srcPhi113 + t3srcPhi113;
    srcPhi122 = t1srcPhi122 + t2srcPhi122 + t3srcPhi122;
    srcPhi123 = t1srcPhi123 + t2srcPhi123 + t3srcPhi123;
    srcPhi133 = t1srcPhi133 + t2srcPhi133 + t3srcPhi133;
                                                       
    srcPhi200 = t1srcPhi200 + t2srcPhi200 + t3srcPhi200;
    srcPhi201 = t1srcPhi201 + t2srcPhi201 + t3srcPhi201;
    srcPhi202 = t1srcPhi202 + t2srcPhi202 + t3srcPhi202;
    srcPhi203 = t1srcPhi203 + t2srcPhi203 + t3srcPhi203;
    srcPhi211 = t1srcPhi211 + t2srcPhi211 + t3srcPhi211;
    srcPhi212 = t1srcPhi212 + t2srcPhi212 + t3srcPhi212;
    srcPhi213 = t1srcPhi213 + t2srcPhi213 + t3srcPhi213;
    srcPhi222 = t1srcPhi222 + t2srcPhi222 + t3srcPhi222;
    srcPhi223 = t1srcPhi223 + t2srcPhi223 + t3srcPhi223;
    srcPhi233 = t1srcPhi233 + t2srcPhi233 + t3srcPhi233;
                                                       
    srcPhi300 = t1srcPhi300 + t2srcPhi300 + t3srcPhi300;
    srcPhi301 = t1srcPhi301 + t2srcPhi301 + t3srcPhi301;
    srcPhi302 = t1srcPhi302 + t2srcPhi302 + t3srcPhi302;
    srcPhi303 = t1srcPhi303 + t2srcPhi303 + t3srcPhi303;
    srcPhi311 = t1srcPhi311 + t2srcPhi311 + t3srcPhi311;
    srcPhi312 = t1srcPhi312 + t2srcPhi312 + t3srcPhi312;
    srcPhi313 = t1srcPhi313 + t2srcPhi313 + t3srcPhi313;
    srcPhi322 = t1srcPhi322 + t2srcPhi322 + t3srcPhi322;
    srcPhi323 = t1srcPhi323 + t2srcPhi323 + t3srcPhi323;
    srcPhi333 = t1srcPhi333 + t2srcPhi333 + t3srcPhi333;
    
    //put values of src in values[];
    values[0] = srcPsi00; values[1] = srcPsi01; values[2] = srcPsi02; values[3] = srcPsi03; values[4] = srcPsi11;
    values[5] = srcPsi12; values[6] = srcPsi13; values[7] = srcPsi22; values[8] = srcPsi23; values[9] = srcPsi33;
    values[10] = srcPi00; values[11] = srcPi01; values[12] = srcPi02; values[13] = srcPi03; values[14] = srcPi11;
    values[15] = srcPi12; values[16] = srcPi13; values[17] = srcPi22; values[18] = srcPi23; values[19] = srcPi33;
    values[20] = srcPhi100; values[21] = srcPhi101; values[22] = srcPhi102; values[23] = srcPhi103; values[24] = srcPhi111;
    values[25] = srcPhi112; values[26] = srcPhi113; values[27] = srcPhi122; values[28] = srcPhi123; values[29] = srcPhi133;
    values[30] = srcPhi200; values[31] = srcPhi201; values[32] = srcPhi202; values[33] = srcPhi203; values[34] = srcPhi211;
    values[35] = srcPhi212; values[36] = srcPhi213; values[37] = srcPhi222; values[38] = srcPhi223; values[39] = srcPhi233;
    values[40] = srcPhi300; values[41] = srcPhi301; values[42] = srcPhi302; values[43] = srcPhi303; values[44] = srcPhi311;
    values[45] = srcPhi312; values[46] = srcPhi313; values[47] = srcPhi322; values[48] = srcPhi323; values[49] = srcPhi333;    
}
