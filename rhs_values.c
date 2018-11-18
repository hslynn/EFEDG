#include <math.h>

#define gamma_0 1
#define gamma_1 -1
#define gamma_2 1

//RHS of Psi
rhsPsi00 = (-gamma1)*(N1*Phi100 + N2*Phi200 + N3*Phi300) - N*Pi00; 
rhsPsi01 = (-gamma1)*(N1*Phi101 + N2*Phi201 + N3*Phi301) - N*Pi01; 
rhsPsi02 = (-gamma1)*(N1*Phi102 + N2*Phi202 + N3*Phi302) - N*Pi02; 
rhsPsi03 = (-gamma1)*(N1*Phi103 + N2*Phi203 + N3*Phi303) - N*Pi03; 
rhsPsi11 = (-gamma1)*(N1*Phi111 + N2*Phi211 + N3*Phi311) - N*Pi11; 
rhsPsi12 = (-gamma1)*(N1*Phi112 + N2*Phi212 + N3*Phi312) - N*Pi12;
rhsPsi13 = (-gamma1)*(N1*Phi113 + N2*Phi213 + N3*Phi313) - N*Pi13; 
rhsPsi22 = (-gamma1)*(N1*Phi122 + N2*Phi222 + N3*Phi322) - N*Pi22; 
rhsPsi23 = (-gamma1)*(N1*Phi123 + N2*Phi223 + N3*Phi323) - N*Pi23; 
rhsPsi33 = (-gamma1)*(N1*Phi133 + N2*Phi233 + N3*Phi333) - N*Pi33;

//RHS of Pi

//TERM 1
t1rhsPi00 = 2*N*(invPsi00*(-(Power(Gamma000,2)*invPsi00) - 2*Gamma000*Gamma001*invPsi01 - 2*Gamma000*Gamma002*invPsi02 - 
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

t1rhsPi01 = 2*N*(invPsi00*(-(Gamma000*Gamma100*invPsi00) - Gamma001*Gamma100*invPsi01 - Gamma000*Gamma101*invPsi01 - 
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

t1rhsPi02 = 2*N*(invPsi00*(-(Gamma000*Gamma200*invPsi00) - Gamma001*Gamma200*invPsi01 - Gamma000*Gamma201*invPsi01 - 
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

t1rhsPi03 = 2*N*(invPsi00*(-(Gamma000*Gamma300*invPsi00) - Gamma001*Gamma300*invPsi01 - Gamma000*Gamma301*invPsi01 - 
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

t1rhsPi11 = 2*N*(invPsi00*(-(Power(Gamma100,2)*invPsi00) - 2*Gamma100*Gamma101*invPsi01 - 2*Gamma100*Gamma102*invPsi02 - 
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

t1rhsPi12 = 2*N*(invPsi00*(-(Gamma100*Gamma200*invPsi00) - Gamma101*Gamma200*invPsi01 - Gamma100*Gamma201*invPsi01 - 
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

t1rhsPi13 = 2*N*(invPsi00*(-(Gamma100*Gamma300*invPsi00) - Gamma101*Gamma300*invPsi01 - Gamma100*Gamma301*invPsi01 - 
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

t1rhsPi22 = 2*N*(invPsi00*(-(Power(Gamma200,2)*invPsi00) - 2*Gamma200*Gamma201*invPsi01 - 2*Gamma200*Gamma202*invPsi02 - 
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

t1rhsPi23 = 2*N*(invPsi00*(-(Gamma200*Gamma300*invPsi00) - Gamma201*Gamma300*invPsi01 - Gamma200*Gamma301*invPsi01 - 
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

t1rhsPi33 = 2*N*(invPsi00*(-(Power(Gamma300,2)*invPsi00) - 2*Gamma300*Gamma301*invPsi01 - 2*Gamma300*Gamma302*invPsi02 - 
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
t2rhsPi00 = -((2*dH00 - 2*Gamma000*H0*invPsi00 - 2*Gamma100*H0*invPsi01 - 2*Gamma000*H1*invPsi01 - 2*Gamma200*H0*invPsi02 - 
       2*Gamma000*H2*invPsi02 - 2*Gamma300*H0*invPsi03 - 2*Gamma000*H3*invPsi03 - 2*Gamma100*H1*invPsi11 - 
       2*Gamma200*H1*invPsi12 - 2*Gamma100*H2*invPsi12 - 2*Gamma300*H1*invPsi13 - 2*Gamma100*H3*invPsi13 - 
       2*Gamma200*H2*invPsi22 - 2*Gamma300*H2*invPsi23 - 2*Gamma200*H3*invPsi23 - 2*Gamma300*H3*invPsi33)*N);

t2rhsPi01 = -((dH01 + dH10 - 2*Gamma001*H0*invPsi00 - 2*Gamma101*H0*invPsi01 - 2*Gamma001*H1*invPsi01 - 2*Gamma201*H0*invPsi02 - 
       2*Gamma001*H2*invPsi02 - 2*Gamma301*H0*invPsi03 - 2*Gamma001*H3*invPsi03 - 2*Gamma101*H1*invPsi11 - 
       2*Gamma201*H1*invPsi12 - 2*Gamma101*H2*invPsi12 - 2*Gamma301*H1*invPsi13 - 2*Gamma101*H3*invPsi13 - 
       2*Gamma201*H2*invPsi22 - 2*Gamma301*H2*invPsi23 - 2*Gamma201*H3*invPsi23 - 2*Gamma301*H3*invPsi33)*N);

t2rhsPi02 = -((dH02 + dH20 - 2*Gamma002*H0*invPsi00 - 2*Gamma102*H0*invPsi01 - 2*Gamma002*H1*invPsi01 - 2*Gamma202*H0*invPsi02 - 
       2*Gamma002*H2*invPsi02 - 2*Gamma302*H0*invPsi03 - 2*Gamma002*H3*invPsi03 - 2*Gamma102*H1*invPsi11 - 
       2*Gamma202*H1*invPsi12 - 2*Gamma102*H2*invPsi12 - 2*Gamma302*H1*invPsi13 - 2*Gamma102*H3*invPsi13 - 
       2*Gamma202*H2*invPsi22 - 2*Gamma302*H2*invPsi23 - 2*Gamma202*H3*invPsi23 - 2*Gamma302*H3*invPsi33)*N);

t2rhsPi03 = -((dH03 + dH30 - 2*Gamma003*H0*invPsi00 - 2*Gamma103*H0*invPsi01 - 2*Gamma003*H1*invPsi01 - 2*Gamma203*H0*invPsi02 - 
       2*Gamma003*H2*invPsi02 - 2*Gamma303*H0*invPsi03 - 2*Gamma003*H3*invPsi03 - 2*Gamma103*H1*invPsi11 - 
       2*Gamma203*H1*invPsi12 - 2*Gamma103*H2*invPsi12 - 2*Gamma303*H1*invPsi13 - 2*Gamma103*H3*invPsi13 - 
       2*Gamma203*H2*invPsi22 - 2*Gamma303*H2*invPsi23 - 2*Gamma203*H3*invPsi23 - 2*Gamma303*H3*invPsi33)*N);

t2rhsPi11 = -((2*dH11 - 2*Gamma011*H0*invPsi00 - 2*Gamma111*H0*invPsi01 - 2*Gamma011*H1*invPsi01 - 2*Gamma211*H0*invPsi02 - 
       2*Gamma011*H2*invPsi02 - 2*Gamma311*H0*invPsi03 - 2*Gamma011*H3*invPsi03 - 2*Gamma111*H1*invPsi11 - 
       2*Gamma211*H1*invPsi12 - 2*Gamma111*H2*invPsi12 - 2*Gamma311*H1*invPsi13 - 2*Gamma111*H3*invPsi13 - 
       2*Gamma211*H2*invPsi22 - 2*Gamma311*H2*invPsi23 - 2*Gamma211*H3*invPsi23 - 2*Gamma311*H3*invPsi33)*N);

t2rhsPi12 = -((dH12 + dH21 - 2*Gamma012*H0*invPsi00 - 2*Gamma112*H0*invPsi01 - 2*Gamma012*H1*invPsi01 - 2*Gamma212*H0*invPsi02 - 
       2*Gamma012*H2*invPsi02 - 2*Gamma312*H0*invPsi03 - 2*Gamma012*H3*invPsi03 - 2*Gamma112*H1*invPsi11 - 
       2*Gamma212*H1*invPsi12 - 2*Gamma112*H2*invPsi12 - 2*Gamma312*H1*invPsi13 - 2*Gamma112*H3*invPsi13 - 
       2*Gamma212*H2*invPsi22 - 2*Gamma312*H2*invPsi23 - 2*Gamma212*H3*invPsi23 - 2*Gamma312*H3*invPsi33)*N);

t2rhsPi13 = -((dH13 + dH31 - 2*Gamma013*H0*invPsi00 - 2*Gamma113*H0*invPsi01 - 2*Gamma013*H1*invPsi01 - 2*Gamma213*H0*invPsi02 - 
       2*Gamma013*H2*invPsi02 - 2*Gamma313*H0*invPsi03 - 2*Gamma013*H3*invPsi03 - 2*Gamma113*H1*invPsi11 - 
       2*Gamma213*H1*invPsi12 - 2*Gamma113*H2*invPsi12 - 2*Gamma313*H1*invPsi13 - 2*Gamma113*H3*invPsi13 - 
       2*Gamma213*H2*invPsi22 - 2*Gamma313*H2*invPsi23 - 2*Gamma213*H3*invPsi23 - 2*Gamma313*H3*invPsi33)*N);

t2rhsPi22 = -((2*dH22 - 2*Gamma022*H0*invPsi00 - 2*Gamma122*H0*invPsi01 - 2*Gamma022*H1*invPsi01 - 2*Gamma222*H0*invPsi02 - 
       2*Gamma022*H2*invPsi02 - 2*Gamma322*H0*invPsi03 - 2*Gamma022*H3*invPsi03 - 2*Gamma122*H1*invPsi11 - 
       2*Gamma222*H1*invPsi12 - 2*Gamma122*H2*invPsi12 - 2*Gamma322*H1*invPsi13 - 2*Gamma122*H3*invPsi13 - 
       2*Gamma222*H2*invPsi22 - 2*Gamma322*H2*invPsi23 - 2*Gamma222*H3*invPsi23 - 2*Gamma322*H3*invPsi33)*N);

t2rhsPi23 = -((dH23 + dH32 - 2*Gamma023*H0*invPsi00 - 2*Gamma123*H0*invPsi01 - 2*Gamma023*H1*invPsi01 - 2*Gamma223*H0*invPsi02 - 
       2*Gamma023*H2*invPsi02 - 2*Gamma323*H0*invPsi03 - 2*Gamma023*H3*invPsi03 - 2*Gamma123*H1*invPsi11 - 
       2*Gamma223*H1*invPsi12 - 2*Gamma123*H2*invPsi12 - 2*Gamma323*H1*invPsi13 - 2*Gamma123*H3*invPsi13 - 
       2*Gamma223*H2*invPsi22 - 2*Gamma323*H2*invPsi23 - 2*Gamma223*H3*invPsi23 - 2*Gamma323*H3*invPsi33)*N);

t2rhsPi33 = -((2*dH33 - 2*Gamma033*H0*invPsi00 - 2*Gamma133*H0*invPsi01 - 2*Gamma033*H1*invPsi01 - 2*Gamma233*H0*invPsi02 - 
       2*Gamma033*H2*invPsi02 - 2*Gamma333*H0*invPsi03 - 2*Gamma033*H3*invPsi03 - 2*Gamma133*H1*invPsi11 - 
       2*Gamma233*H1*invPsi12 - 2*Gamma133*H2*invPsi12 - 2*Gamma333*H1*invPsi13 - 2*Gamma133*H3*invPsi13 - 
       2*Gamma233*H2*invPsi22 - 2*Gamma333*H2*invPsi23 - 2*Gamma233*H3*invPsi23 - 2*Gamma333*H3*invPsi33)*N);

//TERM 3

t3rhsPi00 = -0.5*N*(Power(Pi00,2)*Power(t0,2) + 2*Pi00*Pi01*t0*t1 + Pi00*Pi11*Power(t1,2) + 2*Pi00*Pi02*t0*t2 + 
       2*Pi00*Pi12*t1*t2 + Pi00*Pi22*Power(t2,2) + 2*Pi00*Pi03*t0*t3 + 2*Pi00*Pi13*t1*t3 + 2*Pi00*Pi23*t2*t3 + 
       Pi00*Pi33*Power(t3,2)); 

t3rhsPi01 = -0.5*N*(Pi00*Pi01*Power(t0,2) + 2*Power(Pi01,2)*t0*t1 + Pi01*Pi11*Power(t1,2) + 
       2*Pi01*Pi02*t0*t2 + 2*Pi01*Pi12*t1*t2 + Pi01*Pi22*Power(t2,2) + 2*Pi01*Pi03*t0*t3 + 2*Pi01*Pi13*t1*t3 + 
       2*Pi01*Pi23*t2*t3 + Pi01*Pi33*Power(t3,2));

t3rhsPi02 = -0.5*N*(Pi00*Pi02*Power(t0,2) + 2*Pi01*Pi02*t0*t1 + Pi02*Pi11*Power(t1,2) + 2*Power(Pi02,2)*t0*t2 + 2*Pi02*Pi12*t1*t2 + 
       Pi02*Pi22*Power(t2,2) + 2*Pi02*Pi03*t0*t3 + 2*Pi02*Pi13*t1*t3 + 2*Pi02*Pi23*t2*t3 + Pi02*Pi33*Power(t3,2));

t3rhsPi03 = -0.5*N*(Pi00*Pi03*Power(t0,2) + 2*Pi01*Pi03*t0*t1 + Pi03*Pi11*Power(t1,2) + 2*Pi02*Pi03*t0*t2 + 2*Pi03*Pi12*t1*t2 + 
       Pi03*Pi22*Power(t2,2) + 2*Power(Pi03,2)*t0*t3 + 2*Pi03*Pi13*t1*t3 + 2*Pi03*Pi23*t2*t3 + Pi03*Pi33*Power(t3,2));

t3rhsPi11 = -0.5*N*(Pi00*Pi11*Power(t0,2) + 2*Pi01*Pi11*t0*t1 + Power(Pi11,2)*Power(t1,2) + 2*Pi02*Pi11*t0*t2 + 2*Pi11*Pi12*t1*t2 + 
       Pi11*Pi22*Power(t2,2) + 2*Pi03*Pi11*t0*t3 + 2*Pi11*Pi13*t1*t3 + 2*Pi11*Pi23*t2*t3 + Pi11*Pi33*Power(t3,2));

t3rhsPi12 = -0.5*N*(Pi00*Pi12*Power(t0,2) + 2*Pi01*Pi12*t0*t1 + Pi11*Pi12*Power(t1,2) + 2*Pi02*Pi12*t0*t2 + 2*Power(Pi12,2)*t1*t2 + 
       Pi12*Pi22*Power(t2,2) + 2*Pi03*Pi12*t0*t3 + 2*Pi12*Pi13*t1*t3 + 2*Pi12*Pi23*t2*t3 + Pi12*Pi33*Power(t3,2));

t3rhsPi13 = -0.5*N*(Pi00*Pi13*Power(t0,2) + 2*Pi01*Pi13*t0*t1 + Pi11*Pi13*Power(t1,2) + 2*Pi02*Pi13*t0*t2 + 2*Pi12*Pi13*t1*t2 + 
       Pi13*Pi22*Power(t2,2) + 2*Pi03*Pi13*t0*t3 + 2*Power(Pi13,2)*t1*t3 + 2*Pi13*Pi23*t2*t3 + Pi13*Pi33*Power(t3,2));

t3rhsPi22 = -0.5*N*(Pi00*Pi22*Power(t0,2) + 2*Pi01*Pi22*t0*t1 + Pi11*Pi22*Power(t1,2) + 2*Pi02*Pi22*t0*t2 + 2*Pi12*Pi22*t1*t2 + 
       Power(Pi22,2)*Power(t2,2) + 2*Pi03*Pi22*t0*t3 + 2*Pi13*Pi22*t1*t3 + 2*Pi22*Pi23*t2*t3 + Pi22*Pi33*Power(t3,2));

t3rhsPi23 = -0.5*N*(Pi00*Pi23*Power(t0,2) + 2*Pi01*Pi23*t0*t1 + Pi11*Pi23*Power(t1,2) + 2*Pi02*Pi23*t0*t2 + 2*Pi12*Pi23*t1*t2 + 
       Pi22*Pi23*Power(t2,2) + 2*Pi03*Pi23*t0*t3 + 2*Pi13*Pi23*t1*t3 + 2*Power(Pi23,2)*t2*t3 + Pi23*Pi33*Power(t3,2));

t3rhsPi33 = -0.5*N*(Pi00*Pi33*Power(t0,2) + 2*Pi01*Pi33*t0*t1 + Pi11*Pi33*Power(t1,2) + 2*Pi02*Pi33*t0*t2 + 2*Pi12*Pi33*t1*t2 + 
       Pi22*Pi33*Power(t2,2) + 2*Pi03*Pi33*t0*t3 + 2*Pi13*Pi33*t1*t3 + 2*Pi23*Pi33*t2*t3 + Power(Pi33,2)*Power(t3,2));

//TERM 4
t4rhsPi00 = -(N*(g11*Phi100*Pi00*t0 + g12*Phi200*Pi00*t0 + g13*Phi300*Pi00*t0 + g12*Phi100*Pi01*t0 + g22*Phi200*Pi01*t0 + 
         g23*Phi300*Pi01*t0 + g13*Phi100*Pi02*t0 + g23*Phi200*Pi02*t0 + g33*Phi300*Pi02*t0 + g11*Phi100*Pi01*t1 + 
         g12*Phi200*Pi01*t1 + g13*Phi300*Pi01*t1 + g12*Phi100*Pi11*t1 + g22*Phi200*Pi11*t1 + g23*Phi300*Pi11*t1 + 
         g13*Phi100*Pi12*t1 + g23*Phi200*Pi12*t1 + g33*Phi300*Pi12*t1 + g11*Phi100*Pi02*t2 + g12*Phi200*Pi02*t2 + 
         g13*Phi300*Pi02*t2 + g12*Phi100*Pi12*t2 + g22*Phi200*Pi12*t2 + g23*Phi300*Pi12*t2 + g13*Phi100*Pi22*t2 + 
         g23*Phi200*Pi22*t2 + g33*Phi300*Pi22*t2 + g11*Phi100*Pi03*t3 + g12*Phi200*Pi03*t3 + g13*Phi300*Pi03*t3 + 
         g12*Phi100*Pi13*t3 + g22*Phi200*Pi13*t3 + g23*Phi300*Pi13*t3 + g13*Phi100*Pi23*t3 + g23*Phi200*Pi23*t3 + 
         g33*Phi300*Pi23*t3));

t4rhsPi01 = -(N*(g11*Phi101*Pi00*t0 + g12*Phi201*Pi00*t0 + g13*Phi301*Pi00*t0 + g12*Phi101*Pi01*t0 + 
         g22*Phi201*Pi01*t0 + g23*Phi301*Pi01*t0 + g13*Phi101*Pi02*t0 + g23*Phi201*Pi02*t0 + g33*Phi301*Pi02*t0 + 
         g11*Phi101*Pi01*t1 + g12*Phi201*Pi01*t1 + g13*Phi301*Pi01*t1 + g12*Phi101*Pi11*t1 + g22*Phi201*Pi11*t1 + 
         g23*Phi301*Pi11*t1 + g13*Phi101*Pi12*t1 + g23*Phi201*Pi12*t1 + g33*Phi301*Pi12*t1 + g11*Phi101*Pi02*t2 + 
         g12*Phi201*Pi02*t2 + g13*Phi301*Pi02*t2 + g12*Phi101*Pi12*t2 + g22*Phi201*Pi12*t2 + g23*Phi301*Pi12*t2 + 
         g13*Phi101*Pi22*t2 + g23*Phi201*Pi22*t2 + g33*Phi301*Pi22*t2 + g11*Phi101*Pi03*t3 + g12*Phi201*Pi03*t3 + 
         g13*Phi301*Pi03*t3 + g12*Phi101*Pi13*t3 + g22*Phi201*Pi13*t3 + g23*Phi301*Pi13*t3 + g13*Phi101*Pi23*t3 + 
         g23*Phi201*Pi23*t3 + g33*Phi301*Pi23*t3));
            
t4rhsPi02 = -(N*(g11*Phi102*Pi00*t0 + g12*Phi202*Pi00*t0 + g13*Phi302*Pi00*t0 + g12*Phi102*Pi01*t0 + g22*Phi202*Pi01*t0 + 
         g23*Phi302*Pi01*t0 + g13*Phi102*Pi02*t0 + g23*Phi202*Pi02*t0 + g33*Phi302*Pi02*t0 + g11*Phi102*Pi01*t1 + 
         g12*Phi202*Pi01*t1 + g13*Phi302*Pi01*t1 + g12*Phi102*Pi11*t1 + g22*Phi202*Pi11*t1 + g23*Phi302*Pi11*t1 + 
         g13*Phi102*Pi12*t1 + g23*Phi202*Pi12*t1 + g33*Phi302*Pi12*t1 + g11*Phi102*Pi02*t2 + g12*Phi202*Pi02*t2 + 
         g13*Phi302*Pi02*t2 + g12*Phi102*Pi12*t2 + g22*Phi202*Pi12*t2 + g23*Phi302*Pi12*t2 + g13*Phi102*Pi22*t2 + 
         g23*Phi202*Pi22*t2 + g33*Phi302*Pi22*t2 + g11*Phi102*Pi03*t3 + g12*Phi202*Pi03*t3 + g13*Phi302*Pi03*t3 + 
         g12*Phi102*Pi13*t3 + g22*Phi202*Pi13*t3 + g23*Phi302*Pi13*t3 + g13*Phi102*Pi23*t3 + g23*Phi202*Pi23*t3 + 
         g33*Phi302*Pi23*t3));
            
t4rhsPi03 = -(N*(g11*Phi103*Pi00*t0 + g12*Phi203*Pi00*t0 + g13*Phi303*Pi00*t0 + g12*Phi103*Pi01*t0 + 
         g22*Phi203*Pi01*t0 + g23*Phi303*Pi01*t0 + g13*Phi103*Pi02*t0 + g23*Phi203*Pi02*t0 + g33*Phi303*Pi02*t0 + 
         g11*Phi103*Pi01*t1 + g12*Phi203*Pi01*t1 + g13*Phi303*Pi01*t1 + g12*Phi103*Pi11*t1 + g22*Phi203*Pi11*t1 + 
         g23*Phi303*Pi11*t1 + g13*Phi103*Pi12*t1 + g23*Phi203*Pi12*t1 + g33*Phi303*Pi12*t1 + g11*Phi103*Pi02*t2 + 
         g12*Phi203*Pi02*t2 + g13*Phi303*Pi02*t2 + g12*Phi103*Pi12*t2 + g22*Phi203*Pi12*t2 + g23*Phi303*Pi12*t2 + 
         g13*Phi103*Pi22*t2 + g23*Phi203*Pi22*t2 + g33*Phi303*Pi22*t2 + g11*Phi103*Pi03*t3 + g12*Phi203*Pi03*t3 + 
         g13*Phi303*Pi03*t3 + g12*Phi103*Pi13*t3 + g22*Phi203*Pi13*t3 + g23*Phi303*Pi13*t3 + g13*Phi103*Pi23*t3 + 
         g23*Phi203*Pi23*t3 + g33*Phi303*Pi23*t3));

t4rhsPi11 = -(N*(g11*Phi111*Pi00*t0 + g12*Phi211*Pi00*t0 + g13*Phi311*Pi00*t0 + g12*Phi111*Pi01*t0 + 
         g22*Phi211*Pi01*t0 + g23*Phi311*Pi01*t0 + g13*Phi111*Pi02*t0 + g23*Phi211*Pi02*t0 + g33*Phi311*Pi02*t0 + 
         g11*Phi111*Pi01*t1 + g12*Phi211*Pi01*t1 + g13*Phi311*Pi01*t1 + g12*Phi111*Pi11*t1 + g22*Phi211*Pi11*t1 + 
         g23*Phi311*Pi11*t1 + g13*Phi111*Pi12*t1 + g23*Phi211*Pi12*t1 + g33*Phi311*Pi12*t1 + g11*Phi111*Pi02*t2 + 
         g12*Phi211*Pi02*t2 + g13*Phi311*Pi02*t2 + g12*Phi111*Pi12*t2 + g22*Phi211*Pi12*t2 + g23*Phi311*Pi12*t2 + 
         g13*Phi111*Pi22*t2 + g23*Phi211*Pi22*t2 + g33*Phi311*Pi22*t2 + g11*Phi111*Pi03*t3 + g12*Phi211*Pi03*t3 + 
         g13*Phi311*Pi03*t3 + g12*Phi111*Pi13*t3 + g22*Phi211*Pi13*t3 + g23*Phi311*Pi13*t3 + g13*Phi111*Pi23*t3 + 
         g23*Phi211*Pi23*t3 + g33*Phi311*Pi23*t3));

t4rhsPi12 = -(N*(g11*Phi112*Pi00*t0 + g12*Phi212*Pi00*t0 + g13*Phi312*Pi00*t0 + g12*Phi112*Pi01*t0 + g22*Phi212*Pi01*t0 + 
         g23*Phi312*Pi01*t0 + g13*Phi112*Pi02*t0 + g23*Phi212*Pi02*t0 + g33*Phi312*Pi02*t0 + g11*Phi112*Pi01*t1 + 
         g12*Phi212*Pi01*t1 + g13*Phi312*Pi01*t1 + g12*Phi112*Pi11*t1 + g22*Phi212*Pi11*t1 + g23*Phi312*Pi11*t1 + 
         g13*Phi112*Pi12*t1 + g23*Phi212*Pi12*t1 + g33*Phi312*Pi12*t1 + g11*Phi112*Pi02*t2 + g12*Phi212*Pi02*t2 + 
         g13*Phi312*Pi02*t2 + g12*Phi112*Pi12*t2 + g22*Phi212*Pi12*t2 + g23*Phi312*Pi12*t2 + g13*Phi112*Pi22*t2 + 
         g23*Phi212*Pi22*t2 + g33*Phi312*Pi22*t2 + g11*Phi112*Pi03*t3 + g12*Phi212*Pi03*t3 + g13*Phi312*Pi03*t3 + 
         g12*Phi112*Pi13*t3 + g22*Phi212*Pi13*t3 + g23*Phi312*Pi13*t3 + g13*Phi112*Pi23*t3 + g23*Phi212*Pi23*t3 + 
         g33*Phi312*Pi23*t3)); 

t4rhsPi13 = -(N*(g11*Phi113*Pi00*t0 + g12*Phi213*Pi00*t0 + g13*Phi313*Pi00*t0 + g12*Phi113*Pi01*t0 + 
         g22*Phi213*Pi01*t0 + g23*Phi313*Pi01*t0 + g13*Phi113*Pi02*t0 + g23*Phi213*Pi02*t0 + g33*Phi313*Pi02*t0 + 
         g11*Phi113*Pi01*t1 + g12*Phi213*Pi01*t1 + g13*Phi313*Pi01*t1 + g12*Phi113*Pi11*t1 + g22*Phi213*Pi11*t1 + 
         g23*Phi313*Pi11*t1 + g13*Phi113*Pi12*t1 + g23*Phi213*Pi12*t1 + g33*Phi313*Pi12*t1 + g11*Phi113*Pi02*t2 + 
         g12*Phi213*Pi02*t2 + g13*Phi313*Pi02*t2 + g12*Phi113*Pi12*t2 + g22*Phi213*Pi12*t2 + g23*Phi313*Pi12*t2 + 
         g13*Phi113*Pi22*t2 + g23*Phi213*Pi22*t2 + g33*Phi313*Pi22*t2 + g11*Phi113*Pi03*t3 + g12*Phi213*Pi03*t3 + 
         g13*Phi313*Pi03*t3 + g12*Phi113*Pi13*t3 + g22*Phi213*Pi13*t3 + g23*Phi313*Pi13*t3 + g13*Phi113*Pi23*t3 + 
         g23*Phi213*Pi23*t3 + g33*Phi313*Pi23*t3));

t4rhsPi22 = -(N*(g11*Phi122*Pi00*t0 + g12*Phi222*Pi00*t0 + g13*Phi322*Pi00*t0 + g12*Phi122*Pi01*t0 + g22*Phi222*Pi01*t0 + 
         g23*Phi322*Pi01*t0 + g13*Phi122*Pi02*t0 + g23*Phi222*Pi02*t0 + g33*Phi322*Pi02*t0 + g11*Phi122*Pi01*t1 + 
         g12*Phi222*Pi01*t1 + g13*Phi322*Pi01*t1 + g12*Phi122*Pi11*t1 + g22*Phi222*Pi11*t1 + g23*Phi322*Pi11*t1 + 
         g13*Phi122*Pi12*t1 + g23*Phi222*Pi12*t1 + g33*Phi322*Pi12*t1 + g11*Phi122*Pi02*t2 + g12*Phi222*Pi02*t2 + 
         g13*Phi322*Pi02*t2 + g12*Phi122*Pi12*t2 + g22*Phi222*Pi12*t2 + g23*Phi322*Pi12*t2 + g13*Phi122*Pi22*t2 + 
         g23*Phi222*Pi22*t2 + g33*Phi322*Pi22*t2 + g11*Phi122*Pi03*t3 + g12*Phi222*Pi03*t3 + g13*Phi322*Pi03*t3 + 
         g12*Phi122*Pi13*t3 + g22*Phi222*Pi13*t3 + g23*Phi322*Pi13*t3 + g13*Phi122*Pi23*t3 + g23*Phi222*Pi23*t3 + 
         g33*Phi322*Pi23*t3));

t4rhsPi23 = -(N*(g11*Phi123*Pi00*t0 + g12*Phi223*Pi00*t0 + g13*Phi323*Pi00*t0 + g12*Phi123*Pi01*t0 + 
         g22*Phi223*Pi01*t0 + g23*Phi323*Pi01*t0 + g13*Phi123*Pi02*t0 + g23*Phi223*Pi02*t0 + g33*Phi323*Pi02*t0 + 
         g11*Phi123*Pi01*t1 + g12*Phi223*Pi01*t1 + g13*Phi323*Pi01*t1 + g12*Phi123*Pi11*t1 + g22*Phi223*Pi11*t1 + 
         g23*Phi323*Pi11*t1 + g13*Phi123*Pi12*t1 + g23*Phi223*Pi12*t1 + g33*Phi323*Pi12*t1 + g11*Phi123*Pi02*t2 + 
         g12*Phi223*Pi02*t2 + g13*Phi323*Pi02*t2 + g12*Phi123*Pi12*t2 + g22*Phi223*Pi12*t2 + g23*Phi323*Pi12*t2 + 
         g13*Phi123*Pi22*t2 + g23*Phi223*Pi22*t2 + g33*Phi323*Pi22*t2 + g11*Phi123*Pi03*t3 + g12*Phi223*Pi03*t3 + 
         g13*Phi323*Pi03*t3 + g12*Phi123*Pi13*t3 + g22*Phi223*Pi13*t3 + g23*Phi323*Pi13*t3 + g13*Phi123*Pi23*t3 + 
         g23*Phi223*Pi23*t3 + g33*Phi323*Pi23*t3));

t4rhsPi33 = -(N*(g11*Phi133*Pi00*t0 + g12*Phi233*Pi00*t0 + g13*Phi333*Pi00*t0 + g12*Phi133*Pi01*t0 + 
         g22*Phi233*Pi01*t0 + g23*Phi333*Pi01*t0 + g13*Phi133*Pi02*t0 + g23*Phi233*Pi02*t0 + g33*Phi333*Pi02*t0 + 
         g11*Phi133*Pi01*t1 + g12*Phi233*Pi01*t1 + g13*Phi333*Pi01*t1 + g12*Phi133*Pi11*t1 + g22*Phi233*Pi11*t1 + 
         g23*Phi333*Pi11*t1 + g13*Phi133*Pi12*t1 + g23*Phi233*Pi12*t1 + g33*Phi333*Pi12*t1 + g11*Phi133*Pi02*t2 + 
         g12*Phi233*Pi02*t2 + g13*Phi333*Pi02*t2 + g12*Phi133*Pi12*t2 + g22*Phi233*Pi12*t2 + g23*Phi333*Pi12*t2 + 
         g13*Phi133*Pi22*t2 + g23*Phi233*Pi22*t2 + g33*Phi333*Pi22*t2 + g11*Phi133*Pi03*t3 + g12*Phi233*Pi03*t3 + 
         g13*Phi333*Pi03*t3 + g12*Phi133*Pi13*t3 + g22*Phi233*Pi13*t3 + g23*Phi333*Pi13*t3 + g13*Phi133*Pi23*t3 + 
         g23*Phi233*Pi23*t3 + g33*Phi333*Pi23*t3));

//TERM 5

t5rhsPi00 = gamma0*N*((2*t0 - Psi00*t0)*(H0 + vecGamma0) - Psi00*t1*(H1 + vecGamma1) - Psi00*t2*(H2 + vecGamma2) - 
       Psi00*t3*(H3 + vecGamma3));

t5rhsPi01 = gamma0*N*((-(Psi01*t0) + t1)*(H0 + vecGamma0) + (t0 - Psi01*t1)*(H1 + vecGamma1) - 
       Psi01*t2*(H2 + vecGamma2) - Psi01*t3*(H3 + vecGamma3));

t5rhsPi02 = gamma0*N*((-(Psi02*t0) + t2)*(H0 + vecGamma0) - Psi02*t1*(H1 + vecGamma1) + (t0 - Psi02*t2)*(H2 + vecGamma2) - 
       Psi02*t3*(H3 + vecGamma3));

t5rhsPi03 = gamma0*N*((-(Psi03*t0) + t3)*(H0 + vecGamma0) - Psi03*t1*(H1 + vecGamma1) - 
       Psi03*t2*(H2 + vecGamma2) + (t0 - Psi03*t3)*(H3 + vecGamma3));


t5rhsPi11 = gamma0*N*(-(Psi11*t0*(H0 + vecGamma0)) + (2*t1 - Psi11*t1)*(H1 + vecGamma1) - 
       Psi11*t2*(H2 + vecGamma2) - Psi11*t3*(H3 + vecGamma3));
 

t5rhsPi12 = gamma0*N*(-(Psi12*t0*(H0 + vecGamma0)) + (-(Psi12*t1) + t2)*(H1 + vecGamma1) + (t1 - Psi12*t2)*(H2 + vecGamma2) - 
       Psi12*t3*(H3 + vecGamma3));


t5rhsPi13 = gamma0*N*(-(Psi13*t0*(H0 + vecGamma0)) + (-(Psi13*t1) + t3)*(H1 + vecGamma1) - 
       Psi13*t2*(H2 + vecGamma2) + (t1 - Psi13*t3)*(H3 + vecGamma3));

t5rhsPi22 = gamma0*N*(-(Psi22*t0*(H0 + vecGamma0)) - Psi22*t1*(H1 + vecGamma1) + (2*t2 - Psi22*t2)*(H2 + vecGamma2) - 
       Psi22*t3*(H3 + vecGamma3));

t5rhsPi23 = gamma0*N*(-(Psi23*t0*(H0 + vecGamma0)) - Psi23*t1*(H1 + vecGamma1) + 
       (-(Psi23*t2) + t3)*(H2 + vecGamma2) + (t2 - Psi23*t3)*(H3 + vecGamma3));
 
t5rhsPi33 = gamma0*N*(-(Psi33*t0*(H0 + vecGamma0)) - Psi33*t1*(H1 + vecGamma1) - 
       Psi33*t2*(H2 + vecGamma2) + (2*t3 - Psi33*t3)*(H3 + vecGamma3));

//TERM 6
t6rhsPi00 = -(gamma1*gamma2*(N1*Phi100 + N2*Phi200 + N3*Phi300));

t6rhsPi01 = -(gamma1*gamma2*(N1*Phi101 + N2*Phi201 + N3*Phi301));

t6rhsPi02 = -(gamma1*gamma2*(N1*Phi102 + N2*Phi202 + N3*Phi302));

t6rhsPi03 = -(gamma1*gamma2*(N1*Phi103 + N2*Phi203 + N3*Phi303));

t6rhsPi11 = -(gamma1*gamma2*(N1*Phi111 + N2*Phi211 + N3*Phi311));

t6rhsPi12 = -(gamma1*gamma2*(N1*Phi112 + N2*Phi212 + N3*Phi312));

t6rhsPi13 = -(gamma1*gamma2*(N1*Phi113 + N2*Phi213 + N3*Phi313));

t6rhsPi22 = -(gamma1*gamma2*(N1*Phi122 + N2*Phi222 + N3*Phi322));

t6rhsPi23 = -(gamma1*gamma2*(N1*Phi123 + N2*Phi223 + N3*Phi323));

t6rhsPi33 = -(gamma1*gamma2*(N1*Phi133 + N2*Phi233 + N3*Phi333));

//Sum all terms of RHS Pi
rhsPi00 = t1rhsPi00 + t2rhsPi00 + t3rhsPi00 + t4rhsPi00 + t5rhsPi00 + t6rhsPi00;
rhsPi01 = t1rhsPi01 + t2rhsPi01 + t3rhsPi01 + t4rhsPi01 + t5rhsPi01 + t6rhsPi01;
rhsPi02 = t1rhsPi02 + t2rhsPi02 + t3rhsPi02 + t4rhsPi02 + t5rhsPi02 + t6rhsPi02;
rhsPi03 = t1rhsPi03 + t2rhsPi03 + t3rhsPi03 + t4rhsPi03 + t5rhsPi03 + t6rhsPi03;
rhsPi11 = t1rhsPi11 + t2rhsPi11 + t3rhsPi11 + t4rhsPi11 + t5rhsPi11 + t6rhsPi11;
rhsPi12 = t1rhsPi12 + t2rhsPi12 + t3rhsPi12 + t4rhsPi12 + t5rhsPi12 + t6rhsPi12;
rhsPi13 = t1rhsPi13 + t2rhsPi13 + t3rhsPi13 + t4rhsPi13 + t5rhsPi13 + t6rhsPi13;
rhsPi22 = t1rhsPi22 + t2rhsPi22 + t3rhsPi22 + t4rhsPi22 + t5rhsPi22 + t6rhsPi22;
rhsPi23 = t1rhsPi23 + t2rhsPi23 + t3rhsPi23 + t4rhsPi23 + t5rhsPi23 + t6rhsPi23;
rhsPi33 = t1rhsPi33 + t2rhsPi33 + t3rhsPi33 + t4rhsPi33 + t5rhsPi33 + t6rhsPi33;

//RHS of Phi

//TEMR 1 of Phi

t1rhsPhi100 = 0.5*N*Pi00*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi101 =  0.5*N*Pi01*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi102 = 0.5*N*Pi02*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi103 = 0.5*N*Pi03*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi111 = 0.5*N*Pi11*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi112 = 0.5*N*Pi12*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));
    
t1rhsPhi113 = 0.5*N*Pi13*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi122 = 0.5*N*Pi22*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi123 = 0.5*N*Pi23*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi133 = 0.5*N*Pi33*(Phi100*Power(t0,2) + 2*Phi101*t0*t1 + Phi111*Power(t1,2) + 2*Phi102*t0*t2 + 2*Phi112*t1*t2 + 
        Phi122*Power(t2,2) + 2*Phi103*t0*t3 + 2*Phi113*t1*t3 + 2*Phi123*t2*t3 + Phi133*Power(t3,2));

t1rhsPhi200 = 0.5*N*Pi00*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi201 = 0.5*N*Pi01*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi202 = 0.5*N*Pi02*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi203 = 0.5*N*Pi03*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi211 = 0.5*N*Pi11*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi212 = 0.5*N*Pi12*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi213 = 0.5*N*Pi13*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi222 = 0.5*N*Pi22*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi223 = 0.5*N*Pi23*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi233 = 0.5*N*Pi33*(Phi200*Power(t0,2) + 2*Phi201*t0*t1 + Phi211*Power(t1,2) + 2*Phi202*t0*t2 + 2*Phi212*t1*t2 + 
        Phi222*Power(t2,2) + 2*Phi203*t0*t3 + 2*Phi213*t1*t3 + 2*Phi223*t2*t3 + Phi233*Power(t3,2));

t1rhsPhi300 = 0.5*N*Pi00*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
        Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi301 = 0.5*N*Pi01*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
        Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi302 = 0.5*N*Pi02*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
        Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi303 = 0.5*N*Pi03*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
        Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi311 = 0.5*N*Pi11*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
         Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi312 = 0.5*N*Pi12*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi313 = 0.5*N*Pi13*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi322 = 0.5*N*Pi22*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi323 = 0.5*N*Pi23*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));

t1rhsPhi333 = 0.5*N*Pi33*(Phi300*Power(t0,2) + 2*Phi301*t0*t1 + Phi311*Power(t1,2) + 2*Phi302*t0*t2 + 2*Phi312*t1*t2 + 
                Phi322*Power(t2,2) + 2*Phi303*t0*t3 + 2*Phi313*t1*t3 + 2*Phi323*t2*t3 + Phi333*Power(t3,2));


//TERM 2 of Phi
t2rhsPhi100 = N*(g11*Power(Phi100,2)*t0 + g12*Phi100*Phi101*t0 + g13*Phi100*Phi102*t0 + g12*Phi100*Phi200*t0 + 
        g22*Phi101*Phi200*t0 + g23*Phi102*Phi200*t0 + g13*Phi100*Phi300*t0 + g23*Phi101*Phi300*t0 + g33*Phi102*Phi300*t0 + 
        g11*Phi100*Phi101*t1 + g12*Phi100*Phi111*t1 + g13*Phi100*Phi112*t1 + g12*Phi101*Phi200*t1 + g22*Phi111*Phi200*t1 + 
        g23*Phi112*Phi200*t1 + g13*Phi101*Phi300*t1 + g23*Phi111*Phi300*t1 + g33*Phi112*Phi300*t1 + g11*Phi100*Phi102*t2 + 
        g12*Phi100*Phi112*t2 + g13*Phi100*Phi122*t2 + g12*Phi102*Phi200*t2 + g22*Phi112*Phi200*t2 + g23*Phi122*Phi200*t2 + 
        g13*Phi102*Phi300*t2 + g23*Phi112*Phi300*t2 + g33*Phi122*Phi300*t2 + g11*Phi100*Phi103*t3 + g12*Phi100*Phi113*t3 + 
        g13*Phi100*Phi123*t3 + g12*Phi103*Phi200*t3 + g22*Phi113*Phi200*t3 + g23*Phi123*Phi200*t3 + g13*Phi103*Phi300*t3 + 
        g23*Phi113*Phi300*t3 + g33*Phi123*Phi300*t3);

t2rhsPhi101 = N*(g11*Phi100*Phi101*t0 + g12*Power(Phi101,2)*t0 + g13*Phi101*Phi102*t0 + g12*Phi100*Phi201*t0 + g22*Phi101*Phi201*t0 + 
        g23*Phi102*Phi201*t0 + g13*Phi100*Phi301*t0 + g23*Phi101*Phi301*t0 + g33*Phi102*Phi301*t0 + g11*Power(Phi101,2)*t1 + 
        g12*Phi101*Phi111*t1 + g13*Phi101*Phi112*t1 + g12*Phi101*Phi201*t1 + g22*Phi111*Phi201*t1 + g23*Phi112*Phi201*t1 + 
        g13*Phi101*Phi301*t1 + g23*Phi111*Phi301*t1 + g33*Phi112*Phi301*t1 + g11*Phi101*Phi102*t2 + g12*Phi101*Phi112*t2 + 
        g13*Phi101*Phi122*t2 + g12*Phi102*Phi201*t2 + g22*Phi112*Phi201*t2 + g23*Phi122*Phi201*t2 + g13*Phi102*Phi301*t2 + 
        g23*Phi112*Phi301*t2 + g33*Phi122*Phi301*t2 + g11*Phi101*Phi103*t3 + g12*Phi101*Phi113*t3 + g13*Phi101*Phi123*t3 + 
        g12*Phi103*Phi201*t3 + g22*Phi113*Phi201*t3 + g23*Phi123*Phi201*t3 + g13*Phi103*Phi301*t3 + g23*Phi113*Phi301*t3 + 
        g33*Phi123*Phi301*t3);

t2rhsPhi102 = N*(g11*Phi100*Phi102*t0 + g12*Phi101*Phi102*t0 + g13*Power(Phi102,2)*t0 + g12*Phi100*Phi202*t0 + 
        g22*Phi101*Phi202*t0 + g23*Phi102*Phi202*t0 + g13*Phi100*Phi302*t0 + g23*Phi101*Phi302*t0 + g33*Phi102*Phi302*t0 + 
        g11*Phi101*Phi102*t1 + g12*Phi102*Phi111*t1 + g13*Phi102*Phi112*t1 + g12*Phi101*Phi202*t1 + g22*Phi111*Phi202*t1 + 
        g23*Phi112*Phi202*t1 + g13*Phi101*Phi302*t1 + g23*Phi111*Phi302*t1 + g33*Phi112*Phi302*t1 + g11*Power(Phi102,2)*t2 + 
        g12*Phi102*Phi112*t2 + g13*Phi102*Phi122*t2 + g12*Phi102*Phi202*t2 + g22*Phi112*Phi202*t2 + g23*Phi122*Phi202*t2 + 
        g13*Phi102*Phi302*t2 + g23*Phi112*Phi302*t2 + g33*Phi122*Phi302*t2 + g11*Phi102*Phi103*t3 + g12*Phi102*Phi113*t3 + 
        g13*Phi102*Phi123*t3 + g12*Phi103*Phi202*t3 + g22*Phi113*Phi202*t3 + g23*Phi123*Phi202*t3 + g13*Phi103*Phi302*t3 + 
        g23*Phi113*Phi302*t3 + g33*Phi123*Phi302*t3);

t2rhsPhi103 = N*(g11*Phi100*Phi103*t0 + g12*Phi101*Phi103*t0 + g13*Phi102*Phi103*t0 + g12*Phi100*Phi203*t0 + g22*Phi101*Phi203*t0 + 
        g23*Phi102*Phi203*t0 + g13*Phi100*Phi303*t0 + g23*Phi101*Phi303*t0 + g33*Phi102*Phi303*t0 + g11*Phi101*Phi103*t1 + 
        g12*Phi103*Phi111*t1 + g13*Phi103*Phi112*t1 + g12*Phi101*Phi203*t1 + g22*Phi111*Phi203*t1 + g23*Phi112*Phi203*t1 + 
        g13*Phi101*Phi303*t1 + g23*Phi111*Phi303*t1 + g33*Phi112*Phi303*t1 + g11*Phi102*Phi103*t2 + g12*Phi103*Phi112*t2 + 
        g13*Phi103*Phi122*t2 + g12*Phi102*Phi203*t2 + g22*Phi112*Phi203*t2 + g23*Phi122*Phi203*t2 + g13*Phi102*Phi303*t2 + 
        g23*Phi112*Phi303*t2 + g33*Phi122*Phi303*t2 + g11*Power(Phi103,2)*t3 + g12*Phi103*Phi113*t3 + g13*Phi103*Phi123*t3 + 
        g12*Phi103*Phi203*t3 + g22*Phi113*Phi203*t3 + g23*Phi123*Phi203*t3 + g13*Phi103*Phi303*t3 + g23*Phi113*Phi303*t3 + 
        g33*Phi123*Phi303*t3);
     
t2rhsPhi111 = N*(g11*Phi100*Phi111*t0 + g12*Phi101*Phi111*t0 + g13*Phi102*Phi111*t0 + g12*Phi100*Phi211*t0 + g22*Phi101*Phi211*t0 + 
        g23*Phi102*Phi211*t0 + g13*Phi100*Phi311*t0 + g23*Phi101*Phi311*t0 + g33*Phi102*Phi311*t0 + g11*Phi101*Phi111*t1 + 
        g12*Power(Phi111,2)*t1 + g13*Phi111*Phi112*t1 + g12*Phi101*Phi211*t1 + g22*Phi111*Phi211*t1 + g23*Phi112*Phi211*t1 + 
        g13*Phi101*Phi311*t1 + g23*Phi111*Phi311*t1 + g33*Phi112*Phi311*t1 + g11*Phi102*Phi111*t2 + g12*Phi111*Phi112*t2 + 
        g13*Phi111*Phi122*t2 + g12*Phi102*Phi211*t2 + g22*Phi112*Phi211*t2 + g23*Phi122*Phi211*t2 + g13*Phi102*Phi311*t2 + 
        g23*Phi112*Phi311*t2 + g33*Phi122*Phi311*t2 + g11*Phi103*Phi111*t3 + g12*Phi111*Phi113*t3 + g13*Phi111*Phi123*t3 + 
        g12*Phi103*Phi211*t3 + g22*Phi113*Phi211*t3 + g23*Phi123*Phi211*t3 + g13*Phi103*Phi311*t3 + g23*Phi113*Phi311*t3 + 
        g33*Phi123*Phi311*t3);

t2rhsPhi112 = N*(g11*Phi100*Phi112*t0 + g12*Phi101*Phi112*t0 + g13*Phi102*Phi112*t0 + g12*Phi100*Phi212*t0 + 
        g22*Phi101*Phi212*t0 + g23*Phi102*Phi212*t0 + g13*Phi100*Phi312*t0 + g23*Phi101*Phi312*t0 + g33*Phi102*Phi312*t0 + 
        g11*Phi101*Phi112*t1 + g12*Phi111*Phi112*t1 + g13*Power(Phi112,2)*t1 + g12*Phi101*Phi212*t1 + g22*Phi111*Phi212*t1 + 
        g23*Phi112*Phi212*t1 + g13*Phi101*Phi312*t1 + g23*Phi111*Phi312*t1 + g33*Phi112*Phi312*t1 + g11*Phi102*Phi112*t2 + 
        g12*Power(Phi112,2)*t2 + g13*Phi112*Phi122*t2 + g12*Phi102*Phi212*t2 + g22*Phi112*Phi212*t2 + g23*Phi122*Phi212*t2 + 
        g13*Phi102*Phi312*t2 + g23*Phi112*Phi312*t2 + g33*Phi122*Phi312*t2 + g11*Phi103*Phi112*t3 + g12*Phi112*Phi113*t3 + 
        g13*Phi112*Phi123*t3 + g12*Phi103*Phi212*t3 + g22*Phi113*Phi212*t3 + g23*Phi123*Phi212*t3 + g13*Phi103*Phi312*t3 + 
        g23*Phi113*Phi312*t3 + g33*Phi123*Phi312*t3);

t2rhsPhi113 = N*(g11*Phi100*Phi113*t0 + g12*Phi101*Phi113*t0 + g13*Phi102*Phi113*t0 + g12*Phi100*Phi213*t0 + g22*Phi101*Phi213*t0 + 
        g23*Phi102*Phi213*t0 + g13*Phi100*Phi313*t0 + g23*Phi101*Phi313*t0 + g33*Phi102*Phi313*t0 + g11*Phi101*Phi113*t1 + 
        g12*Phi111*Phi113*t1 + g13*Phi112*Phi113*t1 + g12*Phi101*Phi213*t1 + g22*Phi111*Phi213*t1 + g23*Phi112*Phi213*t1 + 
        g13*Phi101*Phi313*t1 + g23*Phi111*Phi313*t1 + g33*Phi112*Phi313*t1 + g11*Phi102*Phi113*t2 + g12*Phi112*Phi113*t2 + 
        g13*Phi113*Phi122*t2 + g12*Phi102*Phi213*t2 + g22*Phi112*Phi213*t2 + g23*Phi122*Phi213*t2 + g13*Phi102*Phi313*t2 + 
        g23*Phi112*Phi313*t2 + g33*Phi122*Phi313*t2 + g11*Phi103*Phi113*t3 + g12*Power(Phi113,2)*t3 + g13*Phi113*Phi123*t3 + 
        g12*Phi103*Phi213*t3 + g22*Phi113*Phi213*t3 + g23*Phi123*Phi213*t3 + g13*Phi103*Phi313*t3 + g23*Phi113*Phi313*t3 + 
        g33*Phi123*Phi313*t3);

t2rhsPhi122 = N*(g11*Phi100*Phi122*t0 + g12*Phi101*Phi122*t0 + g13*Phi102*Phi122*t0 + g12*Phi100*Phi222*t0 + 
        g22*Phi101*Phi222*t0 + g23*Phi102*Phi222*t0 + g13*Phi100*Phi322*t0 + g23*Phi101*Phi322*t0 + g33*Phi102*Phi322*t0 + 
        g11*Phi101*Phi122*t1 + g12*Phi111*Phi122*t1 + g13*Phi112*Phi122*t1 + g12*Phi101*Phi222*t1 + g22*Phi111*Phi222*t1 + 
        g23*Phi112*Phi222*t1 + g13*Phi101*Phi322*t1 + g23*Phi111*Phi322*t1 + g33*Phi112*Phi322*t1 + g11*Phi102*Phi122*t2 + 
        g12*Phi112*Phi122*t2 + g13*Power(Phi122,2)*t2 + g12*Phi102*Phi222*t2 + g22*Phi112*Phi222*t2 + g23*Phi122*Phi222*t2 + 
        g13*Phi102*Phi322*t2 + g23*Phi112*Phi322*t2 + g33*Phi122*Phi322*t2 + g11*Phi103*Phi122*t3 + g12*Phi113*Phi122*t3 + 
        g13*Phi122*Phi123*t3 + g12*Phi103*Phi222*t3 + g22*Phi113*Phi222*t3 + g23*Phi123*Phi222*t3 + g13*Phi103*Phi322*t3 + 
        g23*Phi113*Phi322*t3 + g33*Phi123*Phi322*t3);

t2rhsPhi123 = N*(g11*Phi100*Phi123*t0 + g12*Phi101*Phi123*t0 + g13*Phi102*Phi123*t0 + g12*Phi100*Phi223*t0 + g22*Phi101*Phi223*t0 + 
        g23*Phi102*Phi223*t0 + g13*Phi100*Phi323*t0 + g23*Phi101*Phi323*t0 + g33*Phi102*Phi323*t0 + g11*Phi101*Phi123*t1 + 
        g12*Phi111*Phi123*t1 + g13*Phi112*Phi123*t1 + g12*Phi101*Phi223*t1 + g22*Phi111*Phi223*t1 + g23*Phi112*Phi223*t1 + 
        g13*Phi101*Phi323*t1 + g23*Phi111*Phi323*t1 + g33*Phi112*Phi323*t1 + g11*Phi102*Phi123*t2 + g12*Phi112*Phi123*t2 + 
        g13*Phi122*Phi123*t2 + g12*Phi102*Phi223*t2 + g22*Phi112*Phi223*t2 + g23*Phi122*Phi223*t2 + g13*Phi102*Phi323*t2 + 
        g23*Phi112*Phi323*t2 + g33*Phi122*Phi323*t2 + g11*Phi103*Phi123*t3 + g12*Phi113*Phi123*t3 + g13*Power(Phi123,2)*t3 + 
        g12*Phi103*Phi223*t3 + g22*Phi113*Phi223*t3 + g23*Phi123*Phi223*t3 + g13*Phi103*Phi323*t3 + g23*Phi113*Phi323*t3 + 
        g33*Phi123*Phi323*t3);

t2rhsPhi133 = N*(g11*Phi100*Phi133*t0 + g12*Phi101*Phi133*t0 + g13*Phi102*Phi133*t0 + g12*Phi100*Phi233*t0 + g22*Phi101*Phi233*t0 + 
        g23*Phi102*Phi233*t0 + g13*Phi100*Phi333*t0 + g23*Phi101*Phi333*t0 + g33*Phi102*Phi333*t0 + g11*Phi101*Phi133*t1 + 
        g12*Phi111*Phi133*t1 + g13*Phi112*Phi133*t1 + g12*Phi101*Phi233*t1 + g22*Phi111*Phi233*t1 + g23*Phi112*Phi233*t1 + 
        g13*Phi101*Phi333*t1 + g23*Phi111*Phi333*t1 + g33*Phi112*Phi333*t1 + g11*Phi102*Phi133*t2 + g12*Phi112*Phi133*t2 + 
        g13*Phi122*Phi133*t2 + g12*Phi102*Phi233*t2 + g22*Phi112*Phi233*t2 + g23*Phi122*Phi233*t2 + g13*Phi102*Phi333*t2 + 
        g23*Phi112*Phi333*t2 + g33*Phi122*Phi333*t2 + g11*Phi103*Phi133*t3 + g12*Phi113*Phi133*t3 + g13*Phi123*Phi133*t3 + 
        g12*Phi103*Phi233*t3 + g22*Phi113*Phi233*t3 + g23*Phi123*Phi233*t3 + g13*Phi103*Phi333*t3 + g23*Phi113*Phi333*t3 + 
        g33*Phi123*Phi333*t3);

t2rhsPhi200 = N*(g11*Phi100*Phi200*t0 + g12*Power(Phi200,2)*t0 + g12*Phi100*Phi201*t0 + 
        g22*Phi200*Phi201*t0 + g13*Phi100*Phi202*t0 + g23*Phi200*Phi202*t0 + g13*Phi200*Phi300*t0 + g23*Phi201*Phi300*t0 + 
        g33*Phi202*Phi300*t0 + g11*Phi100*Phi201*t1 + g12*Phi200*Phi201*t1 + g12*Phi100*Phi211*t1 + g22*Phi200*Phi211*t1 + 
        g13*Phi100*Phi212*t1 + g23*Phi200*Phi212*t1 + g13*Phi201*Phi300*t1 + g23*Phi211*Phi300*t1 + g33*Phi212*Phi300*t1 + 
        g11*Phi100*Phi202*t2 + g12*Phi200*Phi202*t2 + g12*Phi100*Phi212*t2 + g22*Phi200*Phi212*t2 + g13*Phi100*Phi222*t2 + 
        g23*Phi200*Phi222*t2 + g13*Phi202*Phi300*t2 + g23*Phi212*Phi300*t2 + g33*Phi222*Phi300*t2 + g11*Phi100*Phi203*t3 + 
        g12*Phi200*Phi203*t3 + g12*Phi100*Phi213*t3 + g22*Phi200*Phi213*t3 + g13*Phi100*Phi223*t3 + g23*Phi200*Phi223*t3 + 
        g13*Phi203*Phi300*t3 + g23*Phi213*Phi300*t3 + g33*Phi223*Phi300*t3);

t2rhsPhi201 = N*(g11*Phi101*Phi200*t0 + g12*Phi101*Phi201*t0 + g12*Phi200*Phi201*t0 + g22*Power(Phi201,2)*t0 + g13*Phi101*Phi202*t0 + 
        g23*Phi201*Phi202*t0 + g13*Phi200*Phi301*t0 + g23*Phi201*Phi301*t0 + g33*Phi202*Phi301*t0 + g11*Phi101*Phi201*t1 + 
        g12*Power(Phi201,2)*t1 + g12*Phi101*Phi211*t1 + g22*Phi201*Phi211*t1 + g13*Phi101*Phi212*t1 + g23*Phi201*Phi212*t1 + 
        g13*Phi201*Phi301*t1 + g23*Phi211*Phi301*t1 + g33*Phi212*Phi301*t1 + g11*Phi101*Phi202*t2 + g12*Phi201*Phi202*t2 + 
        g12*Phi101*Phi212*t2 + g22*Phi201*Phi212*t2 + g13*Phi101*Phi222*t2 + g23*Phi201*Phi222*t2 + g13*Phi202*Phi301*t2 + 
        g23*Phi212*Phi301*t2 + g33*Phi222*Phi301*t2 + g11*Phi101*Phi203*t3 + g12*Phi201*Phi203*t3 + g12*Phi101*Phi213*t3 + 
        g22*Phi201*Phi213*t3 + g13*Phi101*Phi223*t3 + g23*Phi201*Phi223*t3 + g13*Phi203*Phi301*t3 + g23*Phi213*Phi301*t3 + 
        g33*Phi223*Phi301*t3);

t2rhsPhi202 = N*(g11*Phi102*Phi200*t0 + g12*Phi102*Phi201*t0 + g13*Phi102*Phi202*t0 + g12*Phi200*Phi202*t0 + 
        g22*Phi201*Phi202*t0 + g23*Power(Phi202,2)*t0 + g13*Phi200*Phi302*t0 + g23*Phi201*Phi302*t0 + g33*Phi202*Phi302*t0 + 
        g11*Phi102*Phi201*t1 + g12*Phi201*Phi202*t1 + g12*Phi102*Phi211*t1 + g22*Phi202*Phi211*t1 + g13*Phi102*Phi212*t1 + 
        g23*Phi202*Phi212*t1 + g13*Phi201*Phi302*t1 + g23*Phi211*Phi302*t1 + g33*Phi212*Phi302*t1 + g11*Phi102*Phi202*t2 + 
        g12*Power(Phi202,2)*t2 + g12*Phi102*Phi212*t2 + g22*Phi202*Phi212*t2 + g13*Phi102*Phi222*t2 + g23*Phi202*Phi222*t2 + 
        g13*Phi202*Phi302*t2 + g23*Phi212*Phi302*t2 + g33*Phi222*Phi302*t2 + g11*Phi102*Phi203*t3 + g12*Phi202*Phi203*t3 + 
        g12*Phi102*Phi213*t3 + g22*Phi202*Phi213*t3 + g13*Phi102*Phi223*t3 + g23*Phi202*Phi223*t3 + g13*Phi203*Phi302*t3 + 
        g23*Phi213*Phi302*t3 + g33*Phi223*Phi302*t3);

t2rhsPhi203 = N*(g11*Phi103*Phi200*t0 + g12*Phi103*Phi201*t0 + g13*Phi103*Phi202*t0 + g12*Phi200*Phi203*t0 + g22*Phi201*Phi203*t0 + 
        g23*Phi202*Phi203*t0 + g13*Phi200*Phi303*t0 + g23*Phi201*Phi303*t0 + g33*Phi202*Phi303*t0 + g11*Phi103*Phi201*t1 + 
        g12*Phi201*Phi203*t1 + g12*Phi103*Phi211*t1 + g22*Phi203*Phi211*t1 + g13*Phi103*Phi212*t1 + g23*Phi203*Phi212*t1 + 
        g13*Phi201*Phi303*t1 + g23*Phi211*Phi303*t1 + g33*Phi212*Phi303*t1 + g11*Phi103*Phi202*t2 + g12*Phi202*Phi203*t2 + 
        g12*Phi103*Phi212*t2 + g22*Phi203*Phi212*t2 + g13*Phi103*Phi222*t2 + g23*Phi203*Phi222*t2 + g13*Phi202*Phi303*t2 + 
        g23*Phi212*Phi303*t2 + g33*Phi222*Phi303*t2 + g11*Phi103*Phi203*t3 + g12*Power(Phi203,2)*t3 + g12*Phi103*Phi213*t3 + 
        g22*Phi203*Phi213*t3 + g13*Phi103*Phi223*t3 + g23*Phi203*Phi223*t3 + g13*Phi203*Phi303*t3 + g23*Phi213*Phi303*t3 + 
        g33*Phi223*Phi303*t3);

t2rhsPhi211 = N*(g11*Phi111*Phi200*t0 + g12*Phi111*Phi201*t0 + g13*Phi111*Phi202*t0 + g12*Phi200*Phi211*t0 + g22*Phi201*Phi211*t0 + 
        g23*Phi202*Phi211*t0 + g13*Phi200*Phi311*t0 + g23*Phi201*Phi311*t0 + g33*Phi202*Phi311*t0 + g11*Phi111*Phi201*t1 + 
        g12*Phi111*Phi211*t1 + g12*Phi201*Phi211*t1 + g22*Power(Phi211,2)*t1 + g13*Phi111*Phi212*t1 + g23*Phi211*Phi212*t1 + 
        g13*Phi201*Phi311*t1 + g23*Phi211*Phi311*t1 + g33*Phi212*Phi311*t1 + g11*Phi111*Phi202*t2 + g12*Phi202*Phi211*t2 + 
        g12*Phi111*Phi212*t2 + g22*Phi211*Phi212*t2 + g13*Phi111*Phi222*t2 + g23*Phi211*Phi222*t2 + g13*Phi202*Phi311*t2 + 
        g23*Phi212*Phi311*t2 + g33*Phi222*Phi311*t2 + g11*Phi111*Phi203*t3 + g12*Phi203*Phi211*t3 + g12*Phi111*Phi213*t3 + 
        g22*Phi211*Phi213*t3 + g13*Phi111*Phi223*t3 + g23*Phi211*Phi223*t3 + g13*Phi203*Phi311*t3 + g23*Phi213*Phi311*t3 + 
        g33*Phi223*Phi311*t3);

t2rhsPhi212 = N*(g11*Phi112*Phi200*t0 + g12*Phi112*Phi201*t0 + g13*Phi112*Phi202*t0 + g12*Phi200*Phi212*t0 + 
        g22*Phi201*Phi212*t0 + g23*Phi202*Phi212*t0 + g13*Phi200*Phi312*t0 + g23*Phi201*Phi312*t0 + g33*Phi202*Phi312*t0 + 
        g11*Phi112*Phi201*t1 + g12*Phi112*Phi211*t1 + g13*Phi112*Phi212*t1 + g12*Phi201*Phi212*t1 + g22*Phi211*Phi212*t1 + 
        g23*Power(Phi212,2)*t1 + g13*Phi201*Phi312*t1 + g23*Phi211*Phi312*t1 + g33*Phi212*Phi312*t1 + g11*Phi112*Phi202*t2 + 
        g12*Phi112*Phi212*t2 + g12*Phi202*Phi212*t2 + g22*Power(Phi212,2)*t2 + g13*Phi112*Phi222*t2 + g23*Phi212*Phi222*t2 + 
        g13*Phi202*Phi312*t2 + g23*Phi212*Phi312*t2 + g33*Phi222*Phi312*t2 + g11*Phi112*Phi203*t3 + g12*Phi203*Phi212*t3 + 
        g12*Phi112*Phi213*t3 + g22*Phi212*Phi213*t3 + g13*Phi112*Phi223*t3 + g23*Phi212*Phi223*t3 + g13*Phi203*Phi312*t3 + 
        g23*Phi213*Phi312*t3 + g33*Phi223*Phi312*t3);

t2rhsPhi213 = N*(g11*Phi113*Phi200*t0 + g12*Phi113*Phi201*t0 + g13*Phi113*Phi202*t0 + g12*Phi200*Phi213*t0 + g22*Phi201*Phi213*t0 + 
        g23*Phi202*Phi213*t0 + g13*Phi200*Phi313*t0 + g23*Phi201*Phi313*t0 + g33*Phi202*Phi313*t0 + g11*Phi113*Phi201*t1 + 
        g12*Phi113*Phi211*t1 + g13*Phi113*Phi212*t1 + g12*Phi201*Phi213*t1 + g22*Phi211*Phi213*t1 + g23*Phi212*Phi213*t1 + 
        g13*Phi201*Phi313*t1 + g23*Phi211*Phi313*t1 + g33*Phi212*Phi313*t1 + g11*Phi113*Phi202*t2 + g12*Phi113*Phi212*t2 + 
        g12*Phi202*Phi213*t2 + g22*Phi212*Phi213*t2 + g13*Phi113*Phi222*t2 + g23*Phi213*Phi222*t2 + g13*Phi202*Phi313*t2 + 
        g23*Phi212*Phi313*t2 + g33*Phi222*Phi313*t2 + g11*Phi113*Phi203*t3 + g12*Phi113*Phi213*t3 + g12*Phi203*Phi213*t3 + 
        g22*Power(Phi213,2)*t3 + g13*Phi113*Phi223*t3 + g23*Phi213*Phi223*t3 + g13*Phi203*Phi313*t3 + g23*Phi213*Phi313*t3 + 
        g33*Phi223*Phi313*t3);

t2rhsPhi222 = N*(g11*Phi122*Phi200*t0 + g12*Phi122*Phi201*t0 + g13*Phi122*Phi202*t0 + g12*Phi200*Phi222*t0 + 
        g22*Phi201*Phi222*t0 + g23*Phi202*Phi222*t0 + g13*Phi200*Phi322*t0 + g23*Phi201*Phi322*t0 + g33*Phi202*Phi322*t0 + 
        g11*Phi122*Phi201*t1 + g12*Phi122*Phi211*t1 + g13*Phi122*Phi212*t1 + g12*Phi201*Phi222*t1 + g22*Phi211*Phi222*t1 + 
        g23*Phi212*Phi222*t1 + g13*Phi201*Phi322*t1 + g23*Phi211*Phi322*t1 + g33*Phi212*Phi322*t1 + g11*Phi122*Phi202*t2 + 
        g12*Phi122*Phi212*t2 + g13*Phi122*Phi222*t2 + g12*Phi202*Phi222*t2 + g22*Phi212*Phi222*t2 + g23*Power(Phi222,2)*t2 + 
        g13*Phi202*Phi322*t2 + g23*Phi212*Phi322*t2 + g33*Phi222*Phi322*t2 + g11*Phi122*Phi203*t3 + g12*Phi122*Phi213*t3 + 
        g12*Phi203*Phi222*t3 + g22*Phi213*Phi222*t3 + g13*Phi122*Phi223*t3 + g23*Phi222*Phi223*t3 + g13*Phi203*Phi322*t3 + 
        g23*Phi213*Phi322*t3 + g33*Phi223*Phi322*t3);

t2rhsPhi223 = N*(g11*Phi123*Phi200*t0 + g12*Phi123*Phi201*t0 + g13*Phi123*Phi202*t0 + g12*Phi200*Phi223*t0 + g22*Phi201*Phi223*t0 + 
        g23*Phi202*Phi223*t0 + g13*Phi200*Phi323*t0 + g23*Phi201*Phi323*t0 + g33*Phi202*Phi323*t0 + g11*Phi123*Phi201*t1 + 
        g12*Phi123*Phi211*t1 + g13*Phi123*Phi212*t1 + g12*Phi201*Phi223*t1 + g22*Phi211*Phi223*t1 + g23*Phi212*Phi223*t1 + 
        g13*Phi201*Phi323*t1 + g23*Phi211*Phi323*t1 + g33*Phi212*Phi323*t1 + g11*Phi123*Phi202*t2 + g12*Phi123*Phi212*t2 + 
        g13*Phi123*Phi222*t2 + g12*Phi202*Phi223*t2 + g22*Phi212*Phi223*t2 + g23*Phi222*Phi223*t2 + g13*Phi202*Phi323*t2 + 
        g23*Phi212*Phi323*t2 + g33*Phi222*Phi323*t2 + g11*Phi123*Phi203*t3 + g12*Phi123*Phi213*t3 + g13*Phi123*Phi223*t3 + 
        g12*Phi203*Phi223*t3 + g22*Phi213*Phi223*t3 + g23*Power(Phi223,2)*t3 + g13*Phi203*Phi323*t3 + g23*Phi213*Phi323*t3 + 
        g33*Phi223*Phi323*t3);

t2rhsPhi233 = N*(g11*Phi133*Phi200*t0 + g12*Phi133*Phi201*t0 + g13*Phi133*Phi202*t0 + g12*Phi200*Phi233*t0 + g22*Phi201*Phi233*t0 + 
        g23*Phi202*Phi233*t0 + g13*Phi200*Phi333*t0 + g23*Phi201*Phi333*t0 + g33*Phi202*Phi333*t0 + g11*Phi133*Phi201*t1 + 
        g12*Phi133*Phi211*t1 + g13*Phi133*Phi212*t1 + g12*Phi201*Phi233*t1 + g22*Phi211*Phi233*t1 + g23*Phi212*Phi233*t1 + 
        g13*Phi201*Phi333*t1 + g23*Phi211*Phi333*t1 + g33*Phi212*Phi333*t1 + g11*Phi133*Phi202*t2 + g12*Phi133*Phi212*t2 + 
        g13*Phi133*Phi222*t2 + g12*Phi202*Phi233*t2 + g22*Phi212*Phi233*t2 + g23*Phi222*Phi233*t2 + g13*Phi202*Phi333*t2 + 
        g23*Phi212*Phi333*t2 + g33*Phi222*Phi333*t2 + g11*Phi133*Phi203*t3 + g12*Phi133*Phi213*t3 + g13*Phi133*Phi223*t3 + 
        g12*Phi203*Phi233*t3 + g22*Phi213*Phi233*t3 + g23*Phi223*Phi233*t3 + g13*Phi203*Phi333*t3 + g23*Phi213*Phi333*t3 + 
        g33*Phi223*Phi333*t3);

t2rhsPhi300 = N*(g11*Phi100*Phi300*t0 + g12*Phi200*Phi300*t0 + g13*Power(Phi300,2)*t0 + 
        g12*Phi100*Phi301*t0 + g22*Phi200*Phi301*t0 + g23*Phi300*Phi301*t0 + g13*Phi100*Phi302*t0 + g23*Phi200*Phi302*t0 + 
        g33*Phi300*Phi302*t0 + g11*Phi100*Phi301*t1 + g12*Phi200*Phi301*t1 + g13*Phi300*Phi301*t1 + g12*Phi100*Phi311*t1 + 
        g22*Phi200*Phi311*t1 + g23*Phi300*Phi311*t1 + g13*Phi100*Phi312*t1 + g23*Phi200*Phi312*t1 + g33*Phi300*Phi312*t1 + 
        g11*Phi100*Phi302*t2 + g12*Phi200*Phi302*t2 + g13*Phi300*Phi302*t2 + g12*Phi100*Phi312*t2 + g22*Phi200*Phi312*t2 + 
        g23*Phi300*Phi312*t2 + g13*Phi100*Phi322*t2 + g23*Phi200*Phi322*t2 + g33*Phi300*Phi322*t2 + g11*Phi100*Phi303*t3 + 
        g12*Phi200*Phi303*t3 + g13*Phi300*Phi303*t3 + g12*Phi100*Phi313*t3 + g22*Phi200*Phi313*t3 + g23*Phi300*Phi313*t3 + 
        g13*Phi100*Phi323*t3 + g23*Phi200*Phi323*t3 + g33*Phi300*Phi323*t3);

t2rhsPhi301 = N*(g11*Phi101*Phi300*t0 + g12*Phi201*Phi300*t0 + g12*Phi101*Phi301*t0 + g22*Phi201*Phi301*t0 + g13*Phi300*Phi301*t0 + 
        g23*Power(Phi301,2)*t0 + g13*Phi101*Phi302*t0 + g23*Phi201*Phi302*t0 + g33*Phi301*Phi302*t0 + g11*Phi101*Phi301*t1 + 
        g12*Phi201*Phi301*t1 + g13*Power(Phi301,2)*t1 + g12*Phi101*Phi311*t1 + g22*Phi201*Phi311*t1 + g23*Phi301*Phi311*t1 + 
        g13*Phi101*Phi312*t1 + g23*Phi201*Phi312*t1 + g33*Phi301*Phi312*t1 + g11*Phi101*Phi302*t2 + g12*Phi201*Phi302*t2 + 
        g13*Phi301*Phi302*t2 + g12*Phi101*Phi312*t2 + g22*Phi201*Phi312*t2 + g23*Phi301*Phi312*t2 + g13*Phi101*Phi322*t2 + 
        g23*Phi201*Phi322*t2 + g33*Phi301*Phi322*t2 + g11*Phi101*Phi303*t3 + g12*Phi201*Phi303*t3 + g13*Phi301*Phi303*t3 + 
        g12*Phi101*Phi313*t3 + g22*Phi201*Phi313*t3 + g23*Phi301*Phi313*t3 + g13*Phi101*Phi323*t3 + g23*Phi201*Phi323*t3 + 
        g33*Phi301*Phi323*t3);

t2rhsPhi302 = N*(g11*Phi102*Phi300*t0 + g12*Phi202*Phi300*t0 + g12*Phi102*Phi301*t0 + g22*Phi202*Phi301*t0 + 
        g13*Phi102*Phi302*t0 + g23*Phi202*Phi302*t0 + g13*Phi300*Phi302*t0 + g23*Phi301*Phi302*t0 + g33*Power(Phi302,2)*t0 + 
        g11*Phi102*Phi301*t1 + g12*Phi202*Phi301*t1 + g13*Phi301*Phi302*t1 + g12*Phi102*Phi311*t1 + g22*Phi202*Phi311*t1 + 
        g23*Phi302*Phi311*t1 + g13*Phi102*Phi312*t1 + g23*Phi202*Phi312*t1 + g33*Phi302*Phi312*t1 + g11*Phi102*Phi302*t2 + 
        g12*Phi202*Phi302*t2 + g13*Power(Phi302,2)*t2 + g12*Phi102*Phi312*t2 + g22*Phi202*Phi312*t2 + g23*Phi302*Phi312*t2 + 
        g13*Phi102*Phi322*t2 + g23*Phi202*Phi322*t2 + g33*Phi302*Phi322*t2 + g11*Phi102*Phi303*t3 + g12*Phi202*Phi303*t3 + 
        g13*Phi302*Phi303*t3 + g12*Phi102*Phi313*t3 + g22*Phi202*Phi313*t3 + g23*Phi302*Phi313*t3 + g13*Phi102*Phi323*t3 + 
        g23*Phi202*Phi323*t3 + g33*Phi302*Phi323*t3);

t2rhsPhi303 = N*(g11*Phi103*Phi300*t0 + g12*Phi203*Phi300*t0 + g12*Phi103*Phi301*t0 + g22*Phi203*Phi301*t0 + g13*Phi103*Phi302*t0 + 
        g23*Phi203*Phi302*t0 + g13*Phi300*Phi303*t0 + g23*Phi301*Phi303*t0 + g33*Phi302*Phi303*t0 + g11*Phi103*Phi301*t1 + 
        g12*Phi203*Phi301*t1 + g13*Phi301*Phi303*t1 + g12*Phi103*Phi311*t1 + g22*Phi203*Phi311*t1 + g23*Phi303*Phi311*t1 + 
        g13*Phi103*Phi312*t1 + g23*Phi203*Phi312*t1 + g33*Phi303*Phi312*t1 + g11*Phi103*Phi302*t2 + g12*Phi203*Phi302*t2 + 
        g13*Phi302*Phi303*t2 + g12*Phi103*Phi312*t2 + g22*Phi203*Phi312*t2 + g23*Phi303*Phi312*t2 + g13*Phi103*Phi322*t2 + 
        g23*Phi203*Phi322*t2 + g33*Phi303*Phi322*t2 + g11*Phi103*Phi303*t3 + g12*Phi203*Phi303*t3 + g13*Power(Phi303,2)*t3 + 
        g12*Phi103*Phi313*t3 + g22*Phi203*Phi313*t3 + g23*Phi303*Phi313*t3 + g13*Phi103*Phi323*t3 + g23*Phi203*Phi323*t3 + 
        g33*Phi303*Phi323*t3);

t2rhsPhi311 = N*(g11*Phi111*Phi300*t0 + g12*Phi211*Phi300*t0 + g12*Phi111*Phi301*t0 + g22*Phi211*Phi301*t0 + g13*Phi111*Phi302*t0 + 
        g23*Phi211*Phi302*t0 + g13*Phi300*Phi311*t0 + g23*Phi301*Phi311*t0 + g33*Phi302*Phi311*t0 + g11*Phi111*Phi301*t1 + 
        g12*Phi211*Phi301*t1 + g12*Phi111*Phi311*t1 + g22*Phi211*Phi311*t1 + g13*Phi301*Phi311*t1 + g23*Power(Phi311,2)*t1 + 
        g13*Phi111*Phi312*t1 + g23*Phi211*Phi312*t1 + g33*Phi311*Phi312*t1 + g11*Phi111*Phi302*t2 + g12*Phi211*Phi302*t2 + 
        g13*Phi302*Phi311*t2 + g12*Phi111*Phi312*t2 + g22*Phi211*Phi312*t2 + g23*Phi311*Phi312*t2 + g13*Phi111*Phi322*t2 + 
        g23*Phi211*Phi322*t2 + g33*Phi311*Phi322*t2 + g11*Phi111*Phi303*t3 + g12*Phi211*Phi303*t3 + g13*Phi303*Phi311*t3 + 
        g12*Phi111*Phi313*t3 + g22*Phi211*Phi313*t3 + g23*Phi311*Phi313*t3 + g13*Phi111*Phi323*t3 + g23*Phi211*Phi323*t3 + 
        g33*Phi311*Phi323*t3);

t2rhsPhi312 = N*(g11*Phi112*Phi300*t0 + g12*Phi212*Phi300*t0 + g12*Phi112*Phi301*t0 + g22*Phi212*Phi301*t0 + 
        g13*Phi112*Phi302*t0 + g23*Phi212*Phi302*t0 + g13*Phi300*Phi312*t0 + g23*Phi301*Phi312*t0 + g33*Phi302*Phi312*t0 + 
        g11*Phi112*Phi301*t1 + g12*Phi212*Phi301*t1 + g12*Phi112*Phi311*t1 + g22*Phi212*Phi311*t1 + g13*Phi112*Phi312*t1 + 
        g23*Phi212*Phi312*t1 + g13*Phi301*Phi312*t1 + g23*Phi311*Phi312*t1 + g33*Power(Phi312,2)*t1 + g11*Phi112*Phi302*t2 + 
        g12*Phi212*Phi302*t2 + g12*Phi112*Phi312*t2 + g22*Phi212*Phi312*t2 + g13*Phi302*Phi312*t2 + g23*Power(Phi312,2)*t2 + 
        g13*Phi112*Phi322*t2 + g23*Phi212*Phi322*t2 + g33*Phi312*Phi322*t2 + g11*Phi112*Phi303*t3 + g12*Phi212*Phi303*t3 + 
        g13*Phi303*Phi312*t3 + g12*Phi112*Phi313*t3 + g22*Phi212*Phi313*t3 + g23*Phi312*Phi313*t3 + g13*Phi112*Phi323*t3 + 
        g23*Phi212*Phi323*t3 + g33*Phi312*Phi323*t3);

t2rhsPhi313 = N*(g11*Phi113*Phi300*t0 + g12*Phi213*Phi300*t0 + g12*Phi113*Phi301*t0 + g22*Phi213*Phi301*t0 + g13*Phi113*Phi302*t0 + 
        g23*Phi213*Phi302*t0 + g13*Phi300*Phi313*t0 + g23*Phi301*Phi313*t0 + g33*Phi302*Phi313*t0 + g11*Phi113*Phi301*t1 + 
        g12*Phi213*Phi301*t1 + g12*Phi113*Phi311*t1 + g22*Phi213*Phi311*t1 + g13*Phi113*Phi312*t1 + g23*Phi213*Phi312*t1 + 
        g13*Phi301*Phi313*t1 + g23*Phi311*Phi313*t1 + g33*Phi312*Phi313*t1 + g11*Phi113*Phi302*t2 + g12*Phi213*Phi302*t2 + 
        g12*Phi113*Phi312*t2 + g22*Phi213*Phi312*t2 + g13*Phi302*Phi313*t2 + g23*Phi312*Phi313*t2 + g13*Phi113*Phi322*t2 + 
        g23*Phi213*Phi322*t2 + g33*Phi313*Phi322*t2 + g11*Phi113*Phi303*t3 + g12*Phi213*Phi303*t3 + g12*Phi113*Phi313*t3 + 
        g22*Phi213*Phi313*t3 + g13*Phi303*Phi313*t3 + g23*Power(Phi313,2)*t3 + g13*Phi113*Phi323*t3 + g23*Phi213*Phi323*t3 + 
        g33*Phi313*Phi323*t3);

t2rhsPhi322 = N*(g11*Phi122*Phi300*t0 + g12*Phi222*Phi300*t0 + g12*Phi122*Phi301*t0 + g22*Phi222*Phi301*t0 + 
        g13*Phi122*Phi302*t0 + g23*Phi222*Phi302*t0 + g13*Phi300*Phi322*t0 + g23*Phi301*Phi322*t0 + g33*Phi302*Phi322*t0 + 
        g11*Phi122*Phi301*t1 + g12*Phi222*Phi301*t1 + g12*Phi122*Phi311*t1 + g22*Phi222*Phi311*t1 + g13*Phi122*Phi312*t1 + 
        g23*Phi222*Phi312*t1 + g13*Phi301*Phi322*t1 + g23*Phi311*Phi322*t1 + g33*Phi312*Phi322*t1 + g11*Phi122*Phi302*t2 + 
        g12*Phi222*Phi302*t2 + g12*Phi122*Phi312*t2 + g22*Phi222*Phi312*t2 + g13*Phi122*Phi322*t2 + g23*Phi222*Phi322*t2 + 
        g13*Phi302*Phi322*t2 + g23*Phi312*Phi322*t2 + g33*Power(Phi322,2)*t2 + g11*Phi122*Phi303*t3 + g12*Phi222*Phi303*t3 + 
        g12*Phi122*Phi313*t3 + g22*Phi222*Phi313*t3 + g13*Phi303*Phi322*t3 + g23*Phi313*Phi322*t3 + g13*Phi122*Phi323*t3 + 
        g23*Phi222*Phi323*t3 + g33*Phi322*Phi323*t3);

t2rhsPhi323 = N*(g11*Phi123*Phi300*t0 + g12*Phi223*Phi300*t0 + g12*Phi123*Phi301*t0 + g22*Phi223*Phi301*t0 + g13*Phi123*Phi302*t0 + 
        g23*Phi223*Phi302*t0 + g13*Phi300*Phi323*t0 + g23*Phi301*Phi323*t0 + g33*Phi302*Phi323*t0 + g11*Phi123*Phi301*t1 + 
        g12*Phi223*Phi301*t1 + g12*Phi123*Phi311*t1 + g22*Phi223*Phi311*t1 + g13*Phi123*Phi312*t1 + g23*Phi223*Phi312*t1 + 
        g13*Phi301*Phi323*t1 + g23*Phi311*Phi323*t1 + g33*Phi312*Phi323*t1 + g11*Phi123*Phi302*t2 + g12*Phi223*Phi302*t2 + 
        g12*Phi123*Phi312*t2 + g22*Phi223*Phi312*t2 + g13*Phi123*Phi322*t2 + g23*Phi223*Phi322*t2 + g13*Phi302*Phi323*t2 + 
        g23*Phi312*Phi323*t2 + g33*Phi322*Phi323*t2 + g11*Phi123*Phi303*t3 + g12*Phi223*Phi303*t3 + g12*Phi123*Phi313*t3 + 
        g22*Phi223*Phi313*t3 + g13*Phi123*Phi323*t3 + g23*Phi223*Phi323*t3 + g13*Phi303*Phi323*t3 + g23*Phi313*Phi323*t3 + 
        g33*Power(Phi323,2)*t3);

t2rhsPhi333 = N*(g11*Phi133*Phi300*t0 + g12*Phi233*Phi300*t0 + g12*Phi133*Phi301*t0 + g22*Phi233*Phi301*t0 + g13*Phi133*Phi302*t0 + 
        g23*Phi233*Phi302*t0 + g13*Phi300*Phi333*t0 + g23*Phi301*Phi333*t0 + g33*Phi302*Phi333*t0 + g11*Phi133*Phi301*t1 + 
        g12*Phi233*Phi301*t1 + g12*Phi133*Phi311*t1 + g22*Phi233*Phi311*t1 + g13*Phi133*Phi312*t1 + g23*Phi233*Phi312*t1 + 
        g13*Phi301*Phi333*t1 + g23*Phi311*Phi333*t1 + g33*Phi312*Phi333*t1 + g11*Phi133*Phi302*t2 + g12*Phi233*Phi302*t2 + 
        g12*Phi133*Phi312*t2 + g22*Phi233*Phi312*t2 + g13*Phi133*Phi322*t2 + g23*Phi233*Phi322*t2 + g13*Phi302*Phi333*t2 + 
        g23*Phi312*Phi333*t2 + g33*Phi322*Phi333*t2 + g11*Phi133*Phi303*t3 + g12*Phi233*Phi303*t3 + g12*Phi133*Phi313*t3 + 
        g22*Phi233*Phi313*t3 + g13*Phi133*Phi323*t3 + g23*Phi233*Phi323*t3 + g13*Phi303*Phi333*t3 + g23*Phi313*Phi333*t3 + 
        g33*Phi323*Phi333*t3);


//TERM 3 of Phi
t3rhsphi100 = -(gamma2*N*Phi100);
t3rhsPhi101 = -(gamma2*N*Phi101);
t3rhsPhi102 = -(gamma2*N*Phi102);
t3rhsPhi103 = -(gamma2*N*Phi103);
t3rhsPhi111 = -(gamma2*N*Phi111);
t3rhsPhi112 = -(gamma2*N*Phi112);
t3rhsPhi113 = -(gamma2*N*Phi113);
t3rhsPhi122 = -(gamma2*N*Phi122);
t3rhsPhi123 = -(gamma2*N*Phi123);
t3rhsPhi133 = -(gamma2*N*Phi133);

t3rhsPhi200 = -(gamma2*N*Phi200);
t3rhsPhi201 = -(gamma2*N*Phi201);
t3rhsPhi202 = -(gamma2*N*Phi202);
t3rhsPhi203 = -(gamma2*N*Phi203);
t3rhsPhi211 = -(gamma2*N*Phi211);
t3rhsPhi212 = -(gamma2*N*Phi212);
t3rhsPhi213 = -(gamma2*N*Phi213);
t3rhsPhi222 = -(gamma2*N*Phi222);
t3rhsPhi223 = -(gamma2*N*Phi223);
t3rhsPhi233 = -(gamma2*N*Phi233);

t3rhsPhi300 = -(gamma2*N*Phi300);
t3rhsPhi301 = -(gamma2*N*Phi301);
t3rhsPhi302 = -(gamma2*N*Phi302);
t3rhsPhi303 = -(gamma2*N*Phi303);
t3rhsPhi311 = -(gamma2*N*Phi311);
t3rhsPhi312 = -(gamma2*N*Phi312);
t3rhsPhi313 = -(gamma2*N*Phi313);
t3rhsPhi322 = -(gamma2*N*Phi322);
t3rhsPhi323 = -(gamma2*N*Phi323);
t3rhsPhi333 = -(gamma2*N*Phi333);

//Sum all terms of RHS Phi
rhsphi100 = t1rhsPhi100 + t2rhsPhi + t3rhsPhi100;
rhsPhi101 = t1rhsPhi101 + t2rhsPhi + t3rhsPhi101;
rhsPhi102 = t1rhsPhi102 + t2rhsPhi + t3rhsPhi102;
rhsPhi103 = t1rhsPhi103 + t2rhsPhi + t3rhsPhi103;
rhsPhi111 = t1rhsPhi111 + t2rhsPhi + t3rhsPhi111;
rhsPhi112 = t1rhsPhi112 + t2rhsPhi + t3rhsPhi112;
rhsPhi113 = t1rhsPhi113 + t2rhsPhi + t3rhsPhi113;
rhsPhi122 = t1rhsPhi122 + t2rhsPhi + t3rhsPhi122;
rhsPhi123 = t1rhsPhi123 + t2rhsPhi + t3rhsPhi123;
rhsPhi133 = t1rhsPhi133 + t2rhsPhi + t3rhsPhi133;
                                                
rhsPhi200 = t1rhsPhi200 + t2rhsPhi + t3rhsPhi200;
rhsPhi201 = t1rhsPhi201 + t2rhsPhi + t3rhsPhi201;
rhsPhi202 = t1rhsPhi202 + t2rhsPhi + t3rhsPhi202;
rhsPhi203 = t1rhsPhi203 + t2rhsPhi + t3rhsPhi203;
rhsPhi211 = t1rhsPhi211 + t2rhsPhi + t3rhsPhi211;
rhsPhi212 = t1rhsPhi212 + t2rhsPhi + t3rhsPhi212;
rhsPhi213 = t1rhsPhi213 + t2rhsPhi + t3rhsPhi213;
rhsPhi222 = t1rhsPhi222 + t2rhsPhi + t3rhsPhi222;
rhsPhi223 = t1rhsPhi223 + t2rhsPhi + t3rhsPhi223;
rhsPhi233 = t1rhsPhi233 + t2rhsPhi + t3rhsPhi233;
                                                
rhsPhi300 = t1rhsPhi300 + t2rhsPhi + t3rhsPhi300;
rhsPhi301 = t1rhsPhi301 + t2rhsPhi + t3rhsPhi301;
rhsPhi302 = t1rhsPhi302 + t2rhsPhi + t3rhsPhi302;
rhsPhi303 = t1rhsPhi303 + t2rhsPhi + t3rhsPhi303;
rhsPhi311 = t1rhsPhi311 + t2rhsPhi + t3rhsPhi311;
rhsPhi312 = t1rhsPhi312 + t2rhsPhi + t3rhsPhi312;
rhsPhi313 = t1rhsPhi313 + t2rhsPhi + t3rhsPhi313;
rhsPhi322 = t1rhsPhi322 + t2rhsPhi + t3rhsPhi322;
rhsPhi323 = t1rhsPhi323 + t2rhsPhi + t3rhsPhi323;
rhsPhi333 = t1rhsPhi333 + t2rhsPhi + t3rhsPhi333;
