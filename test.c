#include <stdio.h>



static void
list2tensor(float *list, float *tensor, int dim)
{
    int i, j;
    float *p = tensor;
    for(i=0;i<dim;i++){
        *p = *list;
        for(j=1;j<dim-i;j++){
            *(++tensor) = *(++list);
            p += dim;
            *p = *list;
        }
        list++;
        tensor += i+2;
        p -= dim*(dim-i-2) - 1; 
    }
}


int main()
{
    float  a[4][4];
    float  b[10] = {0.0, 0.1, 0.2, 0.3,
                         1.1, 1.2, 1.3,
                              2.2, 2.3,
                                   3.3}; 
    list2tensor(b, a[0], 4); 
    
    for(int i=0;i<4;i++){
        printf("{%f, %f, %f, %f}\n",a[i][0], a[i][1],a[i][2],a[i][3]);
    }

    return 0;
}
