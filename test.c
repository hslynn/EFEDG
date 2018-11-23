#include <stdio.h>

int main()
{
    int b[3] = {0, 1, 2};
    
    int *a0 = b, *a1 = b + 1, *a2 = b + 2;
    int *a[3] = {a0, a1, a2};
    
    *a[1] = 3;
    
    printf("%d\n", b[1]);
    return 0;
}
