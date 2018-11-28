static void
list2tensor(FLOAT *list, FLOAT *tensor, int dim)
{
    else if(dim>1){
        for(i=0;i<dim;i++){
            tensor += i;
            *(tensor++) = *(list++);
        }
    }
}
