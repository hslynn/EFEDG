void
phgQuadDofGradBas(ELEMENT *e, DOF *p, DOF *v, int m, int order, FLOAT *values)
{
    int i;
    const FLOAT *g1, *g2, *w;
    FLOAT d1, d2, d3, vol;
    QUAD *quad;

    assert(!SpecialDofType(v->type));
    assert(p->dim == 1);

    if (order < 0) 
	    order = DofTypeOrder(p, e);
	 
    quad = phgQuadGetQuad3D(order);
    g1 = phgQuadGetDofValues(e, p, quad);
 
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;
    
    
    d1 = d2 = d3 = 0.;
    for (i = 0; i < quad->npoints; i++) {
        d1 += *(g1) * (*g2++) * (*w);
        d2 += *(g1) * (*g2++) * (*w);
        d3 += *(g1) * (*g2++) * (*w);
        g1++;
        w++;
    }

    vol = phgGeomGetVolume(p->g, e);
    values[0] = d1 * vol;
    values[1] = d2 * vol;
    values[2] = d3 * vol;
}

FLOAT
dgQuadFaceNeighDofDotBas(ELEMENT *e, DOF *u, int face, int N, ELEMENT *neigh_e, 
DOF *v, int neigh_face, int order)
{
    int i, j, e_v[NVert], neigh_v[NVert];
    FLOAT d, lambda[Dim + 1], neigh_lambda[Dim + 1];
    FLOAT *dof, *buffer;
    const FLOAT *bas, *p, *w;
    QUAD *quad;
    DOF_TYPE *type;

    assert(!SpecialDofType(u->type));
    assert(face >= 0 && face <= 3);

    type = (DofIsHP(u) ? u->hp->info->types[u->hp->max_order] : u->type);

    if (order < 0) {
	    i = DofTypeOrder(v, neigh_e);
	    j = DofTypeOrder(u, e);
	    if (i < 0)
	        i = j;
	    order = i + j;
    }
    quad = phgQuadGetQuad2D(order);

    /*有更好的实现方法吗？*/
    /*e_v[i]: face号面上的i号顶点在本单元的编号；
      neigh_v[i]: 该顶点在face面的相邻单元的编号*/
    for (i=0; i < NVert - 1; i++) {
        e_v[i] = GetFaceVertex(face, i);
        for (j=0; j < NVert; j++) {
            if (neigh_e->verts[j] == e->verts[e_v[i]]) {
                neigh_v[i] = j;
                break; 
            } 
        }
    }
    lambda[face] = 0.;
    neigh_lambda[neigh_face] = 0.;

    buffer = phgAlloc(sizeof(*buffer));
    p = quad->points;
    w = quad->weights;

    d = 0.;
    for (i = 0; i < quad->npoints; i++) {
        for (j = 0; j < NVert - 1; j++) {
        lambda[e_v[j]] = *(p++);
        neigh_lambda[neigh_v[j]] = lambda[e_v[j]];
        }
        
        dof = phgDofEval(v, neigh_e, neigh_lambda, buffer);
        bas = (FLOAT *)type->BasFuncs(u, e, N, N + 1, lambda);
        d += *(bas++) * *(dof++) * *(w++);
    } 
    phgFree(buffer);

    return d * phgGeomGetFaceArea(u->g, e, face);
}
