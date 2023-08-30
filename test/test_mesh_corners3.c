#include<p8est_mesh.h>
#include<p8est_extended.h>
#include<p8est_iterate.h>
 
#define FALSE 0
#define TRUE 1
 
static
int refine(p8est_t* p4est, p4est_topidx_t which_tree, p8est_quadrant_t* quadrant)
{
    int refine = TRUE;
    if(quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0)
    {
        refine = FALSE;
    }
    return refine;
}
 
static
int check_corner(p4est_topidx_t qid,
                 p4est_locidx_t c,
                 p4est_locidx_t qtc,
                 p4est_locidx_t expected_qtc)
{
    if(qtc != expected_qtc)
    {
        printf("ERROR: qid %d corner %d should have neighbor %d, but has %d\n",
               qid, c, expected_qtc, qtc);
        return 1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    int num_errors = 0;

    sc_MPI_Comm mpicomm = sc_MPI_COMM_SELF;

    sc_MPI_Init(&argc,&argv);
    sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    p8est_connectivity_t* conn = p8est_connectivity_new_unitcube();
    p8est_t* p8est = p8est_new_ext(mpicomm, conn, 0, 1, 0, 0, NULL, NULL);

    // refine all quadrents to level 2 except quadrent touching 0,0,0
    p8est_refine_ext(p8est, FALSE, 2, refine, NULL, NULL);

    p8est_ghost_t* ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);

    p8est_mesh_t* mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FULL);
    P4EST_ASSERT(mesh->local_num_quadrants == 57);

    // assuming that the coarse quadrent is zero with ordering
    p4est_locidx_t coarse_qid = 0;

    // upper x,y,z index
    p4est_locidx_t ux_idx = P8EST_FACES * coarse_qid + 1;
    p4est_locidx_t uy_idx = P8EST_FACES * coarse_qid + 3;
    p4est_locidx_t uz_idx = P8EST_FACES * coarse_qid + 5;

    P4EST_ASSERT(mesh->quad_to_face[ux_idx] < 0);
    P4EST_ASSERT(mesh->quad_to_face[uy_idx] < 0);
    P4EST_ASSERT(mesh->quad_to_face[uz_idx] < 0);

    // get the quadrents on upper x,y,z faces of coarse quadrent
    p4est_locidx_t qtq_ux, qtq_uy, qtq_uz;
    qtq_ux = mesh->quad_to_quad[ux_idx];
    qtq_uy = mesh->quad_to_quad[uy_idx];
    qtq_uz = mesh->quad_to_quad[uz_idx];

    p4est_locidx_t *qth_ux, *qth_uy, *qth_uz;
    qth_ux = (p4est_locidx_t*) sc_array_index(mesh->quad_to_half,qtq_ux);
    qth_uy = (p4est_locidx_t*) sc_array_index(mesh->quad_to_half,qtq_uy);
    qth_uz = (p4est_locidx_t*) sc_array_index(mesh->quad_to_half,qtq_uz);

    p4est_locidx_t qid1, qid2;
    p4est_locidx_t c1, c2;
    p4est_locidx_t qtc1, qtc2;

    //six pairs of corners with missing information

    //two along upper x upper y edge of coarse quad

    qid1 = qth_ux[1]; //lower quad along edge on x face of coarse quad
    qid2 = qth_uy[3]; //upper quad along edge on y face of coarse quad
    c1 = 6; // qid1 should have neighbor on uz uy lx corner
    c2 = 1; // qid2 should have neighbor on lz ly ux corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);

    qid1 = qth_ux[3]; //upper quad along edge on x face of coarse quad
    qid2 = qth_uy[1]; //lower quad along edge on y face of coarse quad
    c1 = 2; // qid1 should have neighbor on lz uy lx corner
    c2 = 5; // qid2 should have neighbor on uz ly ux corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);


    //two along upper x upper z edge of coarse quad

    qid1 = qth_ux[2]; //lower quad along edge on x face of coarse quad
    qid2 = qth_uz[3]; //upper quad along edge on z face of coarse quad
    c1 = 6; // qid1 should have neighbor on uz uy lx corner
    c2 = 1; // qid2 should have neighbor on lz ly ux corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);

    qid1 = qth_ux[3]; //upper quad along edge on x face
    qid2 = qth_uz[1]; //lower quad along edge on z face
    c1 = 4; // qid1 should have neighbor on uz ly lx corner
    c2 = 3; // qid2 should have neighbor on lz uy ux corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);

    //two along upper y upper z edge of coarse quad

    qid1 = qth_uy[2]; //lower quad along edge on y face of coarse quad
    qid2 = qth_uz[3]; //upper quad along edge on z face of coarse quad
    c1 = 5; // qid1 should have neighbor on uz ly ux corner
    c2 = 2; // qid2 should have neighbor on lz uy lx corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);

    qid1 = qth_uy[3]; //upper quad along edge on y face
    qid2 = qth_uz[2]; //lower quad along edge on z face
    c1 = 4; // qid1 should have neighbor on uz ly lx corner
    c2 = 3; // qid2 should have neighbor on lz uy ux corner
    qtc1 = mesh->quad_to_corner[P8EST_CHILDREN*qid1+c1];
    qtc2 = mesh->quad_to_corner[P8EST_CHILDREN*qid2+c2];
    num_errors += check_corner(qid1,c1,qtc1,qid2);
    num_errors += check_corner(qid2,c2,qtc2,qid1);

    p8est_mesh_destroy(mesh);
    p8est_ghost_destroy(ghost);
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);

    sc_finalize();
    sc_MPI_Finalize();

    return num_errors;
}