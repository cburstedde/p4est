#include<p8est_mesh.h>
#include<p8est_extended.h>
#include<p8est_iterate.h>

static int
refine (p8est_t *p4est, p4est_topidx_t which_tree,
        p8est_quadrant_t *quadrant)
{
  int                 refine = 1;
  if (quadrant->x == 0 && quadrant->y == 0 && quadrant->z == 0 &&
      which_tree == 3) {
    refine = 0;
  }
  return refine;
}

static int
check_corner (p4est_topidx_t qid, p4est_locidx_t c, p4est_locidx_t qtc,
              p4est_locidx_t expected_qtc)
{
  if (qtc != expected_qtc) {
    printf ("ERROR: qid %d corner %d should have neighbor %d, but has %d\n",
            qid, c, expected_qtc, qtc);
    return 1;
  }
  return 0;
}

static int
check_corner_across_trees (p4est_topidx_t qid, p4est_locidx_t c,
                           p4est_locidx_t qtc, p4est_locidx_t expected_qtc,
                           p8est_mesh_t *mesh)
{
  int                 offset;
  int                 cind;
  p4est_locidx_t      new_qtc;

  /* corner neighbor information for quadrants from different trees are stored
   * in the corner_quad array and have to be accessed using corner_offset */
  cind = qtc - mesh->local_num_quadrants - mesh->ghost_num_quadrants;
  P4EST_ASSERT (cind >= 0);
  offset = *(int *) sc_array_index (mesh->corner_offset, cind);
  new_qtc = *(p4est_locidx_t *) sc_array_index (mesh->corner_quad, offset);

  return check_corner (qid, c, new_qtc, expected_qtc);
}

int
main (int argc, char **argv)
{
  int                 num_errors = 0;
  sc_MPI_Comm         mpicomm = sc_MPI_COMM_SELF;
  p8est_connectivity_t *conn;
  p8est_t            *p8est;
  p8est_ghost_t      *ghost;
  p8est_mesh_params_t params;
  p8est_mesh_t       *mesh;
  p4est_locidx_t      coarse_qid;
  p4est_locidx_t      lx_idx, ux_idx, ly_idx, uy_idx, uz_idx;
  p4est_locidx_t      qtq_lx, qtq_ux, qtq_ly, qtq_uy, qtq_uz;
  p4est_locidx_t     *qth_lx, *qth_ux, *qth_ly, *qth_uy, *qth_uz;
  p4est_locidx_t      qid1, qid2;
  p4est_locidx_t      c1, c2;
  p4est_locidx_t      qtc1, qtc2;

  sc_MPI_Init (&argc, &argv);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create a 2x2 brick that is uniformly refined to level 2 except for a level 1
   * quadrant with coordinates (0,0,0) in tree 3. At the edges of the coarse
   * quadrant there are edgehanging-corner neighbors inside the tree as well as
   * across tree faces and edges. */
  conn = p8est_connectivity_new_brick (2, 2, 1, 0, 0, 0);
  p8est = p8est_new_ext (mpicomm, conn, 0, 1, 0, 0, NULL, NULL);

  /* refine all quadrants to level 2 except quadrant in tree 3 touching 0,0,0 */
  p8est_refine_ext (p8est, 0, 2, refine, NULL, NULL);

  ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);

  /* create mesh with edgehanging corner neighbor information */
  p8est_mesh_params_init (&params);
  params.btype = P8EST_CONNECT_FULL;
  params.edgehanging_corners = 1;
  mesh = p8est_mesh_new_params (p8est, ghost, &params);
  P4EST_ASSERT (mesh->local_num_quadrants == 249);

  /* assuming that the coarse quadrant is zero in tree 3g */
  coarse_qid = 192;

  /* upper and lower x,y,z indices */
  lx_idx = P8EST_FACES * coarse_qid + 0;
  ux_idx = P8EST_FACES * coarse_qid + 1;
  ly_idx = P8EST_FACES * coarse_qid + 2;
  uy_idx = P8EST_FACES * coarse_qid + 3;
  uz_idx = P8EST_FACES * coarse_qid + 5;

  P4EST_ASSERT (mesh->quad_to_face[lx_idx] < 0);
  P4EST_ASSERT (mesh->quad_to_face[ux_idx] < 0);
  P4EST_ASSERT (mesh->quad_to_face[ly_idx] < 0);
  P4EST_ASSERT (mesh->quad_to_face[uy_idx] < 0);
  P4EST_ASSERT (mesh->quad_to_face[uz_idx] < 0);

  /* get the quadrants on upper x,y,z faces of the coarse quadrant */
  qtq_lx = mesh->quad_to_quad[lx_idx];
  qtq_ux = mesh->quad_to_quad[ux_idx];
  qtq_ly = mesh->quad_to_quad[ly_idx];
  qtq_uy = mesh->quad_to_quad[uy_idx];
  qtq_uz = mesh->quad_to_quad[uz_idx];

  qth_lx = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq_lx);
  qth_ux = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq_ux);
  qth_ly = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq_ly);
  qth_uy = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq_uy);
  qth_uz = (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, qtq_uz);

  /*** test intra tree corners ***/
  /* six pairs of corners with missing information */

  /* two along upper x upper y edge of coarse quad */
  qid1 = qth_ux[1];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_uy[3];             /* upper quad along edge on y face of coarse quad */
  c1 = 6;                       /* qid1 should have neighbor on uz uy lx corner */
  c2 = 1;                       /* qid2 should have neighbor on lz ly ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  qid1 = qth_ux[3];             /* upper quad along edge on x face of coarse quad */
  qid2 = qth_uy[1];             /* lower quad along edge on y face of coarse quad */
  c1 = 2;                       /* qid1 should have neighbor on lz uy lx corner */
  c2 = 5;                       /* qid2 should have neighbor on uz ly ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  /* two along upper x upper z edge of coarse quad */
  qid1 = qth_ux[2];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_uz[3];             /* upper quad along edge on z face of coarse quad */
  c1 = 6;                       /* qid1 should have neighbor on uz uy lx corner */
  c2 = 1;                       /* qid2 should have neighbor on lz ly ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  qid1 = qth_ux[3];             /* upper quad along edge on x face */
  qid2 = qth_uz[1];             /* lower quad along edge on z face */
  c1 = 4;                       /* qid1 should have neighbor on uz ly lx corner */
  c2 = 3;                       /* qid2 should have neighbor on lz uy ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  /* two along upper y upper z edge of coarse quad */
  qid1 = qth_uy[2];             /* lower quad along edge on y face of coarse quad */
  qid2 = qth_uz[3];             /* upper quad along edge on z face of coarse quad */
  c1 = 5;                       /* qid1 should have neighbor on uz ly ux corner */
  c2 = 2;                       /* qid2 should have neighbor on lz uy lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  qid1 = qth_uy[3];             /* upper quad along edge on y face */
  qid2 = qth_uz[2];             /* lower quad along edge on z face */
  c1 = 4;                       /* qid1 should have neighbor on uz ly lx corner */
  c2 = 3;                       /* qid2 should have neighbor on lz uy ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner (qid1, c1, qtc1, qid2);
  num_errors += check_corner (qid2, c2, qtc2, qid1);

  /*** test inter-tree corners across tree-face ***/
  /* eight pairs of corners with missing information */

  /* two along lower x upper y edge of coarse quad (across face to tree 2) */
  qid1 = qth_lx[1];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_uy[2];             /* upper quad along edge on y face of coarse quad */
  c1 = 7;                       /* qid1 should have neighbor on uz uy ux corner */
  c2 = 0;                       /* qid2 should have neighbor on lz ly lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  qid1 = qth_lx[3];             /* upper quad along edge on x face of coarse quad */
  qid2 = qth_uy[0];             /* lower quad along edge on y face of coarse quad */
  c1 = 3;                       /* qid1 should have neighbor on lz uy ux corner */
  c2 = 4;                       /* qid2 should have neighbor on uz ly lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  /* two along lower x upper z edge of coarse quad (across face to tree 2) */
  qid1 = qth_lx[2];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_uz[2];             /* upper quad along edge on z face of coarse quad */
  c1 = 7;                       /* qid1 should have neighbor on uz uy ux corner */
  c2 = 0;                       /* qid2 should have neighbor on lz ly lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  qid1 = qth_lx[3];             /* upper quad along edge on x face of coarse quad */
  qid2 = qth_uz[0];             /* lower quad along edge on z face of coarse quad */
  c1 = 5;                       /* qid1 should have neighbor on uz ly ux corner */
  c2 = 2;                       /* qid2 should have neighbor on lz uy lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  /* two along upper x lower y edge of coarse quad (across face to tree 1) */
  qid1 = qth_ux[0];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_ly[3];             /* upper quad along edge on y face of coarse quad */
  c1 = 4;                       /* qid1 should have neighbor on uz ly lx corner */
  c2 = 3;                       /* qid2 should have neighbor on lz uy ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  qid1 = qth_ux[2];             /* upper quad along edge on x face of coarse quad */
  qid2 = qth_ly[1];             /* lower quad along edge on y face of coarse quad */
  c1 = 0;                       /* qid1 should have neighbor on lz ly lx corner */
  c2 = 7;                       /* qid2 should have neighbor on uz uy ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  /* two along lower y upper z edge of coarse quad (across face to tree 1) */
  qid1 = qth_ly[2];             /* lower quad along edge on y face of coarse quad */
  qid2 = qth_uz[1];             /* upper quad along edge on z face of coarse quad */
  c1 = 7;                       /* qid1 should have neighbor on uz uy ux corner */
  c2 = 0;                       /* qid2 should have neighbor on lz ly lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  qid1 = qth_ly[3];             /* upper quad along edge on y face of coarse quad */
  qid2 = qth_uz[0];             /* lower quad along edge on z face of coarse quad */
  c1 = 6;                       /* qid1 should have neighbor on uz uy lx corner */
  c2 = 1;                       /* qid2 should have neighbor on lz ly ux corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  /*** test inter-tree corners across tree-edge ***/
  /* two pairs of corners with missing information */

  /* two along lower x lower y edge of coarse quad (across edge to tree 0) */
  qid1 = qth_lx[0];             /* lower quad along edge on x face of coarse quad */
  qid2 = qth_ly[2];             /* upper quad along edge on y face of coarse quad */
  c1 = 5;                       /* qid1 should have neighbor on uz ly ux corner */
  c2 = 2;                       /* qid2 should have neighbor on lz uy lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  qid1 = qth_lx[2];             /* upper quad along edge on x face of coarse quad */
  qid2 = qth_ly[0];             /* lower quad along edge on y face of coarse quad */
  c1 = 1;                       /* qid1 should have neighbor on lz ly ux corner */
  c2 = 6;                       /* qid2 should have neighbor on uz uy lx corner */
  qtc1 = mesh->quad_to_corner[P8EST_CHILDREN * qid1 + c1];
  qtc2 = mesh->quad_to_corner[P8EST_CHILDREN * qid2 + c2];
  num_errors += check_corner_across_trees (qid1, c1, qtc1, qid2, mesh);
  num_errors += check_corner_across_trees (qid2, c2, qtc2, qid1, mesh);

  /* check for errors */
  SC_CHECK_ABORT (num_errors == 0, "Missing edge-hanging corner neighbors");

  /* cleanup */
  p8est_mesh_destroy (mesh);
  p8est_ghost_destroy (ghost);
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  sc_finalize ();
  sc_MPI_Finalize ();

  return num_errors;
}
