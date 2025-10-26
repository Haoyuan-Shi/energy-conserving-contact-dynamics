// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Haoyuan Shi (haoyuan.shi@pnnl.gov)
   Shi, Haoyuan, Christopher J. Mundy, Gregory K. Schenter, and Jaehun Chun. 
   Energy-Conserving Contact Dynamics of Nonspherical Rigid-Body Particles. 
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   Ref: Fraige, Langston, Matchett and Dodds, Particuology 2008, 6:455-466
   Note: The current implementation has not taken into account
         the contact history for friction forces.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polygon.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "body_rounded_polygon.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000
#define EPSILON 1e-3    // dimensionless threshold (dot products, end point checks, contact checks)
#define MAX_CONTACTS 4  // maximum number of contacts for 2D models
#define EFF_CONTACTS 2  // effective contacts for 2D models

//#define _CONVEX_POLYGON
//#define _POLYGON_DEBUG

enum {INVALID=0,NONE=1,VERTEXI=2,VERTEXJ=3,EDGE=4,OPP=5};

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygon::PairBodyRoundedPolygon(LAMMPS *lmp) : Pair(lmp)
{
  dmax = nmax = 0;
  discrete = nullptr;
  dnum = dfirst = nullptr;

  edmax = ednummax = 0;
  edge = nullptr;
  ednum = edfirst = nullptr;

  enclosing_radius = nullptr;
  rounded_radius = nullptr;
  maxerad = nullptr;

  single_enable = 0;
  restartinfo = 0;

  c_n = 0.1;
  c_t = 0.2;
  mu = 0.0;
  delta_ua = 1.0;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolygon::~PairBodyRoundedPolygon()
{
  memory->destroy(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  memory->destroy(edge);
  memory->destroy(ednum);
  memory->destroy(edfirst);

  memory->destroy(enclosing_radius);
  memory->destroy(rounded_radius);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(k_n);
    memory->destroy(k_na);
    memory->destroy(maxerad);
  }
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  int nei,nej,iefirst,jefirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl;
  double rsq,r,radi,radj,k_nij,k_naij;
  double facc[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  double *radius = atom->radius;
  tagint* tag = atom->tag;
  int *body = atom->body;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // grow the per-atom lists if necessary and initialize

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nall; i++)
    dnum[i] = ednum[i] = 0;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      npi = dnum[i];
      ifirst = dfirst[i];
      nei = ednum[i];
      iefirst = edfirst[i];
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      radj = radius[j];

      // body/body interactions

      evdwl = 0.0;
      facc[0] = facc[1] = facc[2] = 0;

      if (body[i] < 0 || body[j] < 0) continue;

      if (dnum[j] == 0) body2space(j);
      npj = dnum[j];
      jfirst = dfirst[j];
      nej = ednum[j];
      jefirst = edfirst[j];

      k_nij = k_n[itype][jtype];
      k_naij = k_na[itype][jtype];

      // no interaction

      r = sqrt(rsq);
      if (r > radi + radj + cut_inner) continue;

      if (npi == 1 && npj == 1) {
        sphere_against_sphere(i, j, delx, dely, delz, rsq, k_nij, k_naij, x, v, f, evflag);
        continue;
      }

      // reset vertex and edge forces

      for (ni = 0; ni < npi; ni++) {
        discrete[ifirst+ni][3] = 0;
        discrete[ifirst+ni][4] = 0;
        discrete[ifirst+ni][5] = 0;
      }

      for (nj = 0; nj < npj; nj++) {
        discrete[jfirst+nj][3] = 0;
        discrete[jfirst+nj][4] = 0;
        discrete[jfirst+nj][5] = 0;
      }

      for (ni = 0; ni < nei; ni++) {
        edge[iefirst+ni][2] = 0;
        edge[iefirst+ni][3] = 0;
        edge[iefirst+ni][4] = 0;
      }

      for (nj = 0; nj < nej; nj++) {
        edge[jefirst+nj][2] = 0;
        edge[jefirst+nj][3] = 0;
        edge[jefirst+nj][4] = 0;
      }

      // check interaction between i's vertices and j' edges

      vertex_against_edge(i, j, k_nij, k_naij, x, v, angmom, f, torque, tag, evdwl, facc);

      // check interaction between j's vertices and i' edges

      vertex_against_edge(j, i, k_nij, k_naij, x, v, angmom, f, torque, tag, evdwl, facc);

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                               facc[0],facc[1],facc[2],delx,dely,delz);

    } // end for jj
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(k_n,n+1,n+1,"pair:k_n");
  memory->create(k_na,n+1,n+1,"pair:k_na");
  memory->create(maxerad,n+1,"pair:maxerad");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::settings(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal pair_style command");

  c_n = utils::numeric(FLERR,arg[0],false,lmp);
  c_t = utils::numeric(FLERR,arg[1],false,lmp);
  mu = utils::numeric(FLERR,arg[2],false,lmp);
  delta_ua = utils::numeric(FLERR,arg[3],false,lmp);
  cut_inner = utils::numeric(FLERR,arg[4],false,lmp);

  if (delta_ua < 0) delta_ua = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double k_n_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k_na_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_n[i][j] = k_n_one;
      k_na[i][j] = k_na_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::init_style()
{
  avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  if (!avec)
    error->all(FLERR,"Pair body/rounded/polygon requires atom style body");
  if (strcmp(avec->bptr->style,"rounded/polygon") != 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "body style rounded/polygon");
  bptr = dynamic_cast<BodyRoundedPolygon *>(avec->bptr);

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style body/rounded/polygon requires "
               "newton pair on");

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair body/rounded/polygon requires "
               "ghost atoms store velocity");

  neighbor->add_request(this);

  // find the maximum enclosing radius for each atom type

  int i, itype;
  double eradi;
  int* body = atom->body;
  int* type = atom->type;
  int ntypes = atom->ntypes;
  int nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->destroy(ednum);
    memory->destroy(edfirst);
    memory->destroy(enclosing_radius);
    memory->destroy(rounded_radius);
    nmax = atom->nmax;
    memory->create(dnum,nmax,"pair:dnum");
    memory->create(dfirst,nmax,"pair:dfirst");
    memory->create(ednum,nmax,"pair:ednum");
    memory->create(edfirst,nmax,"pair:edfirst");
    memory->create(enclosing_radius,nmax,"pair:enclosing_radius");
    memory->create(rounded_radius,nmax,"pair:rounded_radius");
  }

  ndiscrete = nedge = 0;
  for (i = 0; i < nlocal; i++)
    dnum[i] = ednum[i] = 0;

  double *merad = nullptr;
  memory->create(merad,ntypes+1,"pair:merad");
  for (i = 1; i <= ntypes; i++)
    maxerad[i] = merad[i] = 0;

  int ipour;
  for (ipour = 0; ipour < modify->nfix; ipour++)
    if (strcmp(modify->fix[ipour]->style,"pour") == 0) break;
  if (ipour == modify->nfix) ipour = -1;

  int idep;
  for (idep = 0; idep < modify->nfix; idep++)
    if (strcmp(modify->fix[idep]->style,"deposit") == 0) break;
  if (idep == modify->nfix) idep = -1;

  for (i = 1; i <= ntypes; i++) {
    merad[i] = 0.0;
    if (ipour >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[ipour]->extract("radius",itype));
    }
    if (idep >= 0) {
      itype = i;
      merad[i] =
        *((double *) modify->fix[idep]->extract("radius",itype));
    }
  }

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if (body[i] >= 0) {
      if (dnum[i] == 0) body2space(i);
      eradi = enclosing_radius[i];
      if (eradi > merad[itype]) merad[itype] = eradi;
    } else
      merad[itype] = 0;
  }

  MPI_Allreduce(&merad[1],&maxerad[1],ntypes,MPI_DOUBLE,MPI_MAX,world);

  memory->destroy(merad);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBodyRoundedPolygon::init_one(int i, int j)
{
  k_n[j][i] = k_n[i][j];
  k_na[j][i] = k_na[i][j];

  return (maxerad[i]+maxerad[j]);
}

/* ----------------------------------------------------------------------
   convert N sub-particles in body I to space frame using current quaternion
   store sub-particle space-frame displacements from COM in discrete list
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::body2space(int i)
{
  int ibonus = atom->body[i];
  AtomVecBody::Bonus *bonus = &avec->bonus[ibonus];
  int nsub = bptr->nsub(bonus);
  double *coords = bptr->coords(bonus);
  int body_num_edges = bptr->nedges(bonus);
  double* edge_ends = bptr->edges(bonus);
  double eradius = bptr->enclosing_radius(bonus);
  double rradius = bptr->rounded_radius(bonus);

  // get the number of sub-particles (vertices)
  // and the index of the first vertex of my body in the list

  dnum[i] = nsub;
  dfirst[i] = ndiscrete;

  // grow the vertex list if necessary
  // the first 3 columns are for coords, the last 3 for forces

  if (ndiscrete + nsub > dmax) {
    dmax += DELTA;
    memory->grow(discrete,dmax,6,"pair:discrete");
  }

  double p[3][3];
  MathExtra::quat_to_mat(bonus->quat,p);

  for (int m = 0; m < nsub; m++) {
    MathExtra::matvec(p,&coords[3*m],discrete[ndiscrete]);
    discrete[ndiscrete][3] = 0;
    discrete[ndiscrete][4] = 0;
    discrete[ndiscrete][5] = 0;
    ndiscrete++;
  }

  // get the number of edges (vertices)
  // and the index of the first edge of my body in the list

  ednum[i] = body_num_edges;
  edfirst[i] = nedge;

  // grow the edge list if necessary
  // the first 2 columns are for vertex indices within body, the last 3 for forces

  if (nedge + body_num_edges > edmax) {
    edmax += DELTA;
    memory->grow(edge,edmax,5,"pair:edge");
  }

  if ((body_num_edges > 0) && (edge_ends == nullptr))
    error->one(FLERR,"Inconsistent edge data for body of atom {}",
                                 atom->tag[i]);

  for (int m = 0; m < body_num_edges; m++) {
    edge[nedge][0] = static_cast<int>(edge_ends[2*m+0]);
    edge[nedge][1] = static_cast<int>(edge_ends[2*m+1]);
    edge[nedge][2] = 0;
    edge[nedge][3] = 0;
    edge[nedge][4] = 0;
    nedge++;
  }

  enclosing_radius[i] = eradius;
  rounded_radius[i] = rradius;
}

/* ----------------------------------------------------------------------
   Interaction between two spheres with different radii
   according to the 2D model from Fraige et al.
---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::sphere_against_sphere(int i, int j,
                       double delx, double dely, double delz, double rsq,
                       double k_n, double k_na, double** /*x*/, double** v,
                       double** f, int evflag)
{
  double rradi,rradj;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double rij,rsqinv,R,fx,fy,fz,fn[3],ft[3],fpair,shift,energy;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  rradi = rounded_radius[i];
  rradj = rounded_radius[j];

  rsqinv = 1.0/rsq;
  rij = sqrt(rsq);
  R = rij - (rradi + rradj);
  shift = k_na * cut_inner;

  double r_Fmin = shift / (k_na + k_n); 
  double E0 = -0.5 * shift * r_Fmin;
  double E1 = -0.5 * k_n * cut_inner * r_Fmin;
  energy = 0;

  if (R <= 0) {           // deformation occurs
    fpair = -k_n * R;
    energy = (0.5 * k_n * R) * R + E1;
  } else if (R <= r_Fmin) {
    fpair = -k_n * R;
    energy = (0.5 * k_n * R) * R + E1;
  } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
    fpair = k_na * R - shift;
    energy = (-0.5 * k_na * R + shift) * R + E0 + E1;
  } else fpair = 0.0;

  fx = delx*fpair/rij;
  fy = dely*fpair/rij;
  fz = delz*fpair/rij;

  double rmin = MIN(rradi, rradj);
  if (R <= EPSILON*rmin) { // in contact

    // relative translational velocity

    vr1 = v[i][0] - v[j][0];
    vr2 = v[i][1] - v[j][1];
    vr3 = v[i][2] - v[j][2];

    // normal component

    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // normal friction term at contact

    fn[0] = -c_n * vn1;
    fn[1] = -c_n * vn2;
    fn[2] = -c_n * vn3;

    // tangential friction term at contact,
    // excluding the tangential deformation term for now

    ft[0] = -c_t * vt1;
    ft[1] = -c_t * vt2;
    ft[2] = -c_t * vt3;

    fx += fn[0] + ft[0];
    fy += fn[1] + ft[1];
    fz += fn[2] + ft[2];
  }

  f[i][0] += fx;
  f[i][1] += fy;
  f[i][2] += fz;

  if (newton_pair || j < nlocal) {
    f[j][0] -= fx;
    f[j][1] -= fy;
    f[j][2] -= fz;
  }

  if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                           energy,0.0,fx,fy,fz,delx,dely,delz);
}

/* ----------------------------------------------------------------------
   Determine the interaction mode between i's vertices against j's edges

   i = atom i (body i)
   j = atom j (body j)
   x      = atoms' coordinates
   f      = atoms' forces
   torque = atoms' torques
   tag    = atoms' tags
   Return:
     interact = 0 no interaction at all
                1 there's at least one case where i's vertices interacts
                  with j's edges
---------------------------------------------------------------------- */

int PairBodyRoundedPolygon::vertex_against_edge(int i, int j,
                                                double k_n, double k_na,
                                                double** x, double** v,
                                                double** angmom, double** f,
                                                double** torque, tagint* tag,
                                                double &evdwl, double* facc)
{
  int ni, npi, ifirst;
  int nj, jfirst, nej, jefirst;
  double xpi[3], xpj[3], dist, eradj, rradi, rradj;
  double fx, fy, fz, energy;
  int interact;

  npi = dnum[i];
  ifirst = dfirst[i];
  rradi = rounded_radius[i];

  jfirst = dfirst[j];
  nej = ednum[j];
  jefirst = edfirst[j];
  eradj = enclosing_radius[j];
  rradj = rounded_radius[j];

  energy = 0;
  interact = 0;

  // loop through body i's vertices

  for (ni = 0; ni < npi; ni++) {

    // convert body-fixed coordinates to space-fixed, xi

    xpi[0] = x[i][0] + discrete[ifirst+ni][0];
    xpi[1] = x[i][1] + discrete[ifirst+ni][1];
    xpi[2] = x[i][2] + discrete[ifirst+ni][2];

    // compute the distance from the vertex to the COM of body j

    distance(xpi, x[j], dist);

    #ifdef _POLYGON_DEBUG
    printf("Distance between vertex %d of body %d (%0.1f %0.1f %0.1f) "
           "to body %d's COM: %f (cut = %0.1f)\n",
           ni, xpi[0], xpi[1], xpi[2], atom->tag[i], atom->tag[j], dist,
           eradj + rradi + rradj + cut_inner);
    #endif

    // the vertex is within the enclosing circle (sphere) of body j,
    // possibly interacting

    if (dist > eradj + rradj + rradi + cut_inner) continue;
    // printf("eradj: %f, rradj: %f, rradi: %f\n", eradj, rradj, rradi);

    int mode, contact, p2vertex;
    double d, R, hi[3], t, delx, dely, delz, fpair, shift, anglesum=0;
    double rij;
    double h_min[3], r_min=DELTA, DD = 0;

    // loop through body j's edges

    for (nj = 0; nj < nej; nj++) {

      // compute the distance between the edge nj to the vertex xpi

      mode = compute_distance_to_vertex(j, nj, x[j], rradj,
                                        xpi, rradi, cut_inner,
                                        d, hi, t, contact, anglesum);

      if (mode == INVALID || mode == NONE) continue;

      delx = xpi[0] - hi[0];
      dely = xpi[1] - hi[1];
      delz = xpi[2] - hi[2];
      d = sqrt(delx*delx + dely*dely + delz*delz);

      if (d < rradj + rradi + cut_inner) {
        if (d < r_min){
          DD = 1;
          r_min = d;
          // printf("%f\n",r_min);
          h_min[0] = hi[0];
          h_min[1] = hi[1];
          h_min[2] = hi[2];
          // printf("hmin: %f,%f,%f\n",h_min[1],h_min[2],h_min[3]);
        }
      }
    }

    if (fabs(anglesum - MY_2PI) < EPSILON) {
        
      //  penetration: vertex lies inside the shape

    } else if (DD == 1) {
        delx = xpi[0] - h_min[0];
        dely = xpi[1] - h_min[1];
        delz = xpi[2] - h_min[2];

        d = sqrt(delx*delx + dely*dely + delz*delz);

        R = d - (rradi + rradj);
        shift = k_na * cut_inner;

        double r_Fmin = shift / (k_na + k_n);
        double E0 = -0.5 * shift * r_Fmin;
        double E1 = -0.5 * k_n * cut_inner * r_Fmin;

        if (R <= 0) {           // deformation occurs
          fpair = -k_n * R;
          energy = (0.5 * k_n * R) * R + E1;
        } else if (R <= r_Fmin) {
          fpair = -k_n * R;
          energy = (0.5 * k_n * R) * R + E1;
        } else if (R <= cut_inner) {   // not deforming but cohesive ranges overlap
          fpair = k_na * R - shift;
          energy = (-0.5 * k_na * R + shift) * R + E0 + E1;
        } else fpair = 0.0;

        fx = delx*fpair/d;
        fy = dely*fpair/d;
        fz = delz*fpair/d;

        // accumulate force and torque to both bodies directly

        if (R <= 0) { // in contact
          contact_forces(xpi, i, ni, h_min, j, nj, fx, fy, fz, x, v, angmom, f, torque, evdwl, facc);
        }

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        sum_torque(x[i], xpi, fx, fy, fz, torque[i]);

        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
        sum_torque(x[j], h_min, -fx, -fy, -fz, torque[j]);

        facc[0] += fx; facc[1] += fy; facc[2] += fz;

        evdwl += energy;
      }

  } // end for looping through the vertices of body i

  return interact;
}

/* -------------------------------------------------------------------------
  Compute the distance between an edge of body i and a vertex from
  another body
  Input:
    ibody      = body i (i.e. atom i)
    edge_index = edge index of body i
    xmi        = atom i's coordinates (body i's center of mass)
    x0         = coordinate of the tested vertex from another body
    x0_rounded_radius = rounded radius of the tested vertex
    cut_inner  = cutoff for vertex-vertex and vertex-edge interaction
  Output:
    d          = Distance from a point x0 to an edge
    hi         = coordinates of the projection of x0 on the edge
    t          = ratio to determine the relative position of hi
                 wrt xi and xj on the segment
  contact      = 0 no contact between the queried vertex and the edge
                 1 contact detected
  return
    INVALID if the edge index is invalid
    NONE    if there is no interaction
    VERTEXI if the tested vertex interacts with the first vertex of the edge
    VERTEXJ if the tested vertex interacts with the second vertex of the edge
    EDGE    if the tested vertex interacts with the edge
------------------------------------------------------------------------- */

int PairBodyRoundedPolygon::compute_distance_to_vertex(int ibody,
                                                int edge_index,
                                                double *xmi,
                                                double rounded_radius,
                                                double* x0,
                                                double x0_rounded_radius,
                                                double cut_inner,
                                                double &d,
                                                double hi[3],
                                                double &t,
                                                int &contact,
                                                double &anglesum)
{
  if (edge_index >= ednum[ibody]) return INVALID;

  int mode,ifirst,iefirst,npi1,npi2;
  double xi1[3],xi2[3],u[3],v[3],uij[3];
  double u2[3],v2[3],costheta,magu2,magv2;
  double udotv, magv, magucostheta;
  double delx,dely,delz;
  double rmin = MIN(rounded_radius, x0_rounded_radius);

  ifirst = dfirst[ibody];
  iefirst = edfirst[ibody];
  npi1 = static_cast<int>(edge[iefirst+edge_index][0]);
  npi2 = static_cast<int>(edge[iefirst+edge_index][1]);

  // compute the space-fixed coordinates for the vertices of the edge

  xi1[0] = xmi[0] + discrete[ifirst+npi1][0];
  xi1[1] = xmi[1] + discrete[ifirst+npi1][1];
  xi1[2] = xmi[2] + discrete[ifirst+npi1][2];

  xi2[0] = xmi[0] + discrete[ifirst+npi2][0];
  xi2[1] = xmi[1] + discrete[ifirst+npi2][1];
  xi2[2] = xmi[2] + discrete[ifirst+npi2][2];

  // u = x0 - xi1

  u[0] = x0[0] - xi1[0];
  u[1] = x0[1] - xi1[1];
  u[2] = x0[2] - xi1[2];

  // v = xi2 - xi1

  v[0] = xi2[0] - xi1[0];
  v[1] = xi2[1] - xi1[1];
  v[2] = xi2[2] - xi1[2];

  // dot product between u and v = magu * magv * costheta

  udotv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  magv = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  magucostheta = udotv / magv;

  // uij is the unit vector pointing from xi to xj

  uij[0] = v[0] / magv;
  uij[1] = v[1] / magv;
  uij[2] = v[2] / magv;

  // position of the projection of x0 on the line (xi, xj)

  hi[0] = xi1[0] + magucostheta * uij[0];
  hi[1] = xi1[1] + magucostheta * uij[1];
  hi[2] = xi1[2] + magucostheta * uij[2];

  // distance from x0 to the line (xi, xj) = distance from x0 to hi

  distance(hi, x0, d);

  // determine the interaction mode
  // for 2D: a vertex can interact with one edge at most
  // for 3D: a vertex can interact with one face at most

  mode = EDGE;
  contact = 0;

  int m = opposite_sides(xi1, xi2, x0, xmi);

  MathExtra::sub3(xi1,x0,u2);
  MathExtra::sub3(xi2,x0,v2);
  magu2 = MathExtra::len3(u2);
  magv2 = MathExtra::len3(v2);
  costheta = MathExtra::dot3(u2,v2)/(magu2*magv2);
  anglesum += acos(costheta);
  

  int max_dim = 0;
  double max_diff = fabs(xi2[0] - xi1[0]);

  // Find dimension with largest difference
  for (int dim = 1; dim < 3; ++dim) {
    double diff = fabs(xi2[dim] - xi1[dim]);
    if (diff > max_diff) {
      max_diff = diff;
      max_dim = dim;
    }
  }

  // Avoid division by zero
  if (max_diff > EPSILON * rmin) {
    t = (hi[max_dim] - xi1[max_dim]) / (xi2[max_dim] - xi1[max_dim]);
  } else {
    error->all(FLERR, "project_pt_line: line segment is degenerate or too short");
  }

  if (t < 0) {
    hi[0] = xi1[0];
    hi[1] = xi1[1];
    hi[2] = xi1[2];
    // measure the distance from x0 to xi1
    delx = x0[0] - xi1[0];
    dely = x0[1] - xi1[1];
    delz = x0[2] - xi1[2];
    d = sqrt(delx*delx + dely*dely + delz*delz);
  } else if (t > 1) {
    // measure the distance from x0 to xi2
    hi[0] = xi2[0];
    hi[1] = xi2[1];
    hi[2] = xi2[2];
    delx = x0[0] - xi2[0];
    dely = x0[1] - xi2[1];
    delz = x0[2] - xi2[2];
    d = sqrt(delx*delx + dely*dely + delz*delz);
  }

  if (m == 0 && ednum[ibody] > 1) return mode = NONE;

  return mode;
}

/* ----------------------------------------------------------------------
  Compute contact forces between two bodies
  modify the force stored at the vertex and edge in contact by j_a
  sum forces and torque to the corresponding bodies
  fn = normal friction component
  ft = tangential friction component (-c_t * v_t)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::contact_forces(double *xpi, int ibody, int ni, double *xpj, int jbody, int nj, 
                                            double fx, double fy, double fz,
                                            double** x, double** v, double** angmom, double** f,
                                            double** torque, double &/*evdwl*/, double* facc)
{
  int ibonus,jbonus;
  double delx,dely,delz,rsq,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double fn[3],ft[3],vi[3],vj[3];
  double *quat, *inertia;
  AtomVecBody::Bonus *bonus;

  // compute the velocity of the vertex in the space-fixed frame

  ibonus = atom->body[ibody];
  bonus = &avec->bonus[ibonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xpi, x[ibody], v[ibody], angmom[ibody],
                 inertia, quat, vi);

  // compute the velocity of the point on the edge in the space-fixed frame

  jbonus = atom->body[jbody];
  bonus = &avec->bonus[jbonus];
  quat = bonus->quat;
  inertia = bonus->inertia;
  total_velocity(xpj, x[jbody], v[jbody], angmom[jbody],
                 inertia, quat, vj);

  // vector pointing from the vertex to the point on the edge

  delx = xpi[0] - xpj[0];
  dely = xpi[1] - xpj[1];
  delz = xpi[2] - xpj[2];
  rsq = delx*delx + dely*dely + delz*delz;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = vi[0] - vj[0];
  vr2 = vi[1] - vj[1];
  vr3 = vi[2] - vj[2];

  // normal component

  vnnr = vr1*delx + vr2*dely + vr3*delz;
  vn1 = delx*vnnr * rsqinv;
  vn2 = dely*vnnr * rsqinv;
  vn3 = delz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // normal friction term at contact

  fn[0] = -c_n * vn1;
  fn[1] = -c_n * vn2;
  fn[2] = -c_n * vn3;

  // tangential friction term at contact
  // excluding the tangential deformation term for now

  ft[0] = -c_t * vt1;
  ft[1] = -c_t * vt2;
  ft[2] = -c_t * vt3;

  // only the cohesive force is scaled by j_a

  // normal force magnitude
  double fnx = fx + fn[0];
  double fny = fy + fn[1];
  double fnz = fz + fn[2];

  double Fn_mag = sqrt(fnx*fnx + fny*fny + fnz*fnz);
  double Ft_mag = sqrt(ft[0]*ft[0] + ft[1]*ft[1] + ft[2]*ft[2]);

  // apply Coulomb limit
  if (Ft_mag > mu * Fn_mag) {
      double scale = (mu * Fn_mag) / Ft_mag;
      ft[0] *= scale;
      ft[1] *= scale;
      ft[2] *= scale;
  }

  fx = fn[0] + ft[0];
  fy = fn[1] + ft[1];
  fz = fn[2] + ft[2];

  f[ibody][0] += fx;
  f[ibody][1] += fy;
  f[ibody][2] += fz;
  sum_torque(x[ibody], xpi, fx, fy, fz, torque[ibody]);

  // accumulate forces to the vertex only

  facc[0] += fx; facc[1] += fy; facc[2] += fz;

  f[jbody][0] -= fx;
  f[jbody][1] -= fy;
  f[jbody][2] -= fz;
  sum_torque(x[jbody], xpj, -fx, -fy, -fz, torque[jbody]);
}

/* ----------------------------------------------------------------------
  Accumulate torque to body from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::sum_torque(double* xm, double *x, double fx,
                                        double fy, double fz, double* torque)
{
  double rx = x[0] - xm[0];
  double ry = x[1] - xm[1];
  double rz = x[2] - xm[2];
  double tx = ry * fz - rz * fy;
  double ty = rz * fx - rx * fz;
  double tz = rx * fy - ry * fx;
  torque[0] += tx;
  torque[1] += ty;
  torque[2] += tz;
}

/* ----------------------------------------------------------------------
  Test if two points a and b are in opposite sides of the line that
  connects two points x1 and x2
------------------------------------------------------------------------- */

int PairBodyRoundedPolygon::opposite_sides(double* x1, double* x2,
                                           double* a, double* b)
{
  double m_a = (x1[1] - x2[1])*(a[0] - x1[0]) + (x2[0] - x1[0])*(a[1] - x1[1]);
  double m_b = (x1[1] - x2[1])*(b[0] - x1[0]) + (x2[0] - x1[0])*(b[1] - x1[1]);
  // equal to zero when either a or b is inline with the line x1-x2
  if (m_a * m_b <= 0)
    return 1;
  else
    return 0;
}

/* ----------------------------------------------------------------------
  Calculate the total velocity of a point (vertex, a point on an edge):
    vi = vcm + omega ^ (p - xcm)
------------------------------------------------------------------------- */

void PairBodyRoundedPolygon::total_velocity(double* p, double *xcm,
                              double* vcm, double *angmom, double *inertia,
                              double *quat, double* vi)
{
  double r[3],omega[3],ex_space[3],ey_space[3],ez_space[3];
  r[0] = p[0] - xcm[0];
  r[1] = p[1] - xcm[1];
  r[2] = p[2] - xcm[2];
  MathExtra::q_to_exyz(quat,ex_space,ey_space,ez_space);
  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,
                             inertia,omega);
  vi[0] = omega[1]*r[2] - omega[2]*r[1] + vcm[0];
  vi[1] = omega[2]*r[0] - omega[0]*r[2] + vcm[1];
  vi[2] = omega[0]*r[1] - omega[1]*r[0] + vcm[2];
}

/* ---------------------------------------------------------------------- */

void PairBodyRoundedPolygon::distance(const double* x2, const double* x1,
                                      double& r)
{
  r = sqrt((x2[0] - x1[0]) * (x2[0] - x1[0])
    + (x2[1] - x1[1]) * (x2[1] - x1[1])
    + (x2[2] - x1[2]) * (x2[2] - x1[2]));
}