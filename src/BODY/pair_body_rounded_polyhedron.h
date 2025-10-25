/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(body/rounded/polyhedron,PairBodyRoundedPolyhedron);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBodyRoundedPolyhedron : public Pair {
 public:
  PairBodyRoundedPolyhedron(class LAMMPS *);
  ~PairBodyRoundedPolyhedron() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  virtual void kernel_force(double R, int itype, int jtype, double &energy, double &fpair);

 protected:
  double **k_n;        // normal repulsion strength
  double **k_na;       // normal attraction strength
  double c_n;          // normal damping coefficient
  double c_t;          // tangential damping coefficient
  double mu;           // normal friction coefficient during gross sliding
  double A_ua;         // characteristic contact area
  double cut_inner;    // cutoff for interaction between vertex-edge surfaces

  class AtomVecBody *avec;
  class BodyRoundedPolyhedron *bptr;

  double **discrete;    // list of all sub-particles for all bodies
  int ndiscrete;        // number of discretes in list
  int dmax;             // allocated size of discrete list
  int *dnum;            // number of discretes per line, 0 if uninit
  int *dfirst;          // index of first discrete per each line
  int nmax;             // allocated size of dnum,dfirst vectors

  double **edge;    // list of all edge for all bodies
  int nedge;        // number of edge in list
  int edmax;        // allocated size of edge list
  int *ednum;       // number of edges per line, 0 if uninit
  int *edfirst;     // index of first edge per each line
  int ednummax;     // allocated size of ednum,edfirst vectors

  double **face;    // list of all edge for all bodies
  int nface;        // number of faces in list
  int facmax;       // allocated size of face list
  int *facnum;      // number of faces per line, 0 if uninit
  int *facfirst;    // index of first face per each line
  int facnummax;    // allocated size of facnum,facfirst vectors

  double *enclosing_radius;    // enclosing radii for all bodies
  double *rounded_radius;      // rounded radii for all bodies
  double *maxerad;             // per-type maximum enclosing radius

  void allocate();
  void body2space(int);

  // sphere-sphere interaction
  void sphere_against_sphere(int ibody, int jbody, int itype, int jtype, double delx, double dely,
                             double delz, double rsq, double **v, double **f, int evflag);
  // sphere-face interaction
  void sphere_against_face(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                           double **f, double **torque, double **angmom, int evflag);
  // edge-edge interactions
  int edge_against_edge(int ibody, int jbody, int itype, int jtype, double **x, double &evdwl, double *facc);
  // edge-face interactions
  int vertex_against_face(int ibody, int jbody, int itype, int jtype, double **x, double &evdwl, double *facc);

  // a face vs. a single edge
  int interaction_face_to_vertex(int ibody, int face_index, double *xmi, double rounded_radius_i,
                                 int jbody, int vertex_index, double *xmj, double rounded_radius_j,
                                 int itype, int jtype, double cut_inner, double *e1_1, double *e1_2, double &r12,
                                 double &energy, double *facc);
  // a face vs. a single edge
  int interaction_edge_to_vertex(int ibody, double rounded_radius_i,
                                 int jbody, int vertex_index, double *xmj, double rounded_radius_j,
                                 int itype, int jtype, double cut_inner, double *e1_1, double *e1_2, double &r12, 
                                 double &energy, double *facc);
  // an edge vs. an edge from another body
  int interaction_edge_to_edge(int ibody, int edge_index_i, double *xmi, double rounded_radius_i,
                               int jbody, int edge_index_j, double *xmj, double rounded_radius_j,
                               int itype, int jtype, double cut_inner, double &energy, double *facc);

  // compute contact forces if contact points are detected
  void contact_forces(int ibody, int jbody, double *xi, double *xj, double delx, double dely,
                      double delz, double fx, double fy, double fz, double **x, double **v,
                      double **angmom, double **f, double **torque, double *facc);

  // compute force and torque between two bodies given a pair of interacting points
  void pair_force_and_torque(int ibody, int jbody, double *pi, double *pj, double r,
                             double contact_dist, int itype, int jtype, double **x, double **v,
                             double **f, double **torque, double **angmom, int jflag,
                             double &energy, double *facc);

  // accumulate torque to a body given a force at a given point
  void sum_torque(double *xm, double *x, double fx, double fy, double fz, double *torque);

  // helper functions
  int opposite_sides(double *n, double *x0, double *a, double *b);
  void project_pt_plane(const double *q, const double *p, const double *n, double *q_proj,
                        double &d);
  void project_pt_plane(const double *q, const double *x1, const double *x2, const double *x3,
                        double *q_proj, double &d, int &inside);
  void project_pt_line(const double *q, const double *xi1, const double *xi2, double *h, double &d,
                       double &t);
  void inside_polygon(int ibody, int face_index, double *xmi, const double *q1, const double *q2,
                      int &inside1, int &inside2);
  void inside_polygon_d(int ibody, int face_index, double* xmi, const double* q1, const double* q2,
                        int& inside1, int& inside2, double* e_out1_1, double* e_out1_2, double& r_out1,
                        double* e_out2_1, double* e_out2_2, double& r_out2);

  void distance_bt_edges(const double *x1, const double *x2, const double *x3, const double *x4,
                         double *h1, double *h2, double &t1, double &t2, double &r);
  void total_velocity(double *p, double *xcm, double *vcm, double *angmom, double *inertia,
                      double *quat, double *vi);
  void sanity_check();
};

}    // namespace LAMMPS_NS

#endif
#endif
