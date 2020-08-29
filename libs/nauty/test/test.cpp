
#include <vector>

extern "C" {
#include "../traces.h"
}

using namespace std;

int main() {

  DYNALLSTAT(int,lab1,lab1_sz);
  DYNALLSTAT(int,ptn,ptn_sz);
  DYNALLSTAT(int,orbits,orbits_sz);

  SG_DECL(sg1);
  SG_DECL(cg1);

  static DEFAULTOPTIONS_TRACES(options);
  /* Select option for canonical labelling */
  options.getcanon = TRUE;
  TracesStats stats;

  int nv = 3;

  //int edges = 2;

  int m = SETWORDSNEEDED(nv);
  // checks that this file is compiled compatibly with the given parameters
  nauty_check(WORDSIZE,m,nv,NAUTYVERSIONID);
  //setword  *workspace = new setword [10*m];

  // make the graph


  DYNALLOC1(int,lab1,lab1_sz,nv,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,nv,"malloc");
  DYNALLOC1(int,orbits,orbits_sz,nv,"malloc");

  //SG_ALLOC(sg1, verts, edges, "malloc");

  vector<size_t> v_edge_indices(nv);
  vector<int> d_out_degrees(nv);
  vector<int> e_neighbours(4);

  d_out_degrees[0] = 2;
  d_out_degrees[1] = 1;
  d_out_degrees[2] = 1;

  // can be computed from d_out_degrees
  v_edge_indices[0] = 0;
  v_edge_indices[1] = 2;
  v_edge_indices[1] = 3;

  e_neighbours[0] = 1;
  e_neighbours[1] = 2;
  e_neighbours[2] = 0;
  e_neighbours[3] = 0;

  sg1.nde = e_neighbours.size();
  sg1.nv = v_edge_indices.size();

  sg1.d = d_out_degrees.data();
  sg1.dlen = d_out_degrees.size();
  sg1.v = v_edge_indices.data();
  sg1.vlen = v_edge_indices.size();
  sg1.e = e_neighbours.data();
  sg1.elen = e_neighbours.size();


  Traces(&sg1,lab1,ptn,orbits,&options,&stats,&cg1);

  int i;
  for (i = 0; i < nv; ++i) printf(" %d->%d ",i,lab1[i]);
  printf("\n");

  //SG_FREE( sg1 ); - all arrays are allocated inside vector
  SG_FREE( cg1 );

  nausparse_freedyn();

  free(lab1);
  free(ptn);
  free(orbits);
}
