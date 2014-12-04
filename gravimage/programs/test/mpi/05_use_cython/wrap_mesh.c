#include <stdio.h>

typedef int bool;
enum { false, true };

void c_mesh_exp(double *r_min, double *r_max, int *N, double *mesh, bool *pp, char msg[]);

int main() {
  int i,N=5;
  double r_min = 0.0;
  double r_max = 1.0;
  double mesh[N];
  bool pp = true;
  char msg[255] = "some";
  printf("%s",msg);



  c_mesh_exp(&r_min, &r_max, &N, mesh, &pp, msg);
  for(i=0; i<N; i++){
    printf("%f\n",mesh[i]);
  }
}
