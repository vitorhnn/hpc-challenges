#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STABILITY 1.0f / sqrt(3.0f)

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

typedef struct {
  size_t yscale;
  size_t xscale;
  double *elements;
} matrix_t;

void mdf_heat(matrix_t *u0, matrix_t *u1,
              const unsigned int npX, const unsigned int npY,
              const unsigned int npZ, const double deltaH, const double deltaT,
              const unsigned int tsteps, const double boundaries);

void save2Text(matrix_t *u, const unsigned int npX, const unsigned int npY,
               const unsigned int npZ);

void save2Bin(matrix_t *u, const unsigned int npX, const unsigned int npY,
              const unsigned int npZ);


static void mset(matrix_t *mat, size_t x, size_t y, size_t z, double value)
{
  size_t pos;

  pos = x * mat->xscale + y * mat->yscale + z;

  mat->elements[pos] = value;
}

static double mget(matrix_t *mat, size_t x, size_t y, size_t z)
{
  size_t pos;

  pos = x * mat->xscale + y * mat->yscale + z;

  return mat->elements[pos];
}

static matrix_t *matrix_new(size_t xdim, size_t ydim, size_t zdim)
{
    matrix_t *mat = malloc(sizeof(matrix_t));
    if (!mat) {
        return NULL;
    }

    mat->elements = malloc(sizeof(double) * xdim * ydim * zdim);
    if (!mat->elements) {
        free(mat);
        return NULL;
    }

    mat->xscale = ydim * zdim;
    mat->yscale = zdim;

    return mat;
}

int main(int ac, char **av) {
  matrix_t *u0;
  matrix_t *u1;
  double deltaT = 0.01;
  double deltaH = 0.25f;
  double time = atof(av[1]);
  double sizeX = atof(av[2]); // 1.0f;
  double sizeY = atof(av[3]); // 1.0f;
  double sizeZ = atof(av[4]); // 1.0f;
  double temp = atof(av[5]);
  int flag2save = atoi(av[6]);

  unsigned int npX =
      (unsigned int)(sizeX / deltaH); // Number of points in X axis
  unsigned int npY = (unsigned int)(sizeY / deltaH);
  unsigned int npZ = (unsigned int)(sizeZ / deltaH);
  unsigned int tsteps = (unsigned int)(time / deltaT);
  fprintf(stdout, "\nSimulação - Domínio(x = %u, y = %u, z = %u, t = %u\n", npX,
          npY, npZ, tsteps);
  fprintf(stdout, "Dt(%f) Dh(%f)\n", deltaT, deltaH);
  // printf("p(%u, %u, %u)\n", npX, npY, npZ);
  // Allocing memory
  u0 = matrix_new(npX, npY, npZ);
  u1 = matrix_new(npX, npY, npZ);

  for (unsigned int i = 0; i < npZ; i++) {
    for (unsigned int j = 0; j < npY; j++) {
      for (unsigned int k = 0; k < npZ; k++) {
        mset(u0, i, j, k, 0.0);
        mset(u1, i, j, k, 0.0);
      }
    }
  }

  mdf_heat(u0, u1, npX, npY, npZ, deltaH, deltaT, tsteps, temp);

  if (flag2save == 1) {
    save2Bin(u0, npX, npY, npZ);
  }

  return EXIT_SUCCESS;
}

void mdf_heat(matrix_t *u0, matrix_t *u1,
              const unsigned int npX, const unsigned int npY,
              const unsigned int npZ, const double deltaH, const double deltaT,
              const unsigned int tsteps, const double boundaries) {

  const double alpha = deltaT / (deltaH * deltaH);
  assert(alpha < STABILITY);

  unsigned int step = tsteps / 20;

  #pragma omp parallel
  for (unsigned int steps = 0; steps < tsteps; steps++) {
    #pragma omp for nowait
    for (unsigned int i = 0; i < npZ; i++) {
      for (unsigned int j = 0; j < npY; j++) {
        for (unsigned int k = 0; k < npX; k++) {
          double left = boundaries;
          double right = boundaries;
          double up = boundaries;
          double down = boundaries;
          double top = boundaries;
          double bottom = boundaries;

          if (likely((k > 0) && (k < (npX - 1)))) {
            left = mget(u0, i, j, k - 1);
            right = mget(u0, i, j, k + 1);
          } else if (k == 0)
            right = mget(u0, i, j, k + 1);
          else
            left = mget(u0, i, j, k - 1);

          if (likely((j > 0) && (j < (npY - 1)))) {
            up = mget(u0, i, j - 1, k);
            down = mget(u0, i, j + 1, k);
          } else if (j == 0)
            down = mget(u0, i, j + 1, k);
          else
            up = mget(u0, i, j - 1, k);

          if (likely((i > 0) && (i < (npZ - 1)))) {
            top = mget(u0, i - 1, j, k);
            bottom = mget(u0, i + 1, j, k);
          } else if (i == 0)
            bottom = mget(u0, i + 1, j, k);
          else
            top = mget(u0, i - 1, j, k);

          mset(u1, i, j, k, alpha * (top + bottom + up + down + left + right -
                                 (6.0f * mget(u0, i, j, k))) +
                        mget(u0, i, j, k));
        }
      }
    }

    #pragma omp single
    {
      matrix_t *ptr = u0;
      u0 = u1;
      u1 = ptr;

      /*
      if ((steps % step) == 0) {
        fprintf(stdout, ".");
        fflush(stdout);
      }
      */
    }
  }
}
/*
 * Salva a saída para texto. Valores tendem a ficar próximos de 100
 */
void save2Text(matrix_t *u, const unsigned int npX, const unsigned int npY,
               const unsigned int npZ) {
  FILE *ptr = fopen("mdf.txt", "w+");
  fprintf(stdout, "\nSaving mdf.txt");
  fflush(stdout);
  assert(ptr != NULL);

  for (unsigned int i = 0; i < npZ; i++) {
    for (unsigned int j = 0; j < npY; j++) {
      for (unsigned int k = 0; k < npX; k++) {
        fprintf(ptr, "%u %u %u %lf \n", k, j, i, mget(u, i, j, k));
      }
    }
  }

  fprintf(stdout, "\t[OK]");
  fclose(ptr);
}

/*
 * Salva a saída em binário. Uso do comando diff para verificar se a saída está
 * ok
 */

void save2Bin(matrix_t *u, const unsigned int npX, const unsigned int npY,
              const unsigned int npZ) {
  char fileName[256];
  sprintf(fileName, "%s.bin", __FILE__);
  FILE *ptr = fopen(fileName, "w+");
  fprintf(stdout, "\nSaving %s", fileName);
  fflush(stdout);
  assert(ptr != NULL);

  for (unsigned int i = 0; i < npZ; i++) {
    for (unsigned int j = 0; j < npY; j++) {
      fwrite(&u->elements[i * u->xscale + j * u->yscale], sizeof(double), npX, ptr);
    }
  }

  fprintf(stdout, "\t[OK]");
  fclose(ptr);
}
