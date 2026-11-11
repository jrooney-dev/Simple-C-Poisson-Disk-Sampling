#ifndef POISSON_DISK_H

#define POISSON_DISK_H
#include <stdbool.h>
#include "pcg64.h"

struct GRID
{
    int width;
    int height;
    double cellsize;
    double (*data)[2];

    double xMin;
    double xMax;
    double yMin;
    double yMax;

    double radius;
    bool firstpoint;
    double (*pointQueue)[2];
    int pointQueueCapacity;
    int pointqueuelen;
    int k;
};

double parse(double x, double k);
void deletePoint(double (*pointQueue)[2], int *count, int idx);
int initialise_param(struct GRID *g, double *viewport, double minDist, int maxTries);
double dist2(double x1, double y1, double x2, double y2);
int createNewPoint(struct GRID *g, double x, double y, double dict[2]);
bool isValidPoint(struct GRID *g, double x, double y);
int nextPoint(struct GRID *g, pcg64_t *rng, double dict[2]);
bool checkPoints(struct GRID *g, double (*points)[2], int lenpoints, double minDist);
int allPoints(struct GRID *g, double (**result)[2], pcg64_t *rng);
int PoissonDisk(struct GRID *g, double *viewport, double minDist, int maxTries, double (**result)[2], pcg64_t *rng);
void freeGrid(struct GRID *g);

#endif