#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "poisson_disk.h"
#include "pcg64.h"

const double SQRT_HALF = 0.70710678118654752440;
const double PI = 3.14159265358979323846;

double parse(double x, double k)
{
    if (isnan(x))
    {
        return k;
    }
    else
    {
        return x;
    }
}

// count = current number of points
void deletePoint(double (*pointQueue)[2], int *count, int idx)
{
    for (int i = idx; i < *count - 1; i++)
    {
        pointQueue[i][0] = pointQueue[i + 1][0];
        pointQueue[i][1] = pointQueue[i + 1][1];
    }

    (*count)--;
}

int initialise_param(struct GRID *g, double *viewport, double minDist, int maxTries)
{
    g->xMin = parse(viewport[0], 0);
    g->xMax = parse(viewport[1], 0);
    g->yMin = parse(viewport[2], 0);
    g->yMax = parse(viewport[3], 0);

    g->radius = fmax(parse(minDist, 1), 1);
    g->k = (int)(fmax(parse(maxTries, 30), 2));
    g->firstpoint = true;
    g->cellsize = (g->radius) * (SQRT_HALF);
    g->width = (ceil((g->xMax - g->xMin) / g->cellsize));
    g->height = (ceil((g->yMax - g->yMin) / g->cellsize));

    g->data = malloc(g->width * g->height * sizeof(*g->data));
    if (g->data == NULL)
    {
        freeGrid(g);
        return 0;
    }

    for (int i = 0; i < g->width * g->height; i++)
    {
        g->data[i][0] = NAN;
        g->data[i][1] = NAN;
    }
    g->pointQueueCapacity = 128;
    g->pointQueue = malloc(g->pointQueueCapacity * sizeof(*g->pointQueue));
    if (g->pointQueue == NULL)
    {
        freeGrid(g);
        return 0;
    }

    g->pointqueuelen = 0;

    return 1;
}

double dist2(double x1, double y1, double x2, double y2)
{
    return ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1));
}

int createNewPoint(struct GRID *g, double x, double y, double dict[2])
{
    dict[0] = x;
    dict[1] = y;

    double col = fmin(floor((x - g->xMin) / g->cellsize), g->width - 1);
    double row = fmin(floor((y - g->yMin) / g->cellsize), g->height - 1);

    int index = col + row * (g->width);
    g->data[index][0] = x;
    g->data[index][1] = y;

    if (g->pointqueuelen >= g->pointQueueCapacity)
    {
        double (*temp)[2] = realloc(g->pointQueue, 2 * g->pointQueueCapacity * sizeof(*g->pointQueue));
        if (temp == NULL)
        {
            return 0;
        }
        g->pointQueue = temp;
        g->pointQueueCapacity *= 2;
    }

    g->pointQueue[g->pointqueuelen][0] = x;
    g->pointQueue[g->pointqueuelen][1] = y;
    g->pointqueuelen += 1;

    return 1;
}

bool isValidPoint(struct GRID *g, double x, double y)
{
    if (x < g->xMin || x > g->xMax || y < g->yMin || y > g->yMax)
    {
        return false;
    }

    if (!isfinite(x) || !isfinite(y))
    {
        return false;
    }

    double col = floor((x - g->xMin) / g->cellsize);
    double row = floor((y - g->yMin) / g->cellsize);

    for (int i = col - 2; i <= col + 2; i++)
    {
        for (int j = row - 2; j <= row + 2; j++)
        {
            if (i >= 0 && i < g->width && j >= 0 && j < g->height)
            {
                int idx = i + (j * g->width);
                if (!(isnan(g->data[idx][0]) && isnan(g->data[idx][1])) && dist2(x, y, (g->data[idx][0]), (g->data[idx][1])) <= g->radius * g->radius)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

int nextPoint(struct GRID *g, pcg64_t *rng, double dict[2])
{
    double x = 0;
    double y = 0;

    if (g->firstpoint)
    {
        g->firstpoint = false;
        x = g->xMin + (g->xMax - g->xMin) * pcg64_double(rng);
        y = g->yMin + (g->yMax - g->yMin) * pcg64_double(rng);
        int go = createNewPoint(g, x, y, dict);
        if (go == 0)
        {
            return -1;
        }

        return 1;
    }

    while (g->pointqueuelen > 0)
    {
        int idx = (int)(g->pointqueuelen * pcg64_double(rng));
        for (int i = 0; i < g->k; i++)
        {
            double distance = g->radius * (pcg64_double(rng) + 1);
            double angle = 2 * PI * pcg64_double(rng);
            double x = g->pointQueue[idx][0] + distance * cos(angle);
            double y = g->pointQueue[idx][1] + distance * sin(angle);
            if (isValidPoint(g, x, y))
            {
                int go = createNewPoint(g, x, y, dict);
                if (go == 0)
                {
                    return -1;
                }
                return 1;
            }
        }
        deletePoint(g->pointQueue, &g->pointqueuelen, idx);
    }
    return 0; // None
}

bool checkPoints(struct GRID *g, double (*points)[2], int lenpoints, double minDist)
{
    double minDist2 = minDist * minDist;

    for (int i = 0; i < lenpoints; i++)
    {
        for (int j = i + 1; j < lenpoints; j++)
        {
            double dx = points[i][0] - points[j][0];
            double dy = points[i][1] - points[j][1];

            if ((dx * dx) + (dy * dy) < minDist2)
            {
                return false;
            }
        }
    }

    return true;
}

int allPoints(struct GRID *g, double (**result)[2], pcg64_t *rng)
{
    int capacity = 128;
    int lenres = 0;

    *result = malloc(capacity * sizeof(**result));
    if (*result == NULL)
    {
        free(*result);
        *result = NULL;
        freeGrid(g);
        return -1;
    }

    double point[2];
    int s;

    while (s = nextPoint(g, rng, point))
    {
        if (s == -1)
        {
            free(*result);
            *result = NULL;
            freeGrid(g);
            return -1;
        }

        if (lenres >= capacity)
        {
            double (*tmp)[2] = realloc(*result, 2 * capacity * sizeof(**result));
            if (tmp == NULL)
            {
                free(*result);
                freeGrid(g);
                *result = NULL;
                return -1;
            }

            *result = tmp;
            capacity *= 2;
        }

        (*result)[lenres][0] = point[0];
        (*result)[lenres][1] = point[1];

        lenres++;
    }

    return lenres;
}

int PoissonDisk(struct GRID *g, double *viewport, double minDist, int maxTries, double (**result)[2], pcg64_t *rng)
{
    if (g == NULL || result == NULL || rng == NULL || viewport == NULL)
    {
        return 0;
    }

    if (viewport[0] >= viewport[1] || viewport[2] >= viewport[3] ||
        !isfinite(viewport[0]) || !isfinite(viewport[1]) ||
        !isfinite(viewport[2]) || !isfinite(viewport[3]) ||
        !isfinite(minDist))
    {
        return 0;
    }

    if (minDist <= 0 || maxTries < 0)
    {
        return 0;
    }

    g->data = NULL;
    g->pointQueue = NULL;
    g->pointqueuelen = 0;
    g->pointQueueCapacity = 0;
    g->width = 0;
    g->height = 0;
    g->cellsize = 0;
    g->xMin = 0;
    g->xMax = 0;
    g->yMin = 0;
    g->yMax = 0;
    g->radius = 0;
    g->firstpoint = NULL;
    g->k = 0;

    *result = NULL;

    int start = initialise_param(g, viewport, minDist, maxTries);
    if (start == 0)
    {
        return 0;
    }

    int all = allPoints(g, result, rng);
    if (all == -1)
    {
        return 0;
    }

    return all;
}

void freeGrid(struct GRID *g)
{
    if (!g)
    {
        return;
    }
    free(g->data);
    g->data = NULL;
    free(g->pointQueue);
    g->pointQueue = NULL;
    g->pointqueuelen = 0;
    g->pointQueueCapacity = 0;
    g->width = 0;
    g->height = 0;
}