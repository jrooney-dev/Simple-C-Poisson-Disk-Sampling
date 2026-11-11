#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "poisson_disk.h"
#include "pcg64.h"

int main()
{
    struct GRID grid;

    double viewport[4] = {0, 100, 0, 100};

    double minDist = 2.0;
    int maxTries = 30;

    double (*points)[2];

    pcg64_t rng;
    pcg64_seed(&rng, 42, 54);

    printf("Starting\n");
    fflush(stdout);

    int count = PoissonDisk(
        &grid,
        viewport,
        minDist,
        maxTries,
        &points,
        &rng);

    printf("Finished\n");

    printf("Generated %d points\n", count);

    for (int i = 0; i < count; i++)
    {
        printf("%f %f\n", points[i][0], points[i][1]);
    }

    bool a = checkPoints(&grid, points, count, minDist);
    if (a == true)
    {
        printf("%s ", "True");
    }
    else if (a == false)
    {
        printf("%s ", "False");
    }
    else
    {
        printf("%s ", "Error");
    }

    FILE *f = fopen("C://Users//admin/Desktop//points.csv", "w");

    for (int i = 0; i < count; i++)
    {
        fprintf(f, "%f,%f\n",
                points[i][0],
                points[i][1]);
    }

    fclose(f);

    free(points);
    freeGrid(&grid);

    return 0;
}