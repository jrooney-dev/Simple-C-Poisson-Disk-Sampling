# Poisson-Disk Sampling in C

A compact C implementation of **2D Poisson-disk sampling**. It generates points inside a rectangular viewport while maintaining a minimum separation between accepted points.

The sampler uses a spatial grid and an active-point queue to check candidate positions efficiently. It is intended for ordinary use with valid, finite inputs and a fresh grid for each generation.

This implementation is inspired by [ogus/poisson-disk](https://github.com/ogus/poisson-disk), which is based on Robert Bridson's paper, [*Fast Poisson Disk Sampling in Arbitrary Dimensions*](https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf). The approach is designed to run in **O(n)** time.

## Features:

- Generates 2D Poisson-disk point distributions within a rectangular viewport.
- Uses a grid-based neighbourhood search and an active-point queue.
- Accepts a caller-provided PCG64 random-number-generator state.
- Returns the generated points as a dynamically allocated array.

## Requirements:

- A C compiler with C11 support.
- The project's `poisson_disk.c` and `poisson_disk.h` files.
- The PCG64 implementation used by the project (`pcg64.h` and its implementation/library).
- The math library (`libm`) on systems where it is linked separately.

## Quick Start:

```c
#include <stdio.h>
#include <stdlib.h>

#include "poisson_disk.h"
#include "pcg64.h"

int main(void)
{
    /* Use a fresh, zero-initialized grid for each generation. */
    struct GRID grid = {0};

    /* Bounds are ordered as {xMin, xMax, yMin, yMax}. */
    double viewport[4] = {0.0, 100.0, 0.0, 100.0};
    double (*points)[2] = NULL;

    pcg64_t rng;
    pcg64_seed(&rng, 42, 54);

    int count = PoissonDisk(
        &grid,
        viewport,
        2.0,       /* minimum point distance */
        30,        /* candidate attempts per active point */
        &points,
        &rng
    );

    if (count <= 0) {
        fprintf(stderr, "Point generation failed.\\n");
        free(points);
        freeGrid(&grid);
        return 1;
    }

    for (int i = 0; i < count; ++i) {
        printf("%f %f\\n", points[i][0], points[i][1]);
    }

    free(points);
    freeGrid(&grid);
    return 0;
}
```

Link the sampler, the PCG64 implementation, and the math library. For example, if the PCG64 implementation is available as `libpcg64.a`:

```sh
cc -std=c11 -Wall -Wextra -Wpedantic -Iinclude main.c poisson_disk.c -Llib -lpcg64 -lm -o bin/poisson_demo
```

Adjust source filenames and linker flags to match your project.

## API:

```c
int PoissonDisk(
    struct GRID *g,
    double *viewport,
    double minDist,
    int maxTries,
    double (**result)[2],
    pcg64_t *rng
);
```

### Parameters:

| Parameter  | Description                                                                                                  |
| ---------- | ------------------------------------------------------------------------------------------------------------ |
| `g`        | Grid state used internally by the sampler. Pass a fresh, zero-initialized `struct GRID`.                     |
| `viewport` | Viewport bounds in the order `{xMin, xMax, yMin, yMax}`.                                                     |
| `minDist`  | Requested minimum separation between points. The implementation clamps the effective radius to at least `1`. |
| `maxTries` | Candidate attempts per active point. The implementation uses at least `2` attempts.                          |
| `result`   | Receives a dynamically allocated array of `{x, y}` pairs. Initialize the pointer to `NULL` before calling.   |
| `rng`      | Initialized PCG64 random-number-generator state.                                                             |

### Return Value:

The function returns the number of generated points. In the current implementation, **a return value of `0` indicates failure**, so callers should treat non-positive results as unsuccessful.

## Viewport Format:

The C implementation expects the bounds in this exact order:

```c
double viewport[4] = {xMin, xMax, yMin, yMax};
```

This ordering may differ from the interface used by the JavaScript project that inspired this implementation.

## Memory Management and Reuse:

- Initialize the grid before each generation:

  ```c
  struct GRID grid = {0};
  ```

- Initialize the result pointer before calling the sampler:

  ```c
  double (*points)[2] = NULL;
  ```

- Release the returned point array with `free(points)`.

- Release the grid's internal allocations with `freeGrid(&grid)`.

- Use a **fresh grid for every generation**. Do not call `PoissonDisk()` again on a grid that still owns allocations: the current setup resets grid fields and is intended for a fresh grid.

Example cleanup:

```c
free(points);
freeGrid(&grid);
```

## Validating Generated Points:

The project provides `checkPoints()` to check pairwise minimum separation:

```c
if (!checkPoints(&grid, points, count, 2.0)) {
    fprintf(stderr, "A point pair violates the minimum-distance check.\\n");
}
```

This helper checks pairwise distances only. It does not verify viewport membership, input validity, or that all generated coordinates are finite.

## Input Expectations:

Use:

- finite viewport bounds;
- `xMin < xMax` and `yMin < yMax`;
- a positive `minDist`;
- a non-negative `maxTries`; and
- viewport dimensions and point counts within reasonable limits.

The API returns a point count rather than a separate status code. As a result, a return value of `0` does not provide detailed failure information.

## Current Limitations:

This is a practical C implementation, not a fully hardened, production-ready general-purpose library. The following limitations are known.

### API Differences From the JavaScript Project:

The implementation follows the same general algorithmic approach, but it is not a drop-in port of the original `PoissonDisk` object API. Differences include:

- procedural API rather than an object-style API;
- different viewport ordering;
- different lifecycle and reuse semantics; and
- different result handling and memory ownership.

### Minimal Error Reporting:

The function returns a point count and does not provide a detailed error status. A return value of `0` can make it difficult to distinguish between invalid input, allocation failure, an impossible sampling case, or another runtime failure.

### Grid Lifecycle Requirements:

The sampler assumes a fresh, zero-initialized `GRID` for each generation. Reusing a grid without the appropriate cleanup or reset may result in stale pointers, memory leaks, invalid frees, or undefined behaviour.

### Large Dimensions and Allocation Overflow:

Very large viewports may expose unchecked overflow in dimension calculations or allocation-size multiplication, as well as invalid allocation sizes or out-of-range grid indexing. Keep inputs within reasonable limits.

### Limited Validation Helper:

`checkPoints()` is a useful pairwise-spacing check, but it does not validate:

- pointer correctness;
- viewport bounds or viewport membership;
- NaN or infinite coordinates; or
- whether every generated value is finite.

It is a testing aid, not a complete correctness guarantee.

## Suitable Use Cases:

For ordinary, well-formed inputs, the implementation is suitable for:

- small and medium-sized viewports;
- one-off test runs;
- local scientific or graphics prototypes; and
- simple C-based sampling workflows.

For production use with untrusted or extreme inputs, additional input validation, overflow checks, clearer error reporting, and safer grid lifecycle handling are recommended.
