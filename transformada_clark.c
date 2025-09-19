#include "transformada_clark.h"

#define DOIS_DIV_3    0.6666667f     // 2/3
#define HALF          0.5f
#define SQRT3_BY_2    0.8660254f

void ClarkeTransform(float a, float b, float c,  float *alpha,float *beta)
{
    *alpha = DOIS_DIV_3 * (a - HALF * b - HALF * c);
    *beta  = DOIS_DIV_3 * (SQRT3_BY_2 * b - SQRT3_BY_2 * c);
}
