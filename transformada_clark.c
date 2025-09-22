#include "transformada_clark.h"

#define DOIS_DIV_3    0.6666667f     // 2/3
#define HALF          0.5f
#define SQRT3_BY_2    0.8660254f

void ClarkeTransform(volatile  float a,volatile  float b,volatile  float c, volatile  float *alpha, volatile  float *beta)
{
    *alpha = DOIS_DIV_3 * (a - HALF * b - HALF * c);
    *beta  = DOIS_DIV_3 * (SQRT3_BY_2 * b - SQRT3_BY_2 * c);
}
