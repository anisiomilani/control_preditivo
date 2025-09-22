#include <stdint.h>

void atualizar_valor_par_impar_vetor_6(volatile  float *v,   uint16_t inicio,  volatile  float novo_valor) {
    // inicio == 0  pares: 0, 2, 4
    // inicio == 1  impares: 1, 3, 5
    v[inicio + 4] = v[inicio + 2];
    v[inicio + 2] = v[inicio + 0];
    v[inicio + 0] = novo_valor;
}

void atualizar_valor_par_impar_vetor_4( volatile float *v, uint16_t inicio,  volatile float novo_valor) {
    v[inicio + 2] = v[inicio];     // posição 0 ou 1 → vai para 2 ou 3
    v[inicio] = novo_valor;        // atualiza valor mais recente
}

void extrapolar_k1(volatile float V[6], volatile float Vk1[2])
{
    Vk1[0] = 3.0f * V[0] - 3.0f * V[2] + V[4];
    Vk1[1] = 3.0f * V[1] - 3.0f * V[3] + V[5];
}

void extrapolar_k2(volatile float V[6], volatile float Vk2[2])
{
    Vk2[0] = 6.0f * V[0] - 8.0f * V[2] + 3.0f * V[4];
    Vk2[1] = 6.0f * V[1] - 8.0f * V[3] + 3.0f * V[5];
}

