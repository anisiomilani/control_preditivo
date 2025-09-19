#ifndef FUNCOES_H_
#define FUNCOES_H_

void atualizar_valor_par_impar_vetor_6(volatile float *v, int inicio, volatile float novo_valor);
void atualizar_valor_par_impar_vetor_4(volatile float *v, int inicio, volatile float novo_valor);


void extrapolar_k1( float V[6],  float Vk1[2]);
void extrapolar_k2( float V[6],  float Vk2[2]);
float calcular_q_filtrado( float *entrada_filtrada,  float *saida_q_filtrada, int entrada_atual, int entrada_ant, int saida_ant);
float calcular_v_filtrado( float *entrada,  float *saida_filtrada,int entrada_atual, int entrada_ant,int saida_ant_1, int saida_ant_2);
#endif/* FUNCOES_H_ */
