#include "driverlib.h"
#include "device.h"
#include "board.h"
#include "math.h"
#include "transformada_clark.h"
#include "funcoes.h"


//ainda não esta sendo utilizado
#define norm_DAC 4095.0f/(84.0f)
#define norm_DAC_il 4095.0f/(8.4f)
#define LIMIAR_REARME_ADC 40.0f


// varaveis criadas para  PWM
uint32_t ePwm_TimeBase;
uint32_t ePwm_MinDuty;
uint32_t ePwm_MaxDuty;
uint32_t ePwm_curDuty;

volatile bool g_trip_clear = false;

// Definições de Constantes
//
#define F_PWM                  10000.0f     // Frequência de chaveamento (Hz)
#define T_PWM                  (1.0f / F_PWM) // Período de chaveamento (s)
#define DT_SIM                 0.000001f    // Passo de simulação (5 µs)
#define N_STEPS_PER_CYCLE      (uint32_t)(T_PWM / DT_SIM) // Passos por ciclo PWM

volatile uint32_t g_step_counter = 0;  // Contador de passos dentro do ciclo PWM
volatile bool g_switch_on = false;           // Estado da chave (true = ligada)
volatile bool g_new_step_ready = false;     // Flag para novo passo de simulação
//volatile float g_duty_cycle = 0.5f;          // Razão cíclica (entre 0 e 1)



//------------------------------------------INICIO DO CONTROLE-------------------------------------------------------------------------


//                  Constantes pré-calculadas para transformada de clark
//-------------------------------------------------------------------------------------------------------
#define TRES_DIV_2 1.5f                                  //    3/2
#define DOIS_DIV_3 0.6666667f                           //     2/3
#define HALF 0.5f
#define SQRT3_BY_2 0.8660254f                          //      sqrt(3)/2
#define UM_DIV_3  0.333f                              //       1.0/3
#define SQRT3_BY_3 0.57735f                          //        sqrt(3)/3
#define INV_RV 0.147f                               //         (1/RV)   RV=6.8184;
#define INV_2_RV 0.0733f                           //          (1/2*RV)
#define LIMITAR_ZERO 0.00000000001f
#define LIMITADOR_UM 1.0f


// valores definidos para o filtro SOGI com k=1
#define INV_TS  40000.0f    // 1/Ts Ts= 0.000025;
#define Lc 0.00584f
#define Lg 0.00106f
#define rc 0.2f
#define rg 0.17f
#define C 0.0000114f
#define ZERO 0.0f
#define UM 1.0f
#define INICIO_OPERACIONAL 1e15f

// Vetor chaveamento alfa/beta

static const float s_ab[8][2] =     { {ZERO,          ZERO },
                                    {DOIS_DIV_3,    ZERO },
                                    {UM_DIV_3,      SQRT3_BY_3 },
                                    {-UM_DIV_3,     SQRT3_BY_3 },
                                    {-DOIS_DIV_3,   ZERO },
                                    {-UM_DIV_3,     -SQRT3_BY_3 },
                                    {UM_DIV_3,      -SQRT3_BY_3 },
                                    {ZERO,          ZERO } };

// Vetor chaveamento trifásico
static const int16_t s_abc[8][3] = { {ZERO,     ZERO,       ZERO },
                                     {UM,       ZERO,       ZERO },
                                     {UM,       UM,         ZERO },
                                     {ZERO,     UM,         ZERO },
                                     {ZERO,     UM,         UM },
                                     {ZERO,     ZERO,       UM },
                                     {UM,       ZERO,       UM },
                                     {UM,       UM,         UM } };

// custo para utilizar no controle preditivo
#define g_custo_i1 1.0f
#define g_custo_i2 1.0f
#define g_custo_vcap 1.0f
#define g_custo_s 1.0f

// constantes utilizadas para predizer amotras a frente
#define Const 0.602655f          // Const=((1/RV)+(C/Ts))
#define const_5 0.456f          // const_5=(C/Ts)
#define phi_1 0.999145f         // phi_1 = 1 - ((rc*Ts)/Lc)
#define gama_1 0.00428082f      //  gama_1 = Ts/Lc
#define phi_2 0.995991f         // phi_2 = 1 - ((rg*Ts)/Lg)
#define gama_2 0.0235849f       // gama_2 = Ts/Lg
#define phi_3 0.6785f           // phi_3 = 1 - (Ts/(C*RV))
#define gama_3 2.19298f         // gama_3 = Ts/C


// valores definidos para o filtro SOGI com k=1
#define a11 0.00471238898f
#define a22 0.00471238898f
#define b1  1.990531226413290f
#define b2  0.990619634277416f
#define a1  0.004690182861292f
#define a2  0.004690182861292f

//-------------------------------------------------------------------------------------------------------------
//variaveis que seram lidas trocar por leitura adc depois
// Tensões e correntes de PAC e CAP
volatile float va_pac, vb_pac, vc_pac;
volatile float i1_a, i1_b, i1_c;
volatile float va_cap, vb_cap, vc_cap;
volatile float i2_a, i2_b, i2_c;
volatile float Vdc =700.0f;



int16_t linha_op[2] = { 0, 0 };
int16_t passado = 0;
volatile int16_t K = 0;
volatile float P_ref = 5000.0f;
volatile float Q_ref = 0.0f;

// Vetores float inicializados posição por posição

volatile float Vg_ab[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i1_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i2_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float vcap_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float Vg_ab_k1[2] = { 0.0f, 0.0f };
volatile float Vg_ab_k2[2] = { 0.0f, 0.0f };
volatile float Vg_ab_filtrado[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float Vg_ab_q_filtrado[4] = { -3.0f, 5.0f, -2.0f, 1.0f };
volatile float Vg_ab_filtrado_neg[2] = { 0.0f, 0.0f };
volatile float Vg_ab_filtrado_pos[2] = { 0.0f, 0.0f };
volatile float i2_ref_ab[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i2_ref_ab_k2[2] = { 0.0f, 0.0f };
volatile float i2_ref_ab_virtual[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float vcap_ref_ab[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float vcap_ref_ab_k2[2] = { 0.0f, 0.0f };
volatile float derivada_i2[2] = { 0.0f, 0.0f };
volatile float i1_ref_ab[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i1_ref_ab_k2[2] = { 0.0f, 0.0f };
volatile float i1_ab_k1[2] = { 0.0f, 0.0f };
volatile float i2_ab_k1[2] = { 0.0f, 0.0f };
volatile float vcap_ab_k1[2] = { 0.0f, 0.0f };
volatile float i1_ab_k2[2] = { 0.0f, 0.0f };
volatile float i2_ab_k2[2] = { 0.0f, 0.0f };
volatile float vcap_ab_k2[2] = { 0.0f, 0.0f };


void main(void)
{

    // Inicialização dos periféricos

    Device_init();
    Interrupt_initModule();
    Interrupt_initVectorTable();
    Board_init();

    ePwm_TimeBase = EPWM_getTimeBasePeriod(EPWM0_BASE);
    ePwm_MinDuty = (uint32_t) (0.95f * (float) ePwm_TimeBase);
    ePwm_MaxDuty = (uint32_t) (0.05f * (float) ePwm_TimeBase);

    EINT;
    ERTM;

    while (1)
    {

        // Componentes alfa-beta
        volatile float valfa_pac, vbeta_pac;
        volatile float i1_alfa, i1_beta;
        volatile float valfa_cap, vbeta_cap;
        volatile float i2_alfa, i2_beta;
        volatile float i2_ref_alfa, i2_ref_alfa_pos, i2_ref_alfa_neg;
        volatile float i2_ref_beta, i2_ref_beta_pos, i2_ref_beta_neg;

        // Ganhos, potências, controle
        volatile float g_op, D, E, inv_D, inv_E;
        volatile float P_ativa, Q_reativa;
        volatile float gcap_alfa, gcap_beta;

        volatile float g2_alfa, g2_beta, g1_alfa, g1_beta;
        volatile int16_t gsa, gsb, gsc;
        volatile float  g_k;
        volatile float Vg_ab_filtrado_pos_0_quadrado, Vg_ab_filtrado_pos_1_quadrado;
        volatile float Vg_ab_filtrado_neg_0_quadrado, Vg_ab_filtrado_neg_1_quadrado;
        volatile int16_t linha;
        volatile float Vk_ab[2];


        //  cmp_Value = (uint32_t) (g_duty_cycle * ePwm_TimeBase);
        //  EPWM_setCounterCompareValue(EPWM0_BASE, EPWM_COUNTER_COMPARE_A, cmp_Value);
        //   ePwm_curDuty = EPWM_getCounterCompareValue(EPWM0_BASE, EPWM_COUNTER_COMPARE_A);

        //     if (g_new_step_ready)
        //      {
        //         g_new_step_ready = false;
//----------------------------------------------------
        passado = linha_op[0];
        linha_op[1] = linha_op[0];
        linha_op[0] = 1000;
        g_op = INICIO_OPERACIONAL;
// ========================================== 1. Transformadas CLARK ==========================================================================

        ClarkeTransform(va_pac, vb_pac, vc_pac, &valfa_pac, &vbeta_pac);
        ClarkeTransform(i1_a, i1_b, i1_c, &i1_alfa, &i1_beta);
        ClarkeTransform(va_cap, vb_cap, vc_cap, &valfa_cap, &vbeta_cap);
        ClarkeTransform(i2_a, i2_b, i2_c, &i2_alfa, &i2_beta);


// ========================================== 2. Redefinindo as variaveis para vetores =========================================================
        // operador 0 é para PAR
        // operador 1 é para IMPAR

        atualizar_valor_par_impar_vetor_6(Vg_ab, 0, valfa_pac);
        atualizar_valor_par_impar_vetor_6(Vg_ab, 1, vbeta_pac);

        atualizar_valor_par_impar_vetor_4(vcap_ab, 0, valfa_cap);
        atualizar_valor_par_impar_vetor_4(vcap_ab, 1, vbeta_cap);

        atualizar_valor_par_impar_vetor_4(i1_ab, 0, i1_alfa);
        atualizar_valor_par_impar_vetor_4(i1_ab, 1, i1_beta);

        atualizar_valor_par_impar_vetor_4(i2_ab, 0, i2_alfa);
        atualizar_valor_par_impar_vetor_4(i2_ab, 1, i2_beta);

// ========================================== 3. Extrapolações ==========================================================================

        extrapolar_k1(Vg_ab, Vg_ab_k1);
        extrapolar_k2(Vg_ab, Vg_ab_k2);

// =========================================== 4. Filtro SOGI ============================================================================

        atualizar_valor_par_impar_vetor_6(Vg_ab_filtrado, 0, Vg_ab_filtrado[0]);
        atualizar_valor_par_impar_vetor_6(Vg_ab_filtrado, 1, Vg_ab_filtrado[1]);
        Vg_ab_filtrado[0]= a1 * Vg_ab[0] - a2 * Vg_ab[4] + b1 * Vg_ab_filtrado[2]- b2 * Vg_ab_filtrado[4];
        Vg_ab_filtrado[1]= a1 * Vg_ab[1] - a2 * Vg_ab[5] + b1 * Vg_ab_filtrado[3]- b2 * Vg_ab_filtrado[5];

        //   componente em quadratura

        atualizar_valor_par_impar_vetor_4(Vg_ab_q_filtrado, 0, Vg_ab_q_filtrado[0]);
        atualizar_valor_par_impar_vetor_4(Vg_ab_q_filtrado, 1, Vg_ab_q_filtrado[1]);
        Vg_ab_q_filtrado[0]= a11 * Vg_ab_filtrado[0] + a22 * Vg_ab_filtrado[2] + Vg_ab_q_filtrado[2];
        Vg_ab_q_filtrado[1]= a11 * Vg_ab_filtrado[1] + a22 * Vg_ab_filtrado[3] + Vg_ab_q_filtrado[3];

// =========================================== 6. Sequência positiva e negativa =============================================================

        Vg_ab_filtrado_pos[0] = HALF * (Vg_ab_filtrado[0] - Vg_ab_q_filtrado[1]);
        Vg_ab_filtrado_pos[1] = HALF * (Vg_ab_filtrado[1] + Vg_ab_q_filtrado[0]);

        Vg_ab_filtrado_neg[0] = HALF * (Vg_ab_filtrado[0] + Vg_ab_q_filtrado[1]);
        Vg_ab_filtrado_neg[1] = HALF * (Vg_ab_filtrado[1] - Vg_ab_q_filtrado[0]);

        if (Vg_ab_filtrado_pos[0] == ZERO) Vg_ab_filtrado_pos[0] = LIMITAR_ZERO;
        if (Vg_ab_filtrado_pos[1] == ZERO) Vg_ab_filtrado_pos[1] = LIMITAR_ZERO;

// =========================================== 7. Utilizando as sequencia para determinar a referencia conforme i2 ===============================
// =========================================== o valor definido para K =======================================================================
// =========================================== k = 0 sistema senoidal na rede==============================================================
// =========================================== k = 1 controle potencia ativa==============================================================
// =========================================== k = -1 controle potencia reativa==============================================================

        Vg_ab_filtrado_pos_0_quadrado = Vg_ab_filtrado_pos[0]* Vg_ab_filtrado_pos[0];
        Vg_ab_filtrado_pos_1_quadrado = Vg_ab_filtrado_pos[1]* Vg_ab_filtrado_pos[1];
        Vg_ab_filtrado_neg_0_quadrado = Vg_ab_filtrado_neg[0]* Vg_ab_filtrado_neg[0];
        Vg_ab_filtrado_neg_1_quadrado = Vg_ab_filtrado_neg[1] * Vg_ab_filtrado_neg[1];

        D = (Vg_ab_filtrado_pos_0_quadrado + Vg_ab_filtrado_pos_1_quadrado)- K* (Vg_ab_filtrado_neg_0_quadrado+ Vg_ab_filtrado_neg_1_quadrado);
        E = (Vg_ab_filtrado_pos_0_quadrado + Vg_ab_filtrado_pos_1_quadrado)+ K* (Vg_ab_filtrado_neg_0_quadrado+ Vg_ab_filtrado_neg_1_quadrado);

        // Proteger contra divisão por zero
        if (D == ZERO) D = LIMITAR_ZERO;
        if (E == ZERO) E = LIMITAR_ZERO;

        // Fazer só uma divisão para cada
        inv_D = 1.0f / D;
        inv_E = 1.0f / E;

        i2_ref_alfa_pos = DOIS_DIV_3 * (((Vg_ab_filtrado_pos[0] * P_ref) * inv_D) + ((Vg_ab_filtrado_pos[1] * Q_ref) * inv_E));
        i2_ref_beta_pos = DOIS_DIV_3 * (((Vg_ab_filtrado_pos[1] * P_ref) * inv_D) - ((Vg_ab_filtrado_pos[0] * Q_ref) * inv_E));

        i2_ref_alfa_neg = DOIS_DIV_3 * ((-K * ((Vg_ab_filtrado_neg[0] * P_ref) * inv_D)) + (K * ((Vg_ab_filtrado_neg[1] * Q_ref) * inv_E)));
        i2_ref_beta_neg = DOIS_DIV_3 * ((-K * ((Vg_ab_filtrado_neg[1] * P_ref) * inv_D)) - (K * ((Vg_ab_filtrado_neg[0] * Q_ref) * inv_E)));

        atualizar_valor_par_impar_vetor_6(i2_ref_ab, 0, i2_ref_ab[0]);
        atualizar_valor_par_impar_vetor_6(i2_ref_ab, 1, i2_ref_ab[1]);
        i2_ref_ab[0] = i2_ref_alfa_pos + i2_ref_alfa_neg;
        i2_ref_ab[1] = i2_ref_beta_pos + i2_ref_beta_neg;

         // Extrapolação da corrente de referencia i2

         extrapolar_k2(i2_ref_ab, i2_ref_ab_k2);


// =========================================== 8. Calculo corrente virtual ===================================================================


        //foi realizado uma media entre a tensao medida no capacitor e tensao referencia capacitor
        atualizar_valor_par_impar_vetor_6(i2_ref_ab_virtual, 0, i2_ref_ab_virtual[0]);
        atualizar_valor_par_impar_vetor_6(i2_ref_ab_virtual, 1, i2_ref_ab_virtual[1]);
        i2_ref_ab_virtual[0]=i2_ref_ab[0] - (1.55f * vcap_ref_ab[0] + 0.45f * vcap_ab[0]) * INV_2_RV;
        i2_ref_ab_virtual[1]=i2_ref_ab[1] - (1.55f * vcap_ref_ab[1] + 0.45f * vcap_ab[1]) * INV_2_RV;
// =========================================== 9. Calculo potencia ativa e reativa real  ===================================================================

        P_ativa = TRES_DIV_2 * (Vg_ab[0] * i2_ab[0] + Vg_ab[1] * i2_ab[1]);

        Q_reativa = TRES_DIV_2 * (Vg_ab[1] * i2_ab[0] - Vg_ab[0] * i2_ab[1]);

// =========================================== 10. Gerando Ref - Vcap  ===================================================================

        derivada_i2[0] = (i2_ref_ab_virtual[0] - i2_ref_ab_virtual[2]) * INV_TS;
        derivada_i2[1] = (i2_ref_ab_virtual[1] - i2_ref_ab_virtual[3]) * INV_TS;

        //  limitador=1;    //Limitando a derivada para remover picos na transição

        derivada_i2[0] = (derivada_i2[0] > LIMITADOR_UM) ? LIMITADOR_UM : derivada_i2[0];
        derivada_i2[0] = (derivada_i2[0] < -LIMITADOR_UM) ? -LIMITADOR_UM : derivada_i2[0];
        derivada_i2[1] = (derivada_i2[1] > LIMITADOR_UM) ? LIMITADOR_UM : derivada_i2[1];
        derivada_i2[1] = (derivada_i2[1] < -LIMITADOR_UM) ? -LIMITADOR_UM : derivada_i2[1];

        atualizar_valor_par_impar_vetor_6( vcap_ref_ab, 0, vcap_ref_ab[0]);
        atualizar_valor_par_impar_vetor_6( vcap_ref_ab, 1, vcap_ref_ab[1]);
        vcap_ref_ab[0]=Vg_ab[0] + Lg * derivada_i2[0] + rg * i2_ref_ab_virtual[0];
        vcap_ref_ab[1]=Vg_ab[1] + Lg * derivada_i2[1] + rg * i2_ref_ab_virtual[1];

        //extrapolação da tensão no capacitor
        extrapolar_k2(vcap_ref_ab, vcap_ref_ab_k2);

// =========================================== 11. Gerando Ref da corrente i1 com resistor em paralelo com o capacitor  ===================================================================
//----------------------------------------       Amortecimento do filtro LCL---------------------------------------------

        atualizar_valor_par_impar_vetor_6(i1_ref_ab, 0, i1_ref_ab[0]);
        atualizar_valor_par_impar_vetor_6(i1_ref_ab, 1,  i1_ref_ab[1]);
        i1_ref_ab[0]=const_5 * (vcap_ref_ab[0] - vcap_ref_ab[2]) + (vcap_ref_ab[0] * INV_RV) + i2_ref_ab_virtual[0];
        i1_ref_ab[1]=const_5 * (vcap_ref_ab[1] - vcap_ref_ab[3]) + (vcap_ref_ab[1] * INV_RV) + i2_ref_ab_virtual[1];

        //extrapolação da corrente i1
        extrapolar_k2(i1_ref_ab, i1_ref_ab_k2);

// =========================================== 12. Entrada do algoritmo de controle  ===================================================================

//---------------------------                 Estimando os valores em (K+1) do preditivo----------------------------------------------------------

         Vk_ab[0] = Vdc * s_ab[passado][0];
         Vk_ab[1] = Vdc * s_ab[passado][1];

         i1_ab_k1[0] = (phi_1 * i1_ab[0]) + gama_1 * (Vk_ab[0] - vcap_ab[0]);
         i1_ab_k1[1] = (phi_1 * i1_ab[1]) + gama_1 * (Vk_ab[1] - vcap_ab[1]);

         i2_ab_k1[0] = (phi_2 * i2_ab[0]) + gama_2 * (vcap_ab[0] - Vg_ab[0]);
         i2_ab_k1[1] = (phi_2 * i2_ab[1]) + gama_2 * (vcap_ab[1] - Vg_ab[1]);

         //-------------------------------COM AMORTECIMENTO------------------------------------------------------------
         vcap_ab_k1[0] = (phi_3 * vcap_ab[0]) + gama_3 * (i1_ab[0] - i2_ab[0]);
         vcap_ab_k1[1] = (phi_3 * vcap_ab[1]) + gama_3 * (i1_ab[1] - i2_ab[1]);


 //---------------------------                 Estimando os valores em (K+2) do preditivo----------------------------------------------------------

         for (linha = 0; linha < 8; linha++)
         {

             //----------------------------------------Tensão de controle-------------------------------------------------------
             Vk_ab[0] = Vdc * s_ab[linha][0];
             Vk_ab[1] = Vdc * s_ab[linha][1];

             //---------------------------------------- Predição (k+2) -------------------------------------------------------------

             i1_ab_k2[0] = (phi_1 * i1_ab_k1[0]) + gama_1 * (Vk_ab[0] - vcap_ab_k1[0]);
             i1_ab_k2[1] = (phi_1 * i1_ab_k1[1]) + gama_1 * (Vk_ab[1] - vcap_ab_k1[1]);

             i2_ab_k2[0] = (phi_2 * i2_ab_k1[0]) + gama_2 * (vcap_ab_k1[0] - Vg_ab_k1[0]);
             i2_ab_k2[1] = (phi_2 * i2_ab_k1[1]) + gama_2 * (vcap_ab_k1[1] - Vg_ab_k1[1]);

             //--------------------------------------COM AMORTECIMENTO-----------------------------------------------------
             vcap_ab_k2[0] = (phi_3 * vcap_ab_k1[0]) + gama_3 * (i1_ab_k1[0] - i2_ab_k1[0]);
             vcap_ab_k2[1] = (phi_3 * vcap_ab_k1[1]) + gama_3 * (i1_ab_k1[1] - i2_ab_k1[1]);

//--------------------------------- calculando as variaveis do custo para ser comparado no preditivo ----------------------------------
             g2_alfa = (i2_ab_k2[0] - i2_ref_ab_k2[0]) * (i2_ab_k2[0] - i2_ref_ab_k2[0]);
             g2_beta = (i2_ab_k2[1] - i2_ref_ab_k2[1]) * (i2_ab_k2[1] - i2_ref_ab_k2[1]);

             g1_alfa = (i1_ab_k2[0] - i1_ref_ab_k2[0]) * (i1_ab_k2[0] - i1_ref_ab_k2[0]);
             g1_beta = (i1_ab_k2[1] - i1_ref_ab_k2[1]) * (i1_ab_k2[1] - i1_ref_ab_k2[1]);

             gcap_alfa = (vcap_ab_k2[0] - vcap_ref_ab_k2[0]) * (vcap_ab_k2[0] - vcap_ref_ab_k2[0]);
             gcap_beta = (vcap_ab_k2[1] - vcap_ref_ab_k2[1]) * (vcap_ab_k2[1] - vcap_ref_ab_k2[1]);

             gsa = (s_abc[linha][0] - s_abc[passado][0]) * (s_abc[linha][0] - s_abc[passado][0]);
             gsb = (s_abc[linha][1] - s_abc[passado][1]) * (s_abc[linha][1] - s_abc[passado][1]);
             gsc = (s_abc[linha][2] - s_abc[passado][2]) * (s_abc[linha][2] - s_abc[passado][2]);

             g_k = g_custo_i1 * (g1_alfa + g1_beta) + g_custo_vcap * (gcap_alfa + gcap_beta) + g_custo_i2 * (g2_alfa + g2_beta) +  g_custo_s * ((float) (gsa + gsb + gsc));


             if (g_k < g_op)
             {

                 linha_op[0] = linha;  // valor atual
                 g_op = g_k;
             }

         }

        volatile uint32_t pwm1, pwm2, pwm3;
        pwm1= s_abc[linha_op[0]][0];
        pwm2 = s_abc[linha_op[0]][1];
        pwm3= s_abc[linha_op[0]][2];

//--------------------------------------------------------
        /*
         if (g_trip_clear)
         {
         if ((EPWM_getTripZoneFlagStatus(EPWM0_BASE) & EPWM_TZ_FLAG_OST) != 0U)
         {
         EPWM_clearTripZoneFlag(EPWM0_BASE,EPWM_TZ_INTERRUPT | EPWM_TZ_FLAG_OST | EPWM_TZ_FLAG_DCAEVT1);

         }
         //g_trip_clear  = 0;
         }
         */
        //  }
    }
}
// Interrupção externa (XINT1 ou outro XINT ligado ao GPIO que recebe o PWM)
__interrupt void INT_myGPIO0_XINT_ISR(void)
{
    g_switch_on = GPIO_readPin(myGPIO0);

    Interrupt_clearACKGroup(INTERRUPT_ACK_GROUP1);
}

__interrupt void INT_myCPUTIMER0_ISR(void)
{
    // Atualiza contador
    g_step_counter++;

    // Reinicia no fim do ciclo PWM
    if (g_step_counter >= N_STEPS_PER_CYCLE)
        g_step_counter = 0;

    // Sinaliza para o loop principal que deve simular o próximo passo
    g_new_step_ready = true;

    // Libera nova interrupção
    Interrupt_clearACKGroup(INT_myCPUTIMER0_INTERRUPT_ACK_GROUP);
}


