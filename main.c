#include "driverlib.h"
#include "device.h"
#include "board.h"
#include "math.h"
#include "transformada_clark.h"
#include "funcoes.h"

//------------------------------------------INICIO DO CONTROLE-------------------------------------------------------------------------

volatile bool adc_flag = false;      // Flag para novo passo de simulação

// --------------------------------------- Configurações de conversão ------------------------------------------------------------------
// -----------------------------------
#define ADC_RESOLUTION     4095.0f
#define VREF_ADC     2.5f           // referência real do ADC
#define DIV_FACTOR 100.0f
#define DAC_RESOLUTION  4095.0f     // 12 bits
#define VREF_DAC        2.5f        // referência real do DAC (ou 3.3f)
#define VMAX_SINAL   250.0f      // pico do sinal real (ex: ±250V)


#define CALIB_TENSAO 0.92f      // CALIBRADOR PARA 200v
#define CALIB_CORRENTE 0.97f    // calibrador para 30A

#define norm_ADC  (84.0f)/4095.0F
#define norm_DAC 4095.0f/(84.0f)
// -------------------------------- Constantes pré-calculadas para transformada de clark -----------------------------------------------

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


// -------------------------------- valores definidos para o filtro SOGI com k=1 --------------------------------------------------------

#define INV_TS  40000.0f    // 1/Ts Ts= 0.000025;
#define Lc 0.00584f
#define Lg 0.00106f
#define rc 0.2f
#define rg 0.17f
#define C 0.0000114f

//----------------------------------- Vetor chaveamento alfa/beta e ABC ------------------------------------------------------------------

#define ZERO 0
#define UM 1

static const float s_ab[8][2] =     { {ZERO,          ZERO },
                                    {DOIS_DIV_3,    ZERO },
                                    {UM_DIV_3,      SQRT3_BY_3 },
                                    {-UM_DIV_3,     SQRT3_BY_3 },
                                    {-DOIS_DIV_3,   ZERO },
                                    {-UM_DIV_3,     -SQRT3_BY_3 },
                                    {UM_DIV_3,      -SQRT3_BY_3 },
                                    {ZERO,          ZERO } };

static const int16_t s_abc[8][3] = { {ZERO,     ZERO,       ZERO },
                                     {UM,       ZERO,       ZERO },
                                     {UM,       UM,         ZERO },
                                     {ZERO,     UM,         ZERO },
                                     {ZERO,     UM,         UM },
                                     {ZERO,     ZERO,       UM },
                                     {UM,       ZERO,       UM },
                                     {UM,       UM,         UM } };

// ---------------------------------- custo para utilizar no controle preditivo ------------------------------------------------------

#define g_custo_i1 1.0f
#define g_custo_i2 1.0f
#define g_custo_vcap 1.0f
#define g_custo_s 1.0f

// ---------------------------------- constantes utilizadas para predizer amotras a frente -------------------------------------------

#define Const 0.602655f          // Const=((1/RV)+(C/Ts))
#define const_5 0.456f           // const_5=(C/Ts)
#define phi_1 0.999145f          // phi_1 = 1 - ((rc*Ts)/Lc)
#define gama_1 0.00428082f      //  gama_1 = Ts/Lc
#define phi_2 0.995991f          // phi_2 = 1 - ((rg*Ts)/Lg)
#define gama_2 0.0235849f        // gama_2 = Ts/Lg
#define phi_3 0.6785f            // phi_3 = 1 - (Ts/(C*RV))
#define gama_3 2.19298f          // gama_3 = Ts/C


// --------------------------------- valores definidos para o filtro SOGI com k=1 ----------------------------------------------------
#define a11 0.00471238898f
#define a22 0.00471238898f
#define b1  1.990531226413290f
#define b2  0.990619634277416f
#define a1  0.004690182861292f
#define a2  0.004690182861292f

//-------------------------------  variaveis LIDAS FILTRO E REDE ----------------------------------------------------------------------

volatile float va_pac, vb_pac, vc_pac;
volatile float i1_a, i1_b, i1_c;
volatile float va_cap, vb_cap, vc_cap;
volatile float i2_a, i2_b, i2_c;
volatile float Vdc;

// ------------------------------- controle e referencia:  K=0 K=-1 K=1--------------------------------------------------------------------
//       K=0 CONTROLE CORRENTE BALANCEADA
//       K=1 CONTROLE POTENCIA ATIVA
//       K=-1 CONTROLE POTENCIA REATIVA

int16_t linha_op[2] = { 0, 0 };
int16_t passado = 0;
volatile int16_t K = 0;
volatile float P_ref = 5000.0f;
volatile float Q_ref = 0.0f;
volatile float P_ativa, Q_reativa;

// -------------------------------- Vetores para armazenar posições -------------------------------------------------------------------------

volatile float Vg_ab[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i1_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float i2_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float vcap_ab[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
volatile float Vg_ab_k1[2] = { 0.0f, 0.0f };
volatile float Vg_ab_k2[2] = { 0.0f, 0.0f };
volatile float Vg_ab_filtrado[6] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
volatile float Vg_ab_q_filtrado[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
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
volatile float Vk_ab[2];
#define INICIO_OPERACIONAL 1e15f

// --------------------------------------- funçoes esecificas ara cada parte do controle ---------------------------------------------------

void redefinindo_vetores_alfa_beta(void);
void aplicarFiltroSOGI(void);
void calcularSequencias(void);
void gerarCorrenteReferenciaI2(void);
void gerarCorrenteVirtual(void);
void calcularPotencias(void);
void gerarReferenciaVcap(void);
void gerarReferenciaI1(void);
void estimarValores_k1(void);
int16_t calcularLinhaOtimizada(void);
void aplicarPWM(int16_t linha);
uint16_t valor_para_DAC(float valor_volts);
static inline float adc_to_volts(uint16_t adc_code);

//---------------------------------------------- VARIAVEIS DE LEITURA ADC ANTES DE CONVERTER ------------------------------------------------

uint16_t g_i1_a, g_i1_b, g_i1_c, g_i2_a, g_i2_b, g_i2_c, g_cap_a, g_cap_b, g_cap_c, g_vcc;



float corrente = 35.0f, tensao = 200.0f;           // tensão desejada

uint16_t g_switch_on,g_switch_off;
uint32_t pwm1, pwm2, pwm3;


float valor_envio = 123.4f;  // coloque aqui qualquer valor entre -250 e +250
uint16_t codigo_dac;
volatile float valor_retorno,tempo_us;
uint16_t corrente_corrigida;
uint16_t tensao_corrigida;


uint32_t inicio, fim;
volatile float periodo_us = 0.0f;
volatile uint32_t ultimo_contador = 0;
#define TIMER_PERIOD_TICKS 5000.0f
volatile uint32_t ultimo_adc_ticks = 0;    // Guarda o contador da última ISR
volatile float periodo_adc_us = 0.0f;      // Guarda o período em microssegundos

#define MAX_SAMPLES 1000       // número máximo de registros
volatile uint8_t toggle_buffer[MAX_SAMPLES]; // guarda 0 ou 1
volatile uint16_t toggle_index = 0;          // índice atual
volatile bool adc_toggle = false;           // toggle atual

void main(void)
{

    // Inicialização dos periféricos

    Device_init();
    Interrupt_initModule();
    Interrupt_initVectorTable();
    Board_init();

    EINT;
    ERTM;


    while (1)
    {
        if (adc_flag)
        {

            adc_flag = false;


            // Alterna o toggle
            adc_toggle = !adc_toggle;

            // Salva o estado do toggle no vetor (0 ou 1)
            toggle_buffer[toggle_index] = adc_toggle ? 1 : 0;
            toggle_index++;

            // Se chegar no final, volta para o início (ou para de registrar)
            if (toggle_index >= MAX_SAMPLES)
            {
                toggle_index = 0; // sobrescreve os dados antigos (circular)
                // ou, se quiser parar: toggle_index = MAX_SAMPLES - 1;
            }


//-------------teste envio DAC---------------------------------------------------------------------------------

            corrente_corrigida = valor_para_DAC(corrente);  // apenas calcula o código DAC
            DAC_setShadowValue(DAC0_BASE, corrente_corrigida);

            tensao_corrigida = valor_para_DAC(tensao);  // apenas calcula o código DAC
            DAC_setShadowValue(DAC1_BASE, tensao_corrigida);

            codigo_dac = valor_para_DAC(valor_envio);
            valor_retorno = adc_to_volts(codigo_dac);
//------------------------------------------------------------------------------------------------------------


                  // Conversão para tensão real ADC

                  // --------------------------------------
                 i1_a = adc_to_volts(g_i1_a);
                 i1_b = adc_to_volts(g_i1_b);
                 i1_c = adc_to_volts(g_i1_c);


                 i2_a = adc_to_volts(g_i2_a);
                 i2_b = adc_to_volts(g_i2_b);
                 i2_c = adc_to_volts(g_i2_c);

                 va_cap = adc_to_volts(g_cap_a);
                 vb_cap = adc_to_volts(g_cap_b);
                 vc_cap = adc_to_volts(g_cap_c);

                 // tenho que colocar depois para Vdc
                 Vdc = adc_to_volts(g_vcc);

 //--------------------------------------------------------------------------

             passado = linha_op[0];
             linha_op[1] = linha_op[0];
             linha_op[0] = 1000;
             //  g_op = INICIO_OPERACIONAL;

             redefinindo_vetores_alfa_beta();
             aplicarFiltroSOGI();
             calcularSequencias();
             gerarCorrenteReferenciaI2();
             gerarCorrenteVirtual();
             calcularPotencias();
             gerarReferenciaVcap();
             gerarReferenciaI1();
             estimarValores_k1();
             linha_op[0] = calcularLinhaOtimizada();
             aplicarPWM(linha_op[0]);

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
        }

    }
}

void redefinindo_vetores_alfa_beta(void)
{

    // Componentes alfa-beta
    float valfa_pac, vbeta_pac;
    float i1_alfa, i1_beta;
    float valfa_cap, vbeta_cap;
    float i2_alfa, i2_beta;
// ========================================== 1. Transformadas CLARK ==========================================================================

    ClarkeTransform(va_pac, vb_pac, vc_pac, &valfa_pac, &vbeta_pac);
    ClarkeTransform(i1_a, i1_b, i1_c, &i1_alfa, &i1_beta);
    ClarkeTransform(va_cap, vb_cap, vc_cap, &valfa_cap, &vbeta_cap);
    ClarkeTransform(i2_a, i2_b, i2_c, &i2_alfa, &i2_beta);

// ========================================== 2. Redefinindo as variaveis para vetores =========================================================
    // operador 0 é para PAR
    // operador 1 é para IMPAR

    atualizar_valor_par_impar_vetor_6(Vg_ab, 0, Vg_ab[0]);
    atualizar_valor_par_impar_vetor_6(Vg_ab, 1, Vg_ab[1]);
    Vg_ab[0] = valfa_pac;
    Vg_ab[1] = vbeta_pac;

    extrapolar_k1(Vg_ab, Vg_ab_k1);
    extrapolar_k2(Vg_ab, Vg_ab_k2);


    atualizar_valor_par_impar_vetor_4(vcap_ab, 0, vcap_ab[0]);
    atualizar_valor_par_impar_vetor_4(vcap_ab, 1, vcap_ab[1]);
    vcap_ab[0] = valfa_cap;
    vcap_ab[1] = vbeta_cap;

    atualizar_valor_par_impar_vetor_4(i1_ab, 0, i1_ab[0]);
    atualizar_valor_par_impar_vetor_4(i1_ab, 1, i1_ab[1]);
    i1_ab[0] = i1_alfa;
    i1_ab[1] = i1_beta;

    atualizar_valor_par_impar_vetor_4(i2_ab, 0, i2_ab[0]);
    atualizar_valor_par_impar_vetor_4(i2_ab, 1, i2_ab[1]);
    i2_ab[0] = i2_alfa;
    i2_ab[1] = i2_beta;
}

void aplicarFiltroSOGI(void)
{

// =========================================== 4. Filtro SOGI ============================================================================

    atualizar_valor_par_impar_vetor_6(Vg_ab_filtrado, 0, Vg_ab_filtrado[0]);
    atualizar_valor_par_impar_vetor_6(Vg_ab_filtrado, 1, Vg_ab_filtrado[1]);
    Vg_ab_filtrado[0] = a1 * Vg_ab[0] - a2 * Vg_ab[4] + b1 * Vg_ab_filtrado[2] - b2 * Vg_ab_filtrado[4];
    Vg_ab_filtrado[1] = a1 * Vg_ab[1] - a2 * Vg_ab[5] + b1 * Vg_ab_filtrado[3] - b2 * Vg_ab_filtrado[5];

    //   componente em quadratura

    atualizar_valor_par_impar_vetor_4(Vg_ab_q_filtrado, 0, Vg_ab_q_filtrado[0]);
    atualizar_valor_par_impar_vetor_4(Vg_ab_q_filtrado, 1, Vg_ab_q_filtrado[1]);
    Vg_ab_q_filtrado[0] = a11 * Vg_ab_filtrado[0] + a22 * Vg_ab_filtrado[2] + Vg_ab_q_filtrado[2];
    Vg_ab_q_filtrado[1] = a11 * Vg_ab_filtrado[1] + a22 * Vg_ab_filtrado[3] + Vg_ab_q_filtrado[3];

}

void calcularSequencias(void)
{

// =========================================== 5. Sequência positiva e negativa =============================================================

    Vg_ab_filtrado_pos[0] = HALF * (Vg_ab_filtrado[0] - Vg_ab_q_filtrado[1]);
    Vg_ab_filtrado_pos[1] = HALF * (Vg_ab_filtrado[1] + Vg_ab_q_filtrado[0]);

    Vg_ab_filtrado_neg[0] = HALF * (Vg_ab_filtrado[0] + Vg_ab_q_filtrado[1]);
    Vg_ab_filtrado_neg[1] = HALF * (Vg_ab_filtrado[1] - Vg_ab_q_filtrado[0]);

    if (Vg_ab_filtrado_pos[0] == ZERO)
        Vg_ab_filtrado_pos[0] = LIMITAR_ZERO;

    if (Vg_ab_filtrado_pos[1] == ZERO)
        Vg_ab_filtrado_pos[1] = LIMITAR_ZERO;

}

void gerarCorrenteReferenciaI2(void)
{

    float Vg_ab_filtrado_pos_0_quadrado, Vg_ab_filtrado_pos_1_quadrado;
    float Vg_ab_filtrado_neg_0_quadrado, Vg_ab_filtrado_neg_1_quadrado;
    float i2_ref_alfa_pos, i2_ref_alfa_neg;
    float i2_ref_beta_pos, i2_ref_beta_neg;

    // Ganhos, potências, controle
    float D, E, inv_D, inv_E;

    // =========================================== 7. Utilizando as sequencia para determinar a referencia conforme i2 ===============================


    Vg_ab_filtrado_pos_0_quadrado = Vg_ab_filtrado_pos[0] * Vg_ab_filtrado_pos[0];
    Vg_ab_filtrado_pos_1_quadrado = Vg_ab_filtrado_pos[1] * Vg_ab_filtrado_pos[1];
    Vg_ab_filtrado_neg_0_quadrado = Vg_ab_filtrado_neg[0] * Vg_ab_filtrado_neg[0];
    Vg_ab_filtrado_neg_1_quadrado = Vg_ab_filtrado_neg[1] * Vg_ab_filtrado_neg[1];

    D = (Vg_ab_filtrado_pos_0_quadrado + Vg_ab_filtrado_pos_1_quadrado) - K * (Vg_ab_filtrado_neg_0_quadrado + Vg_ab_filtrado_neg_1_quadrado);
    E = (Vg_ab_filtrado_pos_0_quadrado + Vg_ab_filtrado_pos_1_quadrado) + K * (Vg_ab_filtrado_neg_0_quadrado + Vg_ab_filtrado_neg_1_quadrado);

    // Proteger contra divisão por zero
    if (D == ZERO)
        D = LIMITAR_ZERO;
    if (E == ZERO)
        E = LIMITAR_ZERO;

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
}

void gerarCorrenteVirtual(void)
{

    // =========================================== 8. Calculo corrente virtual ===================================================================

    //foi realizado uma media entre a tensao medida no capacitor e tensao referencia capacitor
    atualizar_valor_par_impar_vetor_6(i2_ref_ab_virtual, 0, i2_ref_ab_virtual[0]);
    atualizar_valor_par_impar_vetor_6(i2_ref_ab_virtual, 1, i2_ref_ab_virtual[1]);

    i2_ref_ab_virtual[0] = i2_ref_ab[0] - (1.55f * vcap_ref_ab[0] + 0.45f * vcap_ab[0]) * INV_2_RV;
    i2_ref_ab_virtual[1] = i2_ref_ab[1] - (1.55f * vcap_ref_ab[1] + 0.45f * vcap_ab[1]) * INV_2_RV;
}

void calcularPotencias(void)
{

    // =========================================== 9. Calculo potencia ativa e reativa real  ===================================================================

    P_ativa = TRES_DIV_2 * (Vg_ab[0] * i2_ab[0] + Vg_ab[1] * i2_ab[1]);

    Q_reativa = TRES_DIV_2 * (Vg_ab[1] * i2_ab[0] - Vg_ab[0] * i2_ab[1]);
}

void gerarReferenciaVcap(void)
{

// =========================================== 10. Gerando Ref - Vcap  ===================================================================

    derivada_i2[0] = (i2_ref_ab_virtual[0] - i2_ref_ab_virtual[2]) * INV_TS;
    derivada_i2[1] = (i2_ref_ab_virtual[1] - i2_ref_ab_virtual[3]) * INV_TS;

    //  Limitando a derivada para remover picos na transição

    derivada_i2[0] = (derivada_i2[0] > LIMITADOR_UM)  ?  LIMITADOR_UM : derivada_i2[0];
    derivada_i2[0] = (derivada_i2[0] < -LIMITADOR_UM) ? -LIMITADOR_UM : derivada_i2[0];
    derivada_i2[1] = (derivada_i2[1] > LIMITADOR_UM)  ?  LIMITADOR_UM : derivada_i2[1];
    derivada_i2[1] = (derivada_i2[1] < -LIMITADOR_UM) ? -LIMITADOR_UM : derivada_i2[1];

    atualizar_valor_par_impar_vetor_6(vcap_ref_ab, 0, vcap_ref_ab[0]);
    atualizar_valor_par_impar_vetor_6(vcap_ref_ab, 1, vcap_ref_ab[1]);
    vcap_ref_ab[0] = Vg_ab[0] + Lg * derivada_i2[0] + rg * i2_ref_ab_virtual[0];
    vcap_ref_ab[1] = Vg_ab[1] + Lg * derivada_i2[1] + rg * i2_ref_ab_virtual[1];

    //extrapolação da tensão no capacitor
    extrapolar_k2(vcap_ref_ab, vcap_ref_ab_k2);
}

void gerarReferenciaI1(void)
{

    // =========================================== 11. Gerando Ref da corrente i1 com resistor em paralelo com o capacitor  ===================================================================
    //----------------------------------------       Amortecimento do filtro LCL---------------------------------------------

    atualizar_valor_par_impar_vetor_6(i1_ref_ab, 0, i1_ref_ab[0]);
    atualizar_valor_par_impar_vetor_6(i1_ref_ab, 1, i1_ref_ab[1]);
    i1_ref_ab[0] = const_5 * (vcap_ref_ab[0] - vcap_ref_ab[2]) + (vcap_ref_ab[0] * INV_RV) + i2_ref_ab_virtual[0];
    i1_ref_ab[1] = const_5 * (vcap_ref_ab[1] - vcap_ref_ab[3]) + (vcap_ref_ab[1] * INV_RV) + i2_ref_ab_virtual[1];

    //extrapolação da corrente i1
    extrapolar_k2(i1_ref_ab, i1_ref_ab_k2);
}

void estimarValores_k1(void)
{

// =========================================== 12. Entrada do algoritmo de controle  ===================================================================

    //-------------------------------------------- Estimando os valores em (K+1) do preditivo----------------------------------------------------------

    Vk_ab[0] = Vdc * s_ab[passado][0];
    Vk_ab[1] = Vdc * s_ab[passado][1];

    i1_ab_k1[0] = (phi_1 * i1_ab[0]) + gama_1 * (Vk_ab[0] - vcap_ab[0]);
    i1_ab_k1[1] = (phi_1 * i1_ab[1]) + gama_1 * (Vk_ab[1] - vcap_ab[1]);

    i2_ab_k1[0] = (phi_2 * i2_ab[0]) + gama_2 * (vcap_ab[0] - Vg_ab[0]);
    i2_ab_k1[1] = (phi_2 * i2_ab[1]) + gama_2 * (vcap_ab[1] - Vg_ab[1]);

    //-------------------------------COM AMORTECIMENTO------------------------------------------------------------
    vcap_ab_k1[0] = (phi_3 * vcap_ab[0]) + gama_3 * (i1_ab[0] - i2_ab[0]);
    vcap_ab_k1[1] = (phi_3 * vcap_ab[1]) + gama_3 * (i1_ab[1] - i2_ab[1]);

}

int16_t calcularLinhaOtimizada(void)
{

    float g_op;
    int16_t melhor_linha = 0;
    float gcap_alfa, gcap_beta;
    float g2_alfa, g2_beta, g1_alfa, g1_beta;
    int16_t gsa, gsb, gsc;
    float g_k;
    int16_t linha;

    g_op = INICIO_OPERACIONAL;

    //------------------------------------------- Estimando os valores em (K+2) do preditivo----------------------------------------------------------

    for (linha = 0; linha < 8; linha++)
    {

        //--------------------------------------- Tensão de controle-------------------------------------------------------
        Vk_ab[0] = Vdc * s_ab[linha][0];
        Vk_ab[1] = Vdc * s_ab[linha][1];

        //---------------------------------------- Predição (k+2) ----------------------------------------------------------

        i1_ab_k2[0] = (phi_1 * i1_ab_k1[0]) + gama_1 * (Vk_ab[0] - vcap_ab_k1[0]);
        i1_ab_k2[1] = (phi_1 * i1_ab_k1[1]) + gama_1 * (Vk_ab[1] - vcap_ab_k1[1]);

        i2_ab_k2[0] = (phi_2 * i2_ab_k1[0]) + gama_2 * (vcap_ab_k1[0] - Vg_ab_k1[0]);
        i2_ab_k2[1] = (phi_2 * i2_ab_k1[1]) + gama_2 * (vcap_ab_k1[1] - Vg_ab_k1[1]);

        //-------------------------------------- COM AMORTECIMENTO------------------------------------------------------------

        vcap_ab_k2[0] = (phi_3 * vcap_ab_k1[0]) + gama_3 * (i1_ab_k1[0] - i2_ab_k1[0]);
        vcap_ab_k2[1] = (phi_3 * vcap_ab_k1[1]) + gama_3 * (i1_ab_k1[1] - i2_ab_k1[1]);

        //--------------------------------- calculando as variaveis do custo para ser comparado no preditivo ------------------

        g2_alfa = (i2_ab_k2[0] - i2_ref_ab_k2[0]) * (i2_ab_k2[0] - i2_ref_ab_k2[0]);
        g2_beta = (i2_ab_k2[1] - i2_ref_ab_k2[1]) * (i2_ab_k2[1] - i2_ref_ab_k2[1]);

        g1_alfa = (i1_ab_k2[0] - i1_ref_ab_k2[0]) * (i1_ab_k2[0] - i1_ref_ab_k2[0]);
        g1_beta = (i1_ab_k2[1] - i1_ref_ab_k2[1]) * (i1_ab_k2[1] - i1_ref_ab_k2[1]);

        gcap_alfa = (vcap_ab_k2[0] - vcap_ref_ab_k2[0]) * (vcap_ab_k2[0] - vcap_ref_ab_k2[0]);
        gcap_beta = (vcap_ab_k2[1] - vcap_ref_ab_k2[1]) * (vcap_ab_k2[1] - vcap_ref_ab_k2[1]);

        gsa = (s_abc[linha][0] - s_abc[passado][0]) * (s_abc[linha][0] - s_abc[passado][0]);
        gsb = (s_abc[linha][1] - s_abc[passado][1]) * (s_abc[linha][1] - s_abc[passado][1]);
        gsc = (s_abc[linha][2] - s_abc[passado][2]) * (s_abc[linha][2] - s_abc[passado][2]);

        g_k = g_custo_i1 * (g1_alfa + g1_beta)
              + g_custo_vcap * (gcap_alfa + gcap_beta)
              + g_custo_i2 * (g2_alfa + g2_beta)
              + g_custo_s * ((float) (gsa + gsb + gsc));

        if (g_k < g_op)
        {

            melhor_linha = linha;  // valor atual
            g_op = g_k;
        }

    }

    return melhor_linha;
}

void aplicarPWM(int16_t linha)
{

 //   uint32_t pwm1, pwm2, pwm3;

    pwm1 = s_abc[linha][0]; //braço A
    pwm2 = s_abc[linha][1]; //braço B
    pwm3 = s_abc[linha][2]; //braço C

    // Braço A
    EPWM_setActionQualifierContSWForceAction( myEPWM1_BASE, EPWM_AQ_OUTPUT_A, (EPWM_ActionQualifierSWOutput)
                                                                              (pwm1 ? EPWM_AQ_OUTPUT_HIGH : EPWM_AQ_OUTPUT_LOW));

    EPWM_setActionQualifierContSWForceAction( myEPWM2_BASE, EPWM_AQ_OUTPUT_A, (EPWM_ActionQualifierSWOutput)
                                                                              (pwm2 ? EPWM_AQ_OUTPUT_HIGH : EPWM_AQ_OUTPUT_LOW));

    EPWM_setActionQualifierContSWForceAction( myEPWM3_BASE, EPWM_AQ_OUTPUT_A, (EPWM_ActionQualifierSWOutput)
                                                                              (pwm3 ? EPWM_AQ_OUTPUT_HIGH : EPWM_AQ_OUTPUT_LOW));

}

//---------------- Conversão ADC → Volts ----------------
static inline float adc_to_volts(uint16_t adc_code)
{
    /*
    // Converte ADC → tensão no pino
    float v_adc_pin = ((float)adc_code / ADC_RESOLUTION) * VREF_ADC;

    // Converte para valor real, considerando mapeamento bipolar e divisor
    float valor_volts = ((v_adc_pin / (VREF_ADC / 2.0f)) - 1.0f)*VMAX_SINAL;

    return valor_volts;
*/
    // Corrige offset e aplica ganho em contagem ADC
      int adc_offset = 185;
      int adc_corr = adc_code - adc_offset;

      if (adc_corr >= 2048) { // faixa positiva
          float ganho_pos = 0.857f;
          adc_corr = 2048 + (adc_corr - 2048) * ganho_pos;
      } else { // faixa negativa
          float ganho_neg = 0.766f;
          adc_corr = 2048 + (adc_corr - 2048) * ganho_neg;
      }

      // Converte para volts depois
      float v_adc_pin = ((float)adc_corr / ADC_RESOLUTION) * VREF_ADC;
      float valor_volts = ((v_adc_pin / (VREF_ADC / 2.0f)) - 1.0f) * VMAX_SINAL;

      return valor_volts;

}

//---------------- Conversão Volts → DAC ----------------
uint16_t valor_para_DAC(float valor_volts)
{
    // Limita tensão ±VMAX
    if (valor_volts >  VMAX_SINAL) valor_volts =  VMAX_SINAL;
    if (valor_volts < -VMAX_SINAL) valor_volts = -VMAX_SINAL;

    // Aplica divisor físico
    float v_dac_pin = (valor_volts / VMAX_SINAL + 1.0f) * (VREF_DAC / 2.0f);

    // Código DAC
    uint16_t dac_code = (uint16_t)((v_dac_pin / VREF_DAC) * DAC_RESOLUTION + 0.5f);

    if (dac_code > DAC_RESOLUTION)
        dac_code = (uint16_t)DAC_RESOLUTION;

    return dac_code;
}


// Interrupção externa (XINT1 ou outro XINT ligado ao GPIO que recebe o PWM)
__interrupt void INT_myGPIO0_XINT_ISR(void)
{

    g_switch_on = GPIO_readPin(myGPIO0);
    g_switch_off = GPIO_readPin(myGPIO1);

    Interrupt_clearACKGroup(INTERRUPT_ACK_GROUP1);
}

__interrupt void INT_myCPUTIMER0_ISR(void)
{
/*
uint32_t contador_atual = CPUTimer_getTimerCount(CPUTIMER0_BASE);
if (ultimo_contador != 0)
{
    int32_t delta;
    if (contador_atual > ultimo_contador)
        delta = (int32_t)(TIMER_PERIOD_TICKS - (contador_atual - ultimo_contador));
    else
        delta = (int32_t)(ultimo_contador - contador_atual);

    periodo_us = (float)delta / 200.0f; // CPU = 200 MHz
}
ultimo_contador = contador_atual;

        // Limpa flag da interrupção
        CPUTimer_clearOverflowFlag(CPUTIMER0_BASE);
        Interrupt_clearACKGroup(INTERRUPT_ACK_GROUP1);
*/
    // Atualiza contador
//    g_step_counter++;

    // Reinicia no fim do ciclo PWM
//    if (g_step_counter >= N_STEPS_PER_CYCLE)
//        g_step_counter = 0;

    // Sinaliza para o loop principal que deve simular o próximo passo
//    g_new_step_ready = true;

 //   enviar_vetor_para_DAC();
//   CPUTimer_clearOverflowFlag(myCPUTIMER0_BASE);
    // Libera nova interrupção
    Interrupt_clearACKGroup(INT_myCPUTIMER0_INTERRUPT_ACK_GROUP);
}

// será realizado a leitura de todas as variaveis ADC

__interrupt void INT_ADC_C_1_ISR(void)
{
    uint32_t atual = CPUTimer_getTimerCount(CPUTIMER0_BASE);
    uint32_t delta_ticks;

    // Calcula o tempo entre esta interrupção e a anterior
    if (ultimo_adc_ticks != 0)
    {
        if (atual <= ultimo_adc_ticks)
            delta_ticks = ultimo_adc_ticks - atual;
        else
            delta_ticks = TIMER_PERIOD_TICKS - (atual - ultimo_adc_ticks); // corrige overflow

        // Converte ticks para microssegundos
        periodo_adc_us = (float)delta_ticks / 200.0f;   // CPU = 200 MHz → 1 tick = 5 ns
    }

    ultimo_adc_ticks = atual; // Atualiza para a próxima medição




 //   inicio = CPUTimer_getTimerCount(CPUTIMER0_BASE);

    // ADCA
    g_i1_a = ADC_readResult(ADC_A_RESULT_BASE, ADC_A_SOC0);
    g_i1_b = ADC_readResult(ADC_A_RESULT_BASE, ADC_A_SOC1);
    g_i1_c = ADC_readResult(ADC_A_RESULT_BASE, ADC_A_SOC2);
    g_vcc = ADC_readResult(ADC_A_RESULT_BASE, ADC_A_SOC3);

    // ADCB
    g_i2_a = ADC_readResult(ADC_B_RESULT_BASE, ADC_B_SOC4);
    g_i2_b = ADC_readResult(ADC_B_RESULT_BASE, ADC_B_SOC5);
    g_i2_c = ADC_readResult(ADC_B_RESULT_BASE, ADC_B_SOC6);

    // ADCC
    g_cap_a = ADC_readResult(ADC_C_RESULT_BASE, ADC_C_SOC7);
    g_cap_b = ADC_readResult(ADC_C_RESULT_BASE, ADC_C_SOC8);
    g_cap_c = ADC_readResult(ADC_C_RESULT_BASE, ADC_C_SOC9);

    adc_flag = true;

    // Limpa a interrupção
    ADC_clearInterruptStatus(ADC_C_BASE, ADC_INT_NUMBER1);
    Interrupt_clearACKGroup(INT_ADC_C_1_INTERRUPT_ACK_GROUP);

}

