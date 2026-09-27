/*
    ==========================================================================
    TRABALHO 1 - RESOLUCAO DE SISTEMA NAO LINEAR (METODO DO GRADIENTE)
    ==========================================================================

    Sistema a resolver (3 equacoes, 3 incognitas x, y, z):

        f1(x,y,z) = (x-1)^2 + (y-1)^2 + (z-1)^2 - 1   = 0
        f2(x,y,z) = 2*x^2 + (y-1)^2 - 4*z             = 0
        f3(x,y,z) = 3*x^2 + 2*z^2 - 4*y               = 0

    Ideia do metodo:
        Transformamos o sistema em um problema de MINIMIZACAO, definindo a
        funcao objetivo (soma dos quadrados dos residuos):

            G(x,y,z) = f1^2 + f2^2 + f3^2

        Se conseguirmos encontrar (x,y,z) tal que G = 0, entao f1 = f2 = f3 = 0
        e o sistema original esta resolvido.

        O METODO DO GRADIENTE (steepest descent) caminha, a cada iteracao,
        no sentido contrario ao gradiente de G (direcao de maior decrescimo):

            (x,y,z)_novo = (x,y,z)_velho - alpha * grad(G)

        O passo "alpha" e escolhido a cada iteracao por BUSCA LINEAR (regra de
        Armijo com backtracking): comeca em alpha = 1 e vai reduzindo pela
        metade ate que G realmente diminua o suficiente.

    Criterio de parada:
        ||grad(G)|| < TOL   (gradiente praticamente nulo -> ponto estacionario)
        ou G(x,y,z) < TOL2  (residuos praticamente nulos)
        ou numero maximo de iteracoes atingido.

    Codigo escrito de forma procedural (sem classes), seguindo o mesmo
    estilo dos demais programas do trabalho.

    Multiplos chutes iniciais:
        Como o metodo do gradiente so encontra minimos LOCAIS de G, o ponto
        de partida influencia para qual solucao o metodo converge. Por isso
        o programa testa o chute fixo (0.5, 0.5, 0.5) - que converge em 99
        iteracoes - junto com 10 chutes iniciais aleatorios (x,y,z cada um
        no intervalo [-2, 1]), resolve o sistema a partir de cada um deles
        e, no final, mostra apenas o resultado da tentativa que atingiu o
        menor valor de G (melhor convergencia).
    ==========================================================================
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <string>
#include <chrono>

// ---------------------------------------------------------------------
// CONSTANTES
// ---------------------------------------------------------------------
const int MAX_ITERACOES = 200000;
const int MAX_BACKTRACKING = 60;
const double TOL_GRADIENTE = 1e-10;   // ||grad G|| menor que isso -> parou
const double TOL_OBJETIVO = 1e-16;    // G menor que isso -> residuos ~ 0
const double C_ARMIJO = 1e-4;         // constante da condicao de Armijo
const int NUM_CHUTES = 10;            // quantidade de chutes iniciais aleatorios testados
const double CHUTE_MIN = -1.0;
const double CHUTE_MAX = 1.0;

// ---------------------------------------------------------------------
// FUNCAO: residuos
// Calcula os tres residuos f1, f2, f3 do sistema no ponto (x,y,z).
// ---------------------------------------------------------------------
void residuos(double x, double y, double z, double &f1, double &f2, double &f3)
{
    f1 = (x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1) - 1;
    f2 = 2 * x * x + (y - 1) * (y - 1) - 4 * z;
    f3 = 3 * x * x + 2 * z * z - 4 * y;
}

// ---------------------------------------------------------------------
// FUNCAO: funcaoObjetivo
// G(x,y,z) = f1^2 + f2^2 + f3^2  (queremos minimizar ate G ~ 0)
// ---------------------------------------------------------------------
double funcaoObjetivo(double x, double y, double z)
{
    double f1, f2, f3;
    residuos(x, y, z, f1, f2, f3);
    return f1 * f1 + f2 * f2 + f3 * f3;
}

// ---------------------------------------------------------------------
// FUNCAO: gradiente
// Calcula o gradiente analitico de G em (x,y,z):
//     dG/dv = 2*f1*(df1/dv) + 2*f2*(df2/dv) + 2*f3*(df3/dv),  v = x,y,z
// ---------------------------------------------------------------------
void gradiente(double x, double y, double z, double &gx, double &gy, double &gz)
{
    double f1, f2, f3;
    residuos(x, y, z, f1, f2, f3);

    // derivadas parciais de cada residuo
    double df1dx = 2 * (x - 1), df1dy = 2 * (y - 1), df1dz = 2 * (z - 1);
    double df2dx = 4 * x,       df2dy = 2 * (y - 1), df2dz = -4;
    double df3dx = 6 * x,       df3dy = -4,           df3dz = 4 * z;

    gx = 2 * (f1 * df1dx + f2 * df2dx + f3 * df3dx);
    gy = 2 * (f1 * df1dy + f2 * df2dy + f3 * df3dy);
    gz = 2 * (f1 * df1dz + f2 * df2dz + f3 * df3dz);
}

// ---------------------------------------------------------------------
// FUNCAO: resolverGradiente
// Executa o metodo do gradiente a partir do chute (x0,y0,z0). Devolve o
// ponto final, o numero de iteracoes gastas, se convergiu e o valor final
// de G (funcao objetivo), usado para comparar qual chute deu o melhor resultado.
// Se "arquivo" nao for nulo, grava os residuos f1,f2,f3 de cada iteracao nele.
// ---------------------------------------------------------------------
void resolverGradiente(double x0, double y0, double z0,
                        double &x, double &y, double &z,
                        int &iteracoes, bool &convergiu, double &gFinal,
                        FILE *arquivo = nullptr)
{
    x = x0;
    y = y0;
    z = z0;
    convergiu = false;

    if (arquivo != nullptr)
    {
        fprintf(arquivo, "\n----------------------------------------------------------\n");
        fprintf(arquivo, "Chute inicial: x0 = %.6f, y0 = %.6f, z0 = %.6f\n", x0, y0, z0);
        fprintf(arquivo, "----------------------------------------------------------\n");
        fprintf(arquivo, "%8s %16s %16s %16s %16s\n", "iter", "f1", "f2", "f3", "G");
    }

    for (iteracoes = 1; iteracoes <= MAX_ITERACOES; iteracoes++)
    {
        double gx, gy, gz;
        gradiente(x, y, z, gx, gy, gz);

        double normaGrad = sqrt(gx * gx + gy * gy + gz * gz);
        double gAtual = funcaoObjetivo(x, y, z);

        if (arquivo != nullptr)
        {
            double f1, f2, f3;
            residuos(x, y, z, f1, f2, f3);
            fprintf(arquivo, "%8d %16.6e %16.6e %16.6e %16.6e\n", iteracoes, f1, f2, f3, gAtual);
        }

        if (normaGrad < TOL_GRADIENTE || gAtual < TOL_OBJETIVO)
        {
            convergiu = true;
            gFinal = gAtual;
            return;
        }

        // busca linear (backtracking com condicao de Armijo)
        double alpha = 1.0;
        double xNovo, yNovo, zNovo, gNovo;
        int passoBacktracking;

        for (passoBacktracking = 0; passoBacktracking < MAX_BACKTRACKING; passoBacktracking++)
        {
            xNovo = x - alpha * gx;
            yNovo = y - alpha * gy;
            zNovo = z - alpha * gz;
            gNovo = funcaoObjetivo(xNovo, yNovo, zNovo);

            if (gNovo <= gAtual - C_ARMIJO * alpha * normaGrad * normaGrad)
                break; // decresceu o suficiente, aceita o passo

            alpha *= 0.5; // passo grande demais, reduz e tenta de novo
        }

        x = xNovo;
        y = yNovo;
        z = zNovo;
        gFinal = gNovo;
    }
}

int main()
{
    std::mt19937 gerador(std::random_device{}());
    std::uniform_real_distribution<double> sorteio(CHUTE_MIN, CHUTE_MAX);

    printf("==========================================================\n");
    printf(" Sistema nao linear - Metodo do Gradiente (steepest descent)\n");
    printf(" Testando chute fixo (0.5,0.5,0.5) + %d chutes aleatorios em [%.1f, %.1f]\n", NUM_CHUTES, CHUTE_MIN, CHUTE_MAX);
    printf("==========================================================\n\n");

    double melhorX = 0, melhorY = 0, melhorZ = 0, melhorG = 1e300;
    double melhorX0 = 0, melhorY0 = 0, melhorZ0 = 0;
    int melhorIteracoes = 0;
    bool melhorConvergiu = false;

    // chute fixo conhecido: converge em 99 iteracoes, testado junto com os aleatorios
    const double CHUTE_FIXO_X = 0.5, CHUTE_FIXO_Y = 0.5, CHUTE_FIXO_Z = 0.5;

    const char *NOME_ARQUIVO_RESIDUOS = "residuos_iteracoes.txt";
    FILE *arquivoResiduos = fopen(NOME_ARQUIVO_RESIDUOS, "w");
    if (arquivoResiduos == nullptr)
    {
        printf("AVISO: nao foi possivel abrir '%s' para gravar os residuos.\n", NOME_ARQUIVO_RESIDUOS);
    }
    else
    {
        fprintf(arquivoResiduos, "==========================================================\n");
        fprintf(arquivoResiduos, " Residuos por iteracao - Metodo do Gradiente\n");
        fprintf(arquivoResiduos, " Chute fixo (0.5,0.5,0.5) + %d chutes aleatorios em [%.1f, %.1f]\n", NUM_CHUTES, CHUTE_MIN, CHUTE_MAX);
        fprintf(arquivoResiduos, "==========================================================\n");
    }

    auto inicioCronometro = std::chrono::high_resolution_clock::now();

    for (int tentativa = 0; tentativa <= NUM_CHUTES; tentativa++)
    {
        double x0, y0, z0;
        if (tentativa == 0)
        {
            x0 = CHUTE_FIXO_X;
            y0 = CHUTE_FIXO_Y;
            z0 = CHUTE_FIXO_Z;
        }
        else
        {
            x0 = sorteio(gerador);
            y0 = sorteio(gerador);
            z0 = sorteio(gerador);
        }

        double x, y, z, gFinal;
        int iteracoes;
        bool convergiu;
        resolverGradiente(x0, y0, z0, x, y, z, iteracoes, convergiu, gFinal, arquivoResiduos);

        printf("Tentativa %2s | chute = (% .4f, % .4f, % .4f) | %s | iter = %6d | G = %.3e\n",
               tentativa == 0 ? "fixo" : std::to_string(tentativa).c_str(),
               x0, y0, z0, convergiu ? "convergiu" : "nao convergiu", iteracoes, gFinal);

        if (gFinal < melhorG)
        {
            melhorG = gFinal;
            melhorX = x;
            melhorY = y;
            melhorZ = z;
            melhorX0 = x0;
            melhorY0 = y0;
            melhorZ0 = z0;
            melhorIteracoes = iteracoes;
            melhorConvergiu = convergiu;
        }
    }

    auto fimCronometro = std::chrono::high_resolution_clock::now();
    double tempoTotalMs = std::chrono::duration<double, std::milli>(fimCronometro - inicioCronometro).count();

    if (arquivoResiduos != nullptr)
    {
        fclose(arquivoResiduos);
        printf("Residuos de cada iteracao gravados em '%s'\n\n", NOME_ARQUIVO_RESIDUOS);

        // gera o grafico comparativo a partir do arquivo de residuos recem-criado
        int statusGrafico = system("powershell -ExecutionPolicy Bypass -File gerar_grafico_residuos.ps1");
        if (statusGrafico != 0)
            printf("AVISO: nao foi possivel gerar o grafico (gerar_grafico_residuos.ps1)\n\n");
    }

    double f1, f2, f3;
    residuos(melhorX, melhorY, melhorZ, f1, f2, f3);

    printf("\n----------------------------------------------------------\n");
    printf("MELHOR CONVERGENCIA (menor G entre o chute fixo + %d aleatorios)\n", NUM_CHUTES);
    printf("----------------------------------------------------------\n");
    printf("Chute inicial : x = %.6f, y = %.6f, z = %.6f\n", melhorX0, melhorY0, melhorZ0);
    if (melhorConvergiu)
        printf("Convergiu na iteracao %d\n", melhorIteracoes);
    else
        printf("Nao convergiu (atingiu o numero maximo de iteracoes = %d)\n", MAX_ITERACOES);

    printf("\nSolucao encontrada:\n");
    printf("  x = %.10f\n", melhorX);
    printf("  y = %.10f\n", melhorY);
    printf("  z = %.10f\n", melhorZ);

    printf("\nResiduos finais (devem ser ~ 0):\n");
    printf("  f1 = %.3e\n", f1);
    printf("  f2 = %.3e\n", f2);
    printf("  f3 = %.3e\n", f3);
    printf("  G  = %.3e\n", melhorG);
    printf("\nTempo de processamento (todas as %d tentativas): %.3f ms\n", NUM_CHUTES + 1, tempoTotalMs);
    printf("----------------------------------------------------------\n");

    return 0;
}
