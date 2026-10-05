/*
    ========================================================================
    TRABALHO 1 - SISTEMA NAO LINEAR (METODO DE NEWTON)
    ========================================================================

    Sistema:
        f1(x,y,z) = (x-1)^2 + (y-1)^2 + (z-1)^2 - 1 = 0
        f2(x,y,z) = 2*x^2 + (y-1)^2 - 4*z = 0
        f3(x,y,z) = 3*x^2 + 2*z^2 - 4*y = 0

    Em cada iteracao, o metodo de Newton resolve J * delta = -F e atualiza
    (x,y,z) com (x,y,z) + alpha * delta. O alpha usa backtracking para
    reduzir a soma dos quadrados dos residuos e aumentar a robustez do metodo.
    ========================================================================
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <string>

const int MAX_ITERACOES = 1000;
const int MAX_BACKTRACKING = 60;
const int NUM_CHUTES = 10;
const double TOL_RESIDUO = 1e-10;
const double CHUTE_MIN = -1.0;
const double CHUTE_MAX = 1.0;

void residuos(double x, double y, double z, double &f1, double &f2, double &f3)
{
    f1 = (x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1) - 1;
    f2 = 2 * x * x + (y - 1) * (y - 1) - 4 * z;
    f3 = 3 * x * x + 2 * z * z - 4 * y;
}

double funcaoObjetivo(double x, double y, double z)
{
    double f1, f2, f3;
    residuos(x, y, z, f1, f2, f3);
    return f1 * f1 + f2 * f2 + f3 * f3;
}

bool resolverSistemaLinear(double matriz[3][3], double vetor[3], double solucao[3])
{
    double aumentado[3][4];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            aumentado[i][j] = matriz[i][j];
        aumentado[i][3] = vetor[i];
    }

    for (int coluna = 0; coluna < 3; coluna++)
    {
        int pivo = coluna;
        for (int linha = coluna + 1; linha < 3; linha++)
            if (std::fabs(aumentado[linha][coluna]) > std::fabs(aumentado[pivo][coluna]))
                pivo = linha;

        if (std::fabs(aumentado[pivo][coluna]) < 1e-14)
            return false;

        for (int j = coluna; j < 4; j++)
        {
            double temporario = aumentado[coluna][j];
            aumentado[coluna][j] = aumentado[pivo][j];
            aumentado[pivo][j] = temporario;
        }

        for (int linha = coluna + 1; linha < 3; linha++)
        {
            double fator = aumentado[linha][coluna] / aumentado[coluna][coluna];
            for (int j = coluna; j < 4; j++)
                aumentado[linha][j] -= fator * aumentado[coluna][j];
        }
    }

    for (int linha = 2; linha >= 0; linha--)
    {
        solucao[linha] = aumentado[linha][3];
        for (int j = linha + 1; j < 3; j++)
            solucao[linha] -= aumentado[linha][j] * solucao[j];
        solucao[linha] /= aumentado[linha][linha];
    }
    return true;
}

bool resolverNewton(double x0, double y0, double z0,
                    double &x, double &y, double &z, int &iteracoes, double &gFinal,
                    FILE *arquivo)
{
    x = x0;
    y = y0;
    z = z0;

    if (arquivo != nullptr)
    {
        fprintf(arquivo, "\n----------------------------------------------------------\n");
        fprintf(arquivo, "Chute inicial: x0 = %.6f, y0 = %.6f, z0 = %.6f\n", x0, y0, z0);
        fprintf(arquivo, "----------------------------------------------------------\n");
        fprintf(arquivo, "%8s %16s %16s %16s %16s\n", "iter", "f1", "f2", "f3", "G");
    }

    for (iteracoes = 0; iteracoes < MAX_ITERACOES; iteracoes++)
    {
        double f1, f2, f3;
        residuos(x, y, z, f1, f2, f3);
        gFinal = f1 * f1 + f2 * f2 + f3 * f3;

        if (arquivo != nullptr)
            fprintf(arquivo, "%8d %16.6e %16.6e %16.6e %16.6e\n",
                    iteracoes, f1, f2, f3, gFinal);

        if (std::sqrt(gFinal) < TOL_RESIDUO)
            return true;

        double jacobiana[3][3] = {
            {2 * (x - 1), 2 * (y - 1), 2 * (z - 1)},
            {4 * x, 2 * (y - 1), -4},
            {6 * x, -4, 4 * z}};
        double vetor[3] = {-f1, -f2, -f3};
        double delta[3];
        if (!resolverSistemaLinear(jacobiana, vetor, delta))
            return false;

        double alpha = 1.0;
        bool aceitouPasso = false;
        for (int tentativa = 0; tentativa < MAX_BACKTRACKING; tentativa++)
        {
            double xNovo = x + alpha * delta[0];
            double yNovo = y + alpha * delta[1];
            double zNovo = z + alpha * delta[2];
            double gNovo = funcaoObjetivo(xNovo, yNovo, zNovo);

            if (std::isfinite(gNovo) && gNovo < gFinal)
            {
                x = xNovo;
                y = yNovo;
                z = zNovo;
                gFinal = gNovo;
                aceitouPasso = true;
                break;
            }
            alpha *= 0.5;
        }

        if (!aceitouPasso)
            return false;
    }

    iteracoes = MAX_ITERACOES;
    double f1, f2, f3;
    residuos(x, y, z, f1, f2, f3);
    gFinal = f1 * f1 + f2 * f2 + f3 * f3;
    if (arquivo != nullptr)
        fprintf(arquivo, "%8d %16.6e %16.6e %16.6e %16.6e\n",
                iteracoes, f1, f2, f3, gFinal);
    return std::sqrt(gFinal) < TOL_RESIDUO;
}

int main()
{
    std::mt19937 gerador(std::random_device{}());
    std::uniform_real_distribution<double> sorteio(CHUTE_MIN, CHUTE_MAX);
    const double CHUTE_FIXO_X = 0.5, CHUTE_FIXO_Y = 0.5, CHUTE_FIXO_Z = 0.5;

    printf("==========================================================\n");
    printf(" Sistema nao linear - Metodo de Newton amortecido\n");
    printf(" Testando chute fixo (0.5,0.5,0.5) + %d chutes aleatorios em [%.1f, %.1f]\n",
           NUM_CHUTES, CHUTE_MIN, CHUTE_MAX);
    printf("==========================================================\n\n");

    const char *NOME_ARQUIVO_RESIDUOS = "residuos_iteracoes.txt";
    FILE *arquivoResiduos = std::fopen(NOME_ARQUIVO_RESIDUOS, "w");
    if (arquivoResiduos == nullptr)
        printf("AVISO: nao foi possivel abrir '%s' para gravar os residuos.\n", NOME_ARQUIVO_RESIDUOS);
    else
    {
        fprintf(arquivoResiduos, "==========================================================\n");
        fprintf(arquivoResiduos, " Residuos por iteracao - Metodo de Newton\n");
        fprintf(arquivoResiduos, " Chute fixo (0.5,0.5,0.5) + %d chutes aleatorios em [%.1f, %.1f]\n",
                NUM_CHUTES, CHUTE_MIN, CHUTE_MAX);
        fprintf(arquivoResiduos, "==========================================================\n");
    }

    double melhorX = 0, melhorY = 0, melhorZ = 0, melhorG = 1e300;
    double melhorX0 = 0, melhorY0 = 0, melhorZ0 = 0;
    int melhorIteracoes = 0;
    bool melhorConvergiu = false;

    for (int tentativa = 0; tentativa <= NUM_CHUTES; tentativa++)
    {
        double x0 = tentativa == 0 ? CHUTE_FIXO_X : sorteio(gerador);
        double y0 = tentativa == 0 ? CHUTE_FIXO_Y : sorteio(gerador);
        double z0 = tentativa == 0 ? CHUTE_FIXO_Z : sorteio(gerador);
        double x, y, z, gFinal;
        int iteracoes;
        bool convergiu = resolverNewton(x0, y0, z0, x, y, z, iteracoes, gFinal, arquivoResiduos);

        printf("Tentativa %2s | chute = (% .4f, % .4f, % .4f) | %s | iter = %4d | G = %.3e\n",
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

    if (arquivoResiduos != nullptr)
    {
        std::fclose(arquivoResiduos);
        printf("Residuos de cada iteracao gravados em '%s'\n", NOME_ARQUIVO_RESIDUOS);
        int statusGrafico = std::system("powershell -ExecutionPolicy Bypass -File gerar_grafico_residuos.ps1 -Metodo Newton -AbrirNoNavegador");
        if (statusGrafico != 0)
            printf("AVISO: nao foi possivel gerar ou abrir o grafico (gerar_grafico_residuos.ps1)\n");
    }

    double f1, f2, f3;
    residuos(melhorX, melhorY, melhorZ, f1, f2, f3);
    printf("\n----------------------------------------------------------\n");
    printf("MELHOR RESULTADO (menor G entre o chute fixo + %d aleatorios)\n", NUM_CHUTES);
    printf("----------------------------------------------------------\n");
    printf("Chute inicial : x = %.6f, y = %.6f, z = %.6f\n", melhorX0, melhorY0, melhorZ0);
    printf("Status         : %s na iteracao %d\n",
           melhorConvergiu ? "convergiu" : "nao convergiu", melhorIteracoes);
    printf("Solucao        : x = %.10f, y = %.10f, z = %.10f\n", melhorX, melhorY, melhorZ);
    printf("Residuos       : f1 = %.3e, f2 = %.3e, f3 = %.3e\n", f1, f2, f3);
    printf("G              : %.3e\n", melhorG);
    printf("----------------------------------------------------------\n");
    return melhorConvergiu ? 0 : 1;
}