/*
    ==========================================================================
    TRABALHO 1 - RESOLUCAO DE SISTEMA LINEAR (36 x 36)
    ==========================================================================

    Objetivo didatico:
        Este programa resolve o sistema linear  A * x = b  de tres formas:

            1) Metodo DIRETO   -> Eliminacao de Gauss com pivotamento parcial
                                   (serve de "gabarito" para comparar as
                                   respostas dos metodos iterativos)

            2) Metodo iterativo de JACOBI
            3) Metodo iterativo de GAUSS-SEIDEL

        O codigo foi escrito de forma bem PROCEDURAL (sem classes, sem
        orientacao a objetos, sem funcoes "magicas"), usando apenas
        funcoes simples e vetores/matrizes de tamanho fixo, para facilitar
        o entendimento de todos os integrantes do grupo.

    Criterio de parada dos metodos iterativos (conforme enunciado):

            |x_i^(k+1) - x_i^(k)| < 10^-4        (tolerancia "normal")
            |x_i^(k+1) - x_i^(k)| < 10^-8        (tolerancia mais rigorosa)

        Ou seja: a cada iteracao k, calculamos a diferenca entre o valor
        novo e o valor antigo de CADA incognita x_i. O metodo para quando
        a MAIOR dessas diferencas (em modulo), entre todas as incognitas,
        for menor que a tolerancia escolhida.

        Como o enunciado apresenta as DUAS tolerancias (1e-4 e 1e-8), o
        programa executa cada metodo iterativo duas vezes: uma vez usando
        TOL = 1e-4 e outra vez usando TOL = 1e-8. Assim da para comparar
        quantas iteracoes cada tolerancia exige.

    Arquivos de entrada (devem estar na mesma pasta do executavel):
        Matriz_A.csv -> 36 linhas, 36 numeros por linha, separados por ';'
        Vetor_B.csv  -> 36 numeros, um por linha; a virgula pode ser decimal
    ==========================================================================
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>

// ---------------------------------------------------------------------
// CONSTANTES GLOBAIS
// ---------------------------------------------------------------------
const int N = 36;              // tamanho do sistema (36 equacoes / 36 incognitas)
const int MAX_ITERACOES = 1000; // numero maximo de iteracoes permitido

// as duas tolerancias pedidas no criterio de parada (imagem do enunciado)
const double TOL_1 = 1e-4;
const double TOL_2 = 1e-8;
const double CHUTE_INICIAL = 0.0;

// matrizes/vetores globais (tamanho fixo, evita alocacao dinamica -> mais didatico)
double A[N][N];
double b[N];

// ---------------------------------------------------------------------
// FUNCAO: proximoValorCSV
// Le um valor do CSV regional: ';' separa colunas e ',' representa decimal.
// Retorna true se conseguiu ler um valor, false se o arquivo acabou.
// ---------------------------------------------------------------------
bool proximoValorCSV(FILE *arquivo, double *valor)
{
    char texto[256];
    int quantidade = 0;
    int c;

    do
    {
        quantidade = 0;

        while ((c = fgetc(arquivo)) != EOF && c != ';' && c != '\n' && c != '\r')
        {
            if (quantidade < 255)
                texto[quantidade++] = (char)c;
        }

        texto[quantidade] = '\0';

        // Converte a virgula decimal para o formato aceito por strtod.
        for (int i = 0; i < quantidade; i++)
        {
            if (texto[i] == ',')
                texto[i] = '.';
        }
    } while (quantidade == 0 && c != EOF);

    if (quantidade == 0)
        return false;

    char *fim;
    *valor = strtod(texto, &fim);
    return fim != texto;
}

// ---------------------------------------------------------------------
// FUNCAO: lerMatrizA
// Le a matriz A (N x N) do arquivo CSV "Matriz_A.csv" (36 valores por
// linha, separados por virgula ou ponto e virgula)
// ---------------------------------------------------------------------
bool lerMatrizA(const char *nomeArquivo)
{
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (arquivo == NULL)
    {
        printf("ERRO: nao foi possivel abrir o arquivo %s\n", nomeArquivo);
        return false;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (!proximoValorCSV(arquivo, &A[i][j]))
            {
                printf("ERRO: arquivo %s tem menos de %d x %d valores\n",
                       nomeArquivo, N, N);
                fclose(arquivo);
                return false;
            }
        }
    }

    fclose(arquivo);
    return true;
}

// ---------------------------------------------------------------------
// FUNCAO: lerVetorB
// Le o vetor b (N valores) do arquivo CSV "Vetor_B.csv" (separados por
// virgula e/ou quebra de linha)
// ---------------------------------------------------------------------
bool lerVetorB(const char *nomeArquivo)
{
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (arquivo == NULL)
    {
        printf("ERRO: nao foi possivel abrir o arquivo %s\n", nomeArquivo);
        return false;
    }

    for (int i = 0; i < N; i++)
    {
        if (!proximoValorCSV(arquivo, &b[i]))
        {
            printf("ERRO: arquivo %s tem menos de %d valores\n", nomeArquivo, N);
            fclose(arquivo);
            return false;
        }
    }

    fclose(arquivo);
    return true;
}

// ---------------------------------------------------------------------
// FUNCAO: imprimirVetor
// Imprime um vetor de tamanho N na tela (uso geral, para depuracao)
// ---------------------------------------------------------------------
void imprimirVetor(const char *titulo, double v[N])
{
    printf("%s\n", titulo);
    for (int i = 0; i < N; i++)
    {
        printf("  x[%2d] = %12.6f\n", i, v[i]);
    }
    printf("\n");
}

void registrarResultado(FILE *arquivo, const char *metodo, double tempo,
                        int iteracoes, double solucao[N], double erro)
{
    fprintf(arquivo, "--- %s ---\n", metodo);
    fprintf(arquivo, "Tempo exigido: %.6f ms\n", tempo);
    if (iteracoes >= 0)
    {
        fprintf(arquivo, "Status: CONVERGIU\n");
        fprintf(arquivo, "Iteracoes ate convergir: %d\n", iteracoes);
    }
    else
    {
        fprintf(arquivo, "Status: NAO CONVERGIU\n");
        fprintf(arquivo, "Iteracoes ate convergir: nao atingiu o criterio\n");
    }
    fprintf(arquivo, "Solucao:\n");
    for (int i = 0; i < N; i++)
        fprintf(arquivo, "  x[%2d] = %12.6f\n", i, solucao[i]);
    fprintf(arquivo, "Maior diferenca em relacao ao metodo direto: %e\n\n", erro);
}

// ---------------------------------------------------------------------
// FUNCAO: verificarDiagonalDominante
// Verifica (apenas de forma informativa) se a matriz A e diagonalmente
// dominante. Essa e uma condicao SUFICIENTE (mas nao obrigatoria) para
// garantir a convergencia dos metodos de Jacobi e Gauss-Seidel.
// ---------------------------------------------------------------------
void verificarDiagonalDominante()
{
    bool dominante = true;

    for (int i = 0; i < N; i++)
    {
        double somaForaDiagonal = 0.0;

        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                somaForaDiagonal += fabs(A[i][j]);
            }
        }

        if (fabs(A[i][i]) < somaForaDiagonal)
        {
            dominante = false;
            break;
        }
    }

    if (dominante)
        printf("Observacao: a matriz A E diagonalmente dominante (convergencia esperada).\n\n");
    else
        printf("Observacao: a matriz A NAO e diagonalmente dominante (convergencia nao garantida).\n\n");
}

// ---------------------------------------------------------------------
// FUNCAO: eliminacaoGauss
// Resolve A * x = b pelo metodo DIRETO de eliminacao de Gauss com
// pivotamento parcial. Usado apenas como "gabarito" para comparacao.
//
// IMPORTANTE: como A e b sao usados tambem nos metodos iterativos,
// aqui trabalhamos em COPIAS locais para nao alterar os dados originais.
// ---------------------------------------------------------------------
void eliminacaoGauss(double x[N])
{
    // copias locais da matriz e do vetor (matriz aumentada [A | b])
    double M[N][N];
    double v[N];

    for (int i = 0; i < N; i++)
    {
        v[i] = b[i];
        for (int j = 0; j < N; j++)
            M[i][j] = A[i][j];
    }

    // ---- Etapa 1: escalonamento (transformar em triangular superior) ----
    for (int k = 0; k < N - 1; k++)
    {
        // pivotamento parcial: escolhe a linha com maior valor absoluto na coluna k
        int linhaPivo = k;
        double maiorValor = fabs(M[k][k]);

        for (int i = k + 1; i < N; i++)
        {
            if (fabs(M[i][k]) > maiorValor)
            {
                maiorValor = fabs(M[i][k]);
                linhaPivo = i;
            }
        }

        // troca as linhas k e linhaPivo, se necessario
        if (linhaPivo != k)
        {
            for (int j = 0; j < N; j++)
            {
                double temp = M[k][j];
                M[k][j] = M[linhaPivo][j];
                M[linhaPivo][j] = temp;
            }
            double temp = v[k];
            v[k] = v[linhaPivo];
            v[linhaPivo] = temp;
        }

        // zera os elementos abaixo do pivo na coluna k
        for (int i = k + 1; i < N; i++)
        {
            double fator = M[i][k] / M[k][k];

            for (int j = k; j < N; j++)
                M[i][j] = M[i][j] - fator * M[k][j];

            v[i] = v[i] - fator * v[k];
        }
    }

    // ---- Etapa 2: substituicao regressiva (de baixo para cima) ----
    for (int i = N - 1; i >= 0; i--)
    {
        double soma = v[i];

        for (int j = i + 1; j < N; j++)
            soma -= M[i][j] * x[j];

        x[i] = soma / M[i][i];
    }
}

// ---------------------------------------------------------------------
// FUNCAO: metodoJacobi
// Resolve A * x = b pelo metodo iterativo de JACOBI.
//
// Ideia do metodo: a partir de um chute inicial x^(0), calcula-se
// x^(k+1) usando SOMENTE os valores de x^(k) (iteracao "antiga").
//
// Criterio de parada: paramos quando, para TODO i,
//         |x_i^(k+1) - x_i^(k)| < tolerancia
// (equivalente a dizer que o MAIOR desses valores e menor que a tolerancia)
// ---------------------------------------------------------------------
void metodoJacobi(double x[N], double tolerancia, int *iteracoesRealizadas)
{
    double xAntigo[N];
    double xNovo[N];

    // chute inicial: vetor nulo
    for (int i = 0; i < N; i++)
        xAntigo[i] = CHUTE_INICIAL;

    int k;
    for (k = 0; k < MAX_ITERACOES; k++)
    {
        // calcula x_i^(k+1) para cada linha i, usando somente xAntigo
        for (int i = 0; i < N; i++)
        {
            double soma = b[i];

            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    soma -= A[i][j] * xAntigo[j];
            }

            xNovo[i] = soma / A[i][i];
        }

        // ---- criterio de parada: maior diferenca |x_i^(k+1) - x_i^(k)| ----
        double maiorDiferenca = 0.0;
        for (int i = 0; i < N; i++)
        {
            double diferenca = fabs(xNovo[i] - xAntigo[i]);
            if (diferenca > maiorDiferenca)
                maiorDiferenca = diferenca;
        }

        // atualiza xAntigo para a proxima iteracao
        for (int i = 0; i < N; i++)
            xAntigo[i] = xNovo[i];

        if (maiorDiferenca < tolerancia)
        {
            k++; // conta a iteracao atual antes de sair
            break;
        }
    }

    *iteracoesRealizadas = k;

    for (int i = 0; i < N; i++)
        x[i] = xAntigo[i];
}

// ---------------------------------------------------------------------
// FUNCAO: metodoGaussSeidel
// Resolve A * x = b pelo metodo iterativo de GAUSS-SEIDEL.
//
// Diferenca em relacao ao Jacobi: aqui usamos os valores JA ATUALIZADOS
// de x dentro da mesma iteracao (nao esperamos terminar a iteracao toda
// para usar o valor novo). Isso normalmente faz o metodo convergir mais rapido.
//
// Criterio de parada: igual ao de Jacobi.
// ---------------------------------------------------------------------
void metodoGaussSeidel(double x[N], double tolerancia, int *iteracoesRealizadas)
{
    double xAtual[N];

    // chute inicial: vetor nulo
    for (int i = 0; i < N; i++)
        xAtual[i] = CHUTE_INICIAL;

    int k;
    for (k = 0; k < MAX_ITERACOES; k++)
    {
        double maiorDiferenca = 0.0;

        for (int i = 0; i < N; i++)
        {
            double soma = b[i];

            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    soma -= A[i][j] * xAtual[j]; // usa valores ja atualizados nesta iteracao
            }

            double valorNovo = soma / A[i][i];

            double diferenca = fabs(valorNovo - xAtual[i]);
            if (diferenca > maiorDiferenca)
                maiorDiferenca = diferenca;

            xAtual[i] = valorNovo; // atualiza imediatamente
        }

        if (maiorDiferenca < tolerancia)
        {
            k++;
            break;
        }
    }

    *iteracoesRealizadas = k;

    for (int i = 0; i < N; i++)
        x[i] = xAtual[i];
}

// ---------------------------------------------------------------------
// FUNCOES: metodoJacobiOrdem2 e metodoGaussSeidelOrdem2
// Versoes de ordem 2 que usam as duas aproximacoes anteriores para
// extrapolar a nova aproximacao.
// ---------------------------------------------------------------------
void metodoJacobiOrdem2(double x[N], double tolerancia, int *iteracoesRealizadas)
{
    double xAnterior[N];
    double xAtual[N];
    double xJacobi[N];

    for (int i = 0; i < N; i++)
    {
        xAnterior[i] = CHUTE_INICIAL;
        xAtual[i] = CHUTE_INICIAL;
    }

    int k;
    for (k = 0; k < MAX_ITERACOES; k++)
    {
        for (int i = 0; i < N; i++)
        {
            double soma = b[i];
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    soma -= A[i][j] * xAtual[j];
            }
            xJacobi[i] = soma / A[i][i];
        }

        double maiorDiferenca = 0.0;
        for (int i = 0; i < N; i++)
        {
            double valorNovo = xJacobi[i];
            if (k > 0)
                valorNovo += 0.5 * (xJacobi[i] - xAnterior[i]);

            double diferenca = fabs(valorNovo - xAtual[i]);
            if (diferenca > maiorDiferenca)
                maiorDiferenca = diferenca;

            xAnterior[i] = xAtual[i];
            xAtual[i] = valorNovo;
        }

        if (maiorDiferenca < tolerancia)
        {
            k++;
            break;
        }
    }

    *iteracoesRealizadas = k;
    for (int i = 0; i < N; i++)
        x[i] = xAtual[i];
}

void metodoGaussSeidelOrdem2(double x[N], double tolerancia,
                             int *iteracoesRealizadas)
{
    double xAnterior[N];
    double xAtual[N];
    double xSeidel[N];

    for (int i = 0; i < N; i++)
    {
        xAnterior[i] = CHUTE_INICIAL;
        xAtual[i] = CHUTE_INICIAL;
    }

    int k;
    for (k = 0; k < MAX_ITERACOES; k++)
    {
        for (int i = 0; i < N; i++)
            xSeidel[i] = xAtual[i];

        for (int i = 0; i < N; i++)
        {
            double soma = b[i];
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    soma -= A[i][j] * xSeidel[j];
            }
            xSeidel[i] = soma / A[i][i];
        }

        double maiorDiferenca = 0.0;
        for (int i = 0; i < N; i++)
        {
            double valorNovo = xSeidel[i];
            if (k > 0)
                valorNovo += 0.5 * (xSeidel[i] - xAnterior[i]);

            double diferenca = fabs(valorNovo - xAtual[i]);
            if (diferenca > maiorDiferenca)
                maiorDiferenca = diferenca;

            xAnterior[i] = xAtual[i];
            xAtual[i] = valorNovo;
        }

        if (maiorDiferenca < tolerancia)
        {
            k++;
            break;
        }
    }

    *iteracoesRealizadas = k;
    for (int i = 0; i < N; i++)
        x[i] = xAtual[i];
}

// ---------------------------------------------------------------------
// FUNCAO: calcularErroMaximo
// Calcula a maior diferenca absoluta entre dois vetores (usada para
// comparar a solucao iterativa com a solucao "gabarito" da Eliminacao de Gauss)
// ---------------------------------------------------------------------
double calcularErroMaximo(double v1[N], double v2[N])
{
    double maiorErro = 0.0;
    for (int i = 0; i < N; i++)
    {
        double erro = fabs(v1[i] - v2[i]);
        if (erro > maiorErro)
            maiorErro = erro;
    }
    return maiorErro;
}

// ---------------------------------------------------------------------
// FUNCOES auxiliares usadas pelo Gradiente Conjugado e pelo CGS.
// ---------------------------------------------------------------------
void produtoMatrizVetor(double vetor[N], double resultado[N])
{
    for (int i = 0; i < N; i++)
    {
        resultado[i] = 0.0;
        for (int j = 0; j < N; j++)
            resultado[i] += A[i][j] * vetor[j];
    }
}

double produtoEscalar(double v1[N], double v2[N])
{
    double resultado = 0.0;
    for (int i = 0; i < N; i++)
        resultado += v1[i] * v2[i];
    return resultado;
}

double normaMaxima(double vetor[N])
{
    double maior = 0.0;
    for (int i = 0; i < N; i++)
    {
        if (fabs(vetor[i]) > maior)
            maior = fabs(vetor[i]);
    }
    return maior;
}

bool matrizSimetrica()
{
    const double toleranciaSimetria = 1e-12;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if (fabs(A[i][j] - A[j][i]) > toleranciaSimetria)
                return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// FUNCAO: metodoGradienteConjugado
// Resolve A * x = b para matrizes simetricas (preferencialmente definidas
// positivas), usando o residuo como criterio de parada.
// ---------------------------------------------------------------------
void metodoGradienteConjugado(double x[N], double tolerancia, int *iteracoesRealizadas)
{
    double residuo[N];
    double direcao[N];
    double produtoDirecao[N];

    for (int i = 0; i < N; i++)
    {
        x[i] = CHUTE_INICIAL;
        residuo[i] = b[i];
        direcao[i] = residuo[i];
    }

    double residuoAnterior = produtoEscalar(residuo, residuo);
    int k;

    for (k = 0; k < MAX_ITERACOES; k++)
    {
        produtoMatrizVetor(direcao, produtoDirecao);
        double denominador = produtoEscalar(direcao, produtoDirecao);

        if (fabs(denominador) < 1e-30)
            break;

        double alfa = residuoAnterior / denominador;
        for (int i = 0; i < N; i++)
        {
            x[i] += alfa * direcao[i];
            residuo[i] -= alfa * produtoDirecao[i];
        }

        if (normaMaxima(residuo) < tolerancia)
        {
            k++;
            break;
        }

        double residuoAtual = produtoEscalar(residuo, residuo);
        double beta = residuoAtual / residuoAnterior;
        for (int i = 0; i < N; i++)
            direcao[i] = residuo[i] + beta * direcao[i];

        residuoAnterior = residuoAtual;
    }

    *iteracoesRealizadas = k;
}

// ---------------------------------------------------------------------
// FUNCAO: metodoGradienteConjugadoQuadrado
// Resolve A * x = b para matrizes nao simetricas pelo metodo CGS.
// ---------------------------------------------------------------------
void metodoGradienteConjugadoQuadrado(double x[N], double tolerancia,
                                      int *iteracoesRealizadas)
{
    double residuo[N];
    double residuoBase[N];
    double p[N];
    double q[N];
    double u[N];
    double soma[N];
    double produtoSoma[N];

    for (int i = 0; i < N; i++)
    {
        x[i] = CHUTE_INICIAL;
        residuo[i] = b[i];
        residuoBase[i] = residuo[i];
        p[i] = 0.0;
        q[i] = 0.0;
    }

    double rhoAnterior = 1.0;
    bool convergiu = false;
    int k;

    for (k = 0; k < MAX_ITERACOES; k++)
    {
        double rhoAtual = produtoEscalar(residuoBase, residuo);
        if (fabs(rhoAtual) < 1e-30)
            break;

        double beta = (k == 0) ? 0.0 : rhoAtual / rhoAnterior;
        for (int i = 0; i < N; i++)
        {
            u[i] = residuo[i] + beta * q[i];
            p[i] = u[i] + beta * (q[i] + beta * p[i]);
        }

        double produtoP[N];
        produtoMatrizVetor(p, produtoP);
        double denominador = produtoEscalar(residuoBase, produtoP);
        if (fabs(denominador) < 1e-30)
            break;

        double alfa = rhoAtual / denominador;
        for (int i = 0; i < N; i++)
            q[i] = u[i] - alfa * produtoP[i];

        for (int i = 0; i < N; i++)
            soma[i] = u[i] + q[i];

        produtoMatrizVetor(soma, produtoSoma);
        for (int i = 0; i < N; i++)
        {
            x[i] += alfa * soma[i];
            residuo[i] -= alfa * produtoSoma[i];
        }

        if (normaMaxima(residuo) < tolerancia)
        {
            convergiu = true;
            k++;
            break;
        }

        rhoAnterior = rhoAtual;
    }

    *iteracoesRealizadas = convergiu ? k : -1;
}

// ---------------------------------------------------------------------
// FUNCAO PRINCIPAL
// ---------------------------------------------------------------------
int main()
{
    FILE *arquivoResultados = fopen("resultados.txt", "w");
    if (arquivoResultados == NULL)
    {
        printf("ERRO: nao foi possivel criar o arquivo resultados.txt\n");
        return 1;
    }

    fprintf(arquivoResultados, "RESULTADOS DA RESOLUCAO DO SISTEMA LINEAR\n\n");

    printf("==========================================================\n");
    printf(" RESOLUCAO DE SISTEMA LINEAR 36 x 36 (A * x = b)\n");
    printf("==========================================================\n\n");

    // -------- leitura dos arquivos de entrada --------
    if (!lerMatrizA("Matriz_A.csv"))
    {
        fclose(arquivoResultados);
        return 1;
    }

    if (!lerVetorB("Vetor_B.csv"))
    {
        fclose(arquivoResultados);
        return 1;
    }

    printf("Arquivos lidos com sucesso (A: %dx%d, b: %d).\n\n", N, N, N);

    verificarDiagonalDominante();

    bool ehSimetrica = matrizSimetrica();
    printf("A matriz A %s simetrica.\n\n",
           ehSimetrica ? "e" : "nao e");

    // -------- 1) solucao pelo metodo DIRETO (gabarito) --------
    double xGauss[N];
    std::chrono::high_resolution_clock::time_point inicioGauss =
        std::chrono::high_resolution_clock::now();
    eliminacaoGauss(xGauss);
    std::chrono::high_resolution_clock::time_point fimGauss =
        std::chrono::high_resolution_clock::now();
    double tempoGauss = std::chrono::duration<double, std::milli>(fimGauss - inicioGauss).count();

    printf("---------------------------------------------------------\n");
    printf("1) SOLUCAO PELA ELIMINACAO DE GAUSS (metodo direto)\n");
    printf("---------------------------------------------------------\n");
    printf("Tempo exigido: %.6f ms\n", tempoGauss);
    imprimirVetor("Solucao x:", xGauss);
    registrarResultado(arquivoResultados, "ELIMINACAO DE GAUSS", tempoGauss, 0,
                       xGauss, 0.0);

    // -------- 2) e 3) metodos iterativos, para as duas tolerancias --------
    double tolerancias[2] = {TOL_1, TOL_2};

    for (int t = 0; t < 2; t++)
    {
        double tol = tolerancias[t];

        printf("===========================================================\n");
        printf(" CRITERIO DE PARADA: |x_i^(k+1) - x_i^(k)| < %.0e\n", tol);
        printf("===========================================================\n\n");
        fprintf(arquivoResultados,
            "===========================================================\n");
        fprintf(arquivoResultados,
            "CRITERIO DE PARADA: |x_i^(k+1) - x_i^(k)| < %.0e\n",
            tol);
        fprintf(arquivoResultados,
            "===========================================================\n\n");

        // ---- Jacobi ----
        double xJacobi[N];
        int iterJacobi = 0;
        std::chrono::high_resolution_clock::time_point inicioJacobi =
            std::chrono::high_resolution_clock::now();
        metodoJacobi(xJacobi, tol, &iterJacobi);
        std::chrono::high_resolution_clock::time_point fimJacobi =
            std::chrono::high_resolution_clock::now();
        double tempoJacobi = std::chrono::duration<double, std::milli>(fimJacobi - inicioJacobi).count();

        printf("--- Metodo de JACOBI ---\n");
        printf("Tempo exigido: %.6f ms\n", tempoJacobi);
        printf("Iteracoes ate convergir: %d\n", iterJacobi);
        imprimirVetor("Solucao x (Jacobi):", xJacobi);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xJacobi, xGauss));
        registrarResultado(arquivoResultados, "JACOBI", tempoJacobi,
                   iterJacobi, xJacobi,
                   calcularErroMaximo(xJacobi, xGauss));

        // ---- Gauss-Seidel ----
        double xSeidel[N];
        int iterSeidel = 0;
        std::chrono::high_resolution_clock::time_point inicioSeidel =
            std::chrono::high_resolution_clock::now();
        metodoGaussSeidel(xSeidel, tol, &iterSeidel);
        std::chrono::high_resolution_clock::time_point fimSeidel =
            std::chrono::high_resolution_clock::now();
        double tempoSeidel = std::chrono::duration<double, std::milli>(fimSeidel - inicioSeidel).count();

        printf("--- Metodo de GAUSS-SEIDEL ---\n");
        printf("Tempo exigido: %.6f ms\n", tempoSeidel);
        printf("Iteracoes ate convergir: %d\n", iterSeidel);
        imprimirVetor("Solucao x (Gauss-Seidel):", xSeidel);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xSeidel, xGauss));
        registrarResultado(arquivoResultados, "GAUSS-SEIDEL", tempoSeidel,
                   iterSeidel, xSeidel,
                   calcularErroMaximo(xSeidel, xGauss));

        // ---- Jacobi de ordem 2 ----
        double xJacobiOrdem2[N];
        int iterJacobiOrdem2 = 0;
        std::chrono::high_resolution_clock::time_point inicioJacobiOrdem2 =
            std::chrono::high_resolution_clock::now();
        metodoJacobiOrdem2(xJacobiOrdem2, tol, &iterJacobiOrdem2);
        std::chrono::high_resolution_clock::time_point fimJacobiOrdem2 =
            std::chrono::high_resolution_clock::now();
        double tempoJacobiOrdem2 =
            std::chrono::duration<double, std::milli>(fimJacobiOrdem2 - inicioJacobiOrdem2).count();

        printf("--- Metodo de JACOBI DE ORDEM 2 ---\n");
        printf("Tempo exigido: %.6f ms\n", tempoJacobiOrdem2);
        printf("Iteracoes ate convergir: %d\n", iterJacobiOrdem2);
        imprimirVetor("Solucao x (Jacobi de ordem 2):", xJacobiOrdem2);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xJacobiOrdem2, xGauss));
        registrarResultado(arquivoResultados, "JACOBI DE ORDEM 2",
                   tempoJacobiOrdem2, iterJacobiOrdem2,
                   xJacobiOrdem2,
                   calcularErroMaximo(xJacobiOrdem2, xGauss));

        // ---- Gauss-Seidel de ordem 2 ----
        double xSeidelOrdem2[N];
        int iterSeidelOrdem2 = 0;
        std::chrono::high_resolution_clock::time_point inicioSeidelOrdem2 =
            std::chrono::high_resolution_clock::now();
        metodoGaussSeidelOrdem2(xSeidelOrdem2, tol, &iterSeidelOrdem2);
        std::chrono::high_resolution_clock::time_point fimSeidelOrdem2 =
            std::chrono::high_resolution_clock::now();
        double tempoSeidelOrdem2 =
            std::chrono::duration<double, std::milli>(fimSeidelOrdem2 - inicioSeidelOrdem2).count();

        printf("--- Metodo de GAUSS-SEIDEL DE ORDEM 2 ---\n");
        printf("Tempo exigido: %.6f ms\n", tempoSeidelOrdem2);
        printf("Iteracoes ate convergir: %d\n", iterSeidelOrdem2);
        imprimirVetor("Solucao x (Gauss-Seidel de ordem 2):", xSeidelOrdem2);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xSeidelOrdem2, xGauss));
        registrarResultado(arquivoResultados, "GAUSS-SEIDEL DE ORDEM 2",
                   tempoSeidelOrdem2, iterSeidelOrdem2,
                   xSeidelOrdem2,
                   calcularErroMaximo(xSeidelOrdem2, xGauss));

        if (ehSimetrica)
        {
            double xGradiente[N];
            int iterGradiente = 0;
            std::chrono::high_resolution_clock::time_point inicioGradiente =
                std::chrono::high_resolution_clock::now();
            metodoGradienteConjugado(xGradiente, tol, &iterGradiente);
            std::chrono::high_resolution_clock::time_point fimGradiente =
                std::chrono::high_resolution_clock::now();
            double tempoGradiente = std::chrono::duration<double, std::milli>(fimGradiente - inicioGradiente).count();

            printf("--- Metodo do GRADIENTE CONJUGADO ---\n");
            printf("Tempo exigido: %.6f ms\n", tempoGradiente);
            printf("Iteracoes ate convergir: %d\n", iterGradiente);
            imprimirVetor("Solucao x (Gradiente Conjugado):", xGradiente);
            printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
                   calcularErroMaximo(xGradiente, xGauss));
                 registrarResultado(arquivoResultados, "GRADIENTE CONJUGADO",
                           tempoGradiente, iterGradiente, xGradiente,
                           calcularErroMaximo(xGradiente, xGauss));
        }
        else
        {
            double xCGS[N];
            int iterCGS = 0;
            std::chrono::high_resolution_clock::time_point inicioCGS =
                std::chrono::high_resolution_clock::now();
            metodoGradienteConjugadoQuadrado(xCGS, tol, &iterCGS);
            std::chrono::high_resolution_clock::time_point fimCGS =
                std::chrono::high_resolution_clock::now();
            double tempoCGS = std::chrono::duration<double, std::milli>(fimCGS - inicioCGS).count();

            printf("--- Metodo do GRADIENTE CONJUGADO QUADRADO (CGS) ---\n");
            printf("Tempo exigido: %.6f ms\n", tempoCGS);
            if (iterCGS >= 0)
                printf("Iteracoes ate convergir: %d\n", iterCGS);
            else
                printf("Status: NAO CONVERGIU (criterio nao atingido)\n");
            imprimirVetor("Solucao x (CGS):", xCGS);
            printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
                   calcularErroMaximo(xCGS, xGauss));
            registrarResultado(arquivoResultados,
                               "GRADIENTE CONJUGADO QUADRADO (CGS)",
                               tempoCGS, iterCGS, xCGS,
                               calcularErroMaximo(xCGS, xGauss));
        }
    }

    fclose(arquivoResultados);
    return 0;
}
