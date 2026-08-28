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

// ---------------------------------------------------------------------
// CONSTANTES GLOBAIS
// ---------------------------------------------------------------------
const int N = 36;              // tamanho do sistema (36 equacoes / 36 incognitas)
const int MAX_ITERACOES = 1000; // numero maximo de iteracoes permitido

// as duas tolerancias pedidas no criterio de parada (imagem do enunciado)
const double TOL_1 = 1e-4;
const double TOL_2 = 1e-8;

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
        xAntigo[i] = 0.0;

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
        xAtual[i] = 0.0;

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
// FUNCAO PRINCIPAL
// ---------------------------------------------------------------------
int main()
{
    printf("==========================================================\n");
    printf(" RESOLUCAO DE SISTEMA LINEAR 36 x 36 (A * x = b)\n");
    printf("==========================================================\n\n");

    // -------- leitura dos arquivos de entrada --------
    if (!lerMatrizA("Matriz_A.csv"))
        return 1;

    if (!lerVetorB("Vetor_B.csv"))
        return 1;

    printf("Arquivos lidos com sucesso (A: %dx%d, b: %d).\n\n", N, N, N);

    verificarDiagonalDominante();

    // -------- 1) solucao pelo metodo DIRETO (gabarito) --------
    double xGauss[N];
    eliminacaoGauss(xGauss);

    printf("---------------------------------------------------------\n");
    printf("1) SOLUCAO PELA ELIMINACAO DE GAUSS (metodo direto)\n");
    printf("---------------------------------------------------------\n");
    imprimirVetor("Solucao x:", xGauss);

    // -------- 2) e 3) metodos iterativos, para as duas tolerancias --------
    double tolerancias[2] = {TOL_1, TOL_2};

    for (int t = 0; t < 2; t++)
    {
        double tol = tolerancias[t];

        printf("===========================================================\n");
        printf(" CRITERIO DE PARADA: |x_i^(k+1) - x_i^(k)| < %.0e\n", tol);
        printf("===========================================================\n\n");

        // ---- Jacobi ----
        double xJacobi[N];
        int iterJacobi = 0;
        metodoJacobi(xJacobi, tol, &iterJacobi);

        printf("--- Metodo de JACOBI ---\n");
        printf("Iteracoes ate convergir: %d\n", iterJacobi);
        imprimirVetor("Solucao x (Jacobi):", xJacobi);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xJacobi, xGauss));

        // ---- Gauss-Seidel ----
        double xSeidel[N];
        int iterSeidel = 0;
        metodoGaussSeidel(xSeidel, tol, &iterSeidel);

        printf("--- Metodo de GAUSS-SEIDEL ---\n");
        printf("Iteracoes ate convergir: %d\n", iterSeidel);
        imprimirVetor("Solucao x (Gauss-Seidel):", xSeidel);
        printf("Maior diferenca em relacao ao metodo direto: %e\n\n",
               calcularErroMaximo(xSeidel, xGauss));
    }

    return 0;
}
