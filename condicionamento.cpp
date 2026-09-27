#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <limits>

const int N = 36;
const int MAX_DIM = 2 * N;
const double INF = std::numeric_limits<double>::infinity();

typedef double Matriz[MAX_DIM][MAX_DIM];

double A[N][N];

double lerNumero(FILE *arquivo)
{
    char texto[256];
    int tamanho = 0;
    int c;
    do
    {
        tamanho = 0;
        while ((c = fgetc(arquivo)) != EOF && c != ';' && c != '\n' && c != '\r')
        {
            if (tamanho < 255)
                texto[tamanho++] = (char)c;
        }
        texto[tamanho] = '\0';
    } while (tamanho == 0 && c != EOF);

    for (int i = 0; i < tamanho; i++)
        if (texto[i] == ',')
            texto[i] = '.';

    return strtod(texto, NULL);
}

bool lerA(const char *nome)
{
    FILE *arquivo = fopen(nome, "r");
    if (arquivo == NULL)
        return false;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = lerNumero(arquivo);

    fclose(arquivo);
    return true;
}

void zerar(Matriz m, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m[i][j] = 0.0;
}

void identidade(Matriz m, int n)
{
    zerar(m, n);
    for (int i = 0; i < n; i++)
        m[i][i] = 1.0;
}

bool inverter(const Matriz entrada, Matriz saida, int n)
{
    Matriz trabalho;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            trabalho[i][j] = entrada[i][j];
    identidade(saida, n);

    for (int coluna = 0; coluna < n; coluna++)
    {
        int pivo = coluna;
        for (int i = coluna + 1; i < n; i++)
            if (fabs(trabalho[i][coluna]) > fabs(trabalho[pivo][coluna]))
                pivo = i;

        if (fabs(trabalho[pivo][coluna]) < 1e-14)
            return false;

        for (int j = 0; j < n; j++)
        {
            std::swap(trabalho[coluna][j], trabalho[pivo][j]);
            std::swap(saida[coluna][j], saida[pivo][j]);
        }

        double divisor = trabalho[coluna][coluna];
        for (int j = 0; j < n; j++)
        {
            trabalho[coluna][j] /= divisor;
            saida[coluna][j] /= divisor;
        }

        for (int i = 0; i < n; i++)
        {
            if (i == coluna)
                continue;
            double fator = trabalho[i][coluna];
            for (int j = 0; j < n; j++)
            {
                trabalho[i][j] -= fator * trabalho[coluna][j];
                saida[i][j] -= fator * saida[coluna][j];
            }
        }
    }
    return true;
}

void multiplicarTransposta(const Matriz m, Matriz resultado, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            resultado[i][j] = 0.0;
            for (int k = 0; k < n; k++)
                resultado[i][j] += m[k][i] * m[k][j];
        }
}

double normaEuclidiana(const Matriz m, int n)
{
    Matriz simetrica;
    multiplicarTransposta(m, simetrica, n);

    for (int iteracao = 0; iteracao < 100 * n; iteracao++)
    {
        int p = 0;
        int q = 1;
        double maior = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabs(simetrica[i][j]) > maior)
                {
                    maior = fabs(simetrica[i][j]);
                    p = i;
                    q = j;
                }
        if (maior < 1e-12)
            break;

        double angulo = 0.5 * atan2(2.0 * simetrica[p][q],
                                    simetrica[q][q] - simetrica[p][p]);
        double c = cos(angulo);
        double s = sin(angulo);
        for (int k = 0; k < n; k++)
        {
            double pik = simetrica[p][k];
            double qik = simetrica[q][k];
            simetrica[p][k] = c * pik - s * qik;
            simetrica[q][k] = s * pik + c * qik;
        }
        for (int k = 0; k < n; k++)
        {
            double kip = simetrica[k][p];
            double kiq = simetrica[k][q];
            simetrica[k][p] = c * kip - s * kiq;
            simetrica[k][q] = s * kip + c * kiq;
        }
    }

    double maiorAutovalor = 0.0;
    for (int i = 0; i < n; i++)
        maiorAutovalor = std::max(maiorAutovalor, simetrica[i][i]);
    return sqrt(std::max(0.0, maiorAutovalor));
}

double normaColunas(const Matriz m, int n)
{
    double maior = 0.0;
    for (int j = 0; j < n; j++)
    {
        double soma = 0.0;
        for (int i = 0; i < n; i++)
            soma += fabs(m[i][j]);
        maior = std::max(maior, soma);
    }
    return maior;
}

double normaLinhas(const Matriz m, int n)
{
    double maior = 0.0;
    for (int i = 0; i < n; i++)
    {
        double soma = 0.0;
        for (int j = 0; j < n; j++)
            soma += fabs(m[i][j]);
        maior = std::max(maior, soma);
    }
    return maior;
}

double normaFrobenius(const Matriz m, int n)
{
    double soma = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            soma += m[i][j] * m[i][j];
    return sqrt(soma);
}

double normaMaxima(const Matriz m, int n)
{
    double maior = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            maior = std::max(maior, fabs(m[i][j]));
    return maior;
}

typedef double (*Norma)(const Matriz, int);

const char *nomesNormas[] = {
    "Euclidiana", "Soma maxima de colunas", "Soma maxima de linhas",
    "Frobenius", "Maximo"};
Norma normas[] = {normaEuclidiana, normaColunas, normaLinhas,
                  normaFrobenius, normaMaxima};

void ponderar(const Matriz entrada, Matriz saida, const double pesos[], int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            saida[i][j] = pesos[i] * entrada[i][j] / pesos[j];
}

double condicao(const Matriz m, int n, Norma norma, const double pesos[],
                bool ponderada)
{
    Matriz inversa;
    if (!inverter(m, inversa, n))
        return INF;

    if (!ponderada)
        return norma(m, n) * norma(inversa, n);

    Matriz escalada;
    Matriz inversaEscalada;
    ponderar(m, escalada, pesos, n);
    ponderar(inversa, inversaEscalada, pesos, n);
    return norma(escalada, n) * norma(inversaEscalada, n);
}

void construirJacobi(Matriz m)
{
    zerar(m, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i != j)
                m[i][j] = -A[i][j] / A[i][i];
}

void construirGaussSeidel(Matriz m)
{
    zerar(m, N);
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
        {
            double soma = (i == j) ? 0.0 : -A[i][j];
            for (int k = 0; k < i; k++)
                soma -= A[i][k] * m[k][j];
            m[i][j] = soma / A[i][i];
        }
}

void construirOrdem2(Matriz m, const Matriz base)
{
    zerar(m, 2 * N);
    for (int i = 0; i < N; i++)
    {
        m[i][i] = -0.5;
        m[N + i][i] = 1.0;
        for (int j = 0; j < N; j++)
            m[i][j] += 1.5 * base[i][j];
    }
}

void imprimirResultados(FILE *saida, const char *nome, const Matriz m,
                        int n, const double pesos[])
{
    fprintf(saida, "\n%s (%d x %d)\n", nome, n, n);
    for (int i = 0; i < 5; i++)
    {
        double normal = normas[i](m, n);
        double ponderada = 0.0;
        Matriz escalada;
        ponderar(m, escalada, pesos, n);
        ponderada = normas[i](escalada, n);
        double cond = condicao(m, n, normas[i], pesos, false);
        double condPonderada = condicao(m, n, normas[i], pesos, true);
        fprintf(saida, "%-28s norma=%12.6e  cond=%12.6e\n",
            nomesNormas[i], normal, cond);
        fprintf(saida, "  ponderada                 norma=%12.6e  cond=%12.6e\n",
            ponderada, condPonderada);
    }
    fprintf(saida, "Versao ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.\n");
}

int main()
{
    if (!lerA("Matriz_A.csv"))
    {
        printf("ERRO: nao foi possivel ler Matriz_A.csv\n");
        return 1;
    }

    Matriz m1, m2, m3, m4;
    construirJacobi(m1);
    construirGaussSeidel(m2);
    construirOrdem2(m3, m1);
    construirOrdem2(m4, m2);

    double pesos[MAX_DIM];
    for (int i = 0; i < N; i++)
        pesos[i] = fabs(A[i][i]);
    for (int i = 0; i < N; i++)
    {
        pesos[N + i] = pesos[i];
    }

    Matriz matrizA;
    zerar(matrizA, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrizA[i][j] = A[i][j];

    FILE *saida = fopen("condicionamento.txt", "w");
    if (saida == NULL)
    {
        printf("ERRO: nao foi possivel criar condicionamento.txt\n");
        return 1;
    }

    fprintf(saida, "CONDICIONAMENTO DA MATRIZ A E DAS MATRIZES DE ITERACAO\n");
    fprintf(saida, "M1=Jacobi, M2=Gauss-Seidel, M3=Jacobi ordem 2, M4=Gauss-Seidel ordem 2.\n");
    fprintf(saida, "Condicionamento: kappa(B)=||B||*||B^-1||.\n");

    double pesosA[MAX_DIM];
    for (int i = 0; i < N; i++)
        pesosA[i] = pesos[i];
    imprimirResultados(saida, "A", matrizA, N, pesosA);
    imprimirResultados(saida, "M1 - Jacobi", m1, N, pesos);
    imprimirResultados(saida, "M2 - Gauss-Seidel", m2, N, pesos);
    imprimirResultados(saida, "M3 - Jacobi ordem 2", m3, 2 * N, pesos);
    imprimirResultados(saida, "M4 - Gauss-Seidel ordem 2", m4, 2 * N, pesos);

    fclose(saida);
    printf("Resultados gravados em condicionamento.txt\n");
    return 0;
}
