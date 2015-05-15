/*
 * Estrutura costumer, armazena a coordenada de cada cliente
 * e sua demanda
 */
struct costumer{
	int x;
	int y;
	int demand;
};

/* 
 * Estrutura designa, relaciona um cliente à uma mediana
 * para agilizar o processo de designação.
 * */
struct designa{
	int custom;
	int median;
};

typedef struct designa  designa;
typedef struct costumer costumer;

void cheur(costumer *, int, double **, int, int, int **, int *,int);
void lsheur(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);
void ag(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);

void re_medians(int n, int medians,int *meds, int **x, int *y, double **dist);
void designar(costumer *,int , int , double **, int,int **, int *);
void designar_rate(costumer *costu, int medians, int n, double **dist, int *caps,int cap,int **x, int *meds);
void sort(designa *,int ,double **,costumer *,int *,int);
double calculate(costumer *, int, double **, int, int, int **, int *);
double print_short(int, int, double);
double print_report(int, int, double, costumer *, int, double **, int, int, int **, int *);
int factivel(costumer *costu,int n, int medians, int cap, int **x, int *y);

/* Local Search */
void lsheur2(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);
void lsheur3(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);
void melhor_subs(int *meds,int fi, int *L1, int *L2, int n, int medians, double **d, int *fr,double *lucro);
void atua_l(int *L1, int *L2, int n, int *meds, int medians, double **dist);
void redesignar(costumer *costu, int medians, int n, double **dist, int cap,int **x, int *meds, int fi, int fr, int *L1, int *L2);
void reatua_l(int *L1, int *L2, int n, int *meds, int medians, double **dist,int fi, int fr);
double refine(int **x, int *y, int n, int medians, double answer, double **dist,int cap, int *meds, costumer *costu);
double d_min(double a, double b);

/* tabu */
void ts(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int T, int max_it);
void update_tabu(int *tabu,int n);

/* ACO */
void aco(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);

/* GEO */
void geo(costumer *costu, int n, double **dist, int medians, int cap, int **x, int *y, int max_it);
