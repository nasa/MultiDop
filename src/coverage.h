#ifndef _COVERAGE_H_
#define _COVERAGE_H_

int Coverage_Read(const char *, size_t);
void Coverage_UseDefault(void);
int Get_Cvg_Opt_BG(int, int *);
int Get_Cvg_Sub_BG(int, int *);
int Get_Cvg_Opt_Fil(int, int *);
int Get_Cvg_Sub_Fil(int, int *);
int Get_Cvg_BG(int);
int Get_Cvg_Fil(int);
double Get_SSeq_Trip(int);

#endif
