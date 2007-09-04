MODEL * NewCodonModel ( const int gencode, const double kappa, const double omega, const double *pi, const int codonf, const int freq_type);
MODEL *NewCodonModel_singleDnDs (const int gencode, const double kappa, const double omega, const double *pi, const int codonf, const int freq_type);
MODEL * NewCodonModel_single ( const int gencode, const double kappa, const double omega, const double *pi, const int codonf, const int freq_type);
MODEL * NewCodonModel_full ( const int gencode, const double kappa, const double omega, const double *pi, const int codonf, const int freq_type);
double GetScale_single ( MODEL * model, const double f);
void SetAminoAndCodonFuncs ( const int nucleo_type, const int amino_type, const char * nucleofile, const char * aminofile);
double *  GetEquilibriumDistCodon (const double * pi,const int codonf, const int gencode);
