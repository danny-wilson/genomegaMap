/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the genomegaMap library.
 *
 *  The genomegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The genomegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the genomegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/********************************************/
/*	mutation.h 23rd February 2005			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/
#ifndef _MUTATION_H_
#define _MUTATION_H_
#include <vector>
#include <myutils.h>

using std::vector;
using namespace myutils;

class Mutation_Matrix {
public:
	int n_states;
	vector<double> state_freq;
	vector<string> state_char;
	Random *ran;

	Matrix<double> C;			/*continuous-time rate matrix*/
	Matrix<double> D;			/*discrete-time transition matrix*/
	vector<double> mutation_rate;
	vector<double> mutation_mean;
	
protected:
	/*assumes C_in is a valid rate matrix of size n_states_in*/
	void initialize(const int n_states_in, Matrix<double> *C_in);
public:
	Mutation_Matrix();
	Mutation_Matrix(const Mutation_Matrix &M);
	Mutation_Matrix& operator=(const Mutation_Matrix& M);
	void set_state_freq(vector<double> state_freq_in);
	void set_state_char(vector<string> state_char_in);
	void set_ran(Random *ran_in);
	inline double get_rate(const int state);
	int draw();
	int mutate(const int state);
	int mutate_edge(const int state, const double time);
	int mutate_edge(const int state, const double time, int &nmut);
	double expected_rate();

	// Build a transition probability matrix from the Eigenvectors and Eigenvalues of the symmetric component of a reversible rate matrix
	// (i.e. from calling diagonalize_symmetric)
	static void build_Pt(Matrix<double>& Pt, const Matrix<double>& SymEigenvec, const Vector<double>& SymEigenval, const Vector<double>& pi, const double t);
	// As above, but coalescent-average over the time t (scaled by scale), where the coalescent rate is lambda. Both (scale*t) and lambda are in the same time units
	static void build_Rt(Matrix<double>& Rt, const Matrix<double>& SymEigenvec, const Vector<double>& SymEigenval, const Vector<double>& pi, const double scale, const double lambda);
};

class Nucleotide_Mutation_Matrix : public Mutation_Matrix {
public:
	void set_defaults();
};

class JC69 : public Nucleotide_Mutation_Matrix {
public:
	JC69(const double lambda, Random *ran_in);
	int fast_mutate(const int state, const double time);
};

class F81 : public Nucleotide_Mutation_Matrix {
public:
	F81(const double lambda, const vector<double> state_freq_in, Random *ran_in);
};

class K80 : public Nucleotide_Mutation_Matrix {
public:
	K80(const double lambda, const double kappa, Random *ran_in);
	K80& update(const double lambda, const double kappa);
	int fast_mutate(const int state, const double time);
};

class HKY85 : public Nucleotide_Mutation_Matrix {
public:
	HKY85(const double lambda, const double kappa, const vector<double> state_freq_in, Random *ran_in);
};

class TN93 : public Nucleotide_Mutation_Matrix {
public:
	TN93(const double lambda, const double kappa_R, const double kappa_Y, const vector<double> state_freq_in, Random *ran_in);
};

class Codon_Mutation_Matrix : public Mutation_Matrix {
public:
	void set_defaults();
	virtual Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega) = 0;
	virtual Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi) = 0;
	virtual Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi) = 0;

};

class NY98 : public Codon_Mutation_Matrix {
protected:
	NY98();
public:
	NY98(Random *ran_in);
	NY98(const double mu, const double kappa, const double omega, Random *ran_in);
	NY98(const double mu, const double kappa, const double omega, const vector<double> &pi, Random *ran_in);
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega);
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi);
	Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi);
};

/*class gsl_eigen_symm_workspace_wrapper {
private:
	gsl_eigen_symm_workspace* _workspace;
public:
	gsl_eigen_symm_workspace_wrapper(const size_t n) {
		_workspace = gsl_eigen_symm_alloc(n);
	}
	~gsl_eigen_symm_workspace_wrapper() {
		gsl_eigen_symm_free(_workspace);
	}
	gsl_eigen_symm_workspace* workspace() {
		return _workspace;
	}
}*/

class NY98_61 : public NY98 {
public:
	void set_defaults();
	NY98_61(Random *ran_in);
	NY98_61(const double mu, const double kappa, const double omega, Random *ran_in);
	NY98_61(const double mu, const double kappa, const double omega, const vector<double> &pi, Random *ran_in);
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega);
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi);
	Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi);
	// Build a rate generator matrix
	static void build_C(Matrix<double>& C, const double mu, const double kappa, const double omega, const vector<double> &pi);
	// Build a symmetric rate generator matrix for faster diagonalization
	static void build_A(Matrix<double>& C, const double mu, const double kappa, const double omega, const vector<double> &pi);
	// Build a transition probability matrix
	static void diagonalize_symmetric(Matrix<double>& Eigenvec, Vector<double>& Eigenval, const double kappa, const double omega, const vector<double> &pi);
	// Build a transition probability matrix and return the mean rate
	static void diagonalize_symmetric(Matrix<double>& Eigenvec, double& meanrate, Vector<double>& Eigenval, const double kappa, const double omega, const vector<double> &pi);
};

class FSM_Binary : public Mutation_Matrix {
/****************************************************************/
/*	Mutations occur at rate lambda/2 per unit time.				*/
/*																*/
/*	Transition probability matrix, given time t is				*/
/*																*/
/*	P[0,0] = P[1,1] = 1/2 + 1/2*exp(-lambda*t)					*/
/*	P[0,1] = P[1,0] = 1/2 - 1/2*exp(-lambda*t)					*/
/*																*/
/*	Reversible model, so Pr(observing unordered pair ab) = 		*/
/*			(2-delta[a,b])*pi[a] P[a,b]^(2t)					*/
/*	where delta is the Kronecker delta and pi the equilibrium	*/
/*	frequency which is 1/2.										*/
/*																*/
/*	So Pr(observing unordered pair ab|mrca at t)				*/
/*			=	1/4 + 1/4*exp(-lambda*t)	if a=b				*/
/*			or	1/2 - 1/2*exp(-lambda*t)	otherwise			*/
/*																*/
/*	Expected pairwise diversity in a coalescent model, where	*/
/*	time is measured in units of PNe generations (P is ploidy	*/
/*	Ne is effective population size), is lambda/(1+2*lambda)	*/
/*																*/
/****************************************************************/
public:
	void set_defaults();
	FSM_Binary(Random *ran_in);
	FSM_Binary(const double lambda, Random *ran_in);
	FSM_Binary& update(const double lambda);
	int fast_mutate(const int state, const double time);
};

#endif // _MUTATION_H_
