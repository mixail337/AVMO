#include <iostream>
#include "drobs.h"
#include <vector>
#include <cmath>

using namespace std;

ostream& operator <<(ostream& out, const drobi& dr) {
	if (dr.chisl % dr.znam == 0) out << dr.chisl / dr.znam;
	else {
		out << dr.chisl << "/" << dr.znam;
	}
	return out;
}

bool operator>(drobi f, drobi s) {
	int nok = drobi::NOK(f.znam, s.znam);
	f.chisl *= (nok / f.znam);
	s.chisl *= (nok / s.znam);
	return f.chisl > s.chisl;
}


bool operator<(drobi f, drobi s) {
	int nok = drobi::NOK(f.znam, s.znam);
	f.chisl *= (nok / f.znam);
	s.chisl *= (nok / s.znam);
	return f.chisl < s.chisl;
}

bool operator==(drobi f, drobi s) {
	int nok = drobi::NOK(f.znam, s.znam);
	f.chisl *= (nok / f.znam);
	s.chisl *= (nok / s.znam);
	return f.chisl == s.chisl;
}

bool operator!=(drobi f, drobi s) {
	return !(f == s);
}

bool operator^(drobi f, drobi s) {
	int nok = drobi::NOK(f.znam, s.znam);
	f.chisl *= (nok / f.znam);
	s.chisl *= (nok / s.znam);
	return abs(f.chisl) < abs(s.chisl);
}


drobi operator*(drobi frst, drobi scnd) {
	return drobi(frst.chisl * scnd.chisl, frst.znam * scnd.znam);
}

drobi operator-(drobi frst, drobi scnd) {
	int nok = drobi::NOK(frst.znam, scnd.znam);
	frst.chisl *= (nok / frst.znam);
	scnd.chisl *= (nok / scnd.znam);
	return drobi(frst.chisl - scnd.chisl, nok);
}

drobi operator+(drobi frst, drobi scnd) {
	int nok = drobi::NOK(frst.znam, scnd.znam);
	frst.chisl *= (nok / frst.znam);
	scnd.chisl *= (nok / scnd.znam);
	return drobi(frst.chisl + scnd.chisl, nok);
}

drobi operator/(drobi frst, drobi scnd) {
	return drobi(scnd.znam, scnd.chisl) * frst;
}

istream& operator >> (istream& in,drobi& dr) {
	in >> dr.chisl;
	return in;
}





void basisCO(vector<vector<drobi>> matrix, vector<drobi> Z, vector<int>& op_solution, int& ind, int iter) {
	vector<drobi> co;
	drobi min = drobi(INT_MAX);
	for (int i = 0; i < matrix.size(); i++) {
		co.push_back(matrix[i].back() / matrix[i][ind]);
	}
	for (int i = 0; i < co.size(); i++) {
		if (co[i] < min && co[i]>0) {
			min = co[i];
			ind = i;
		}
	}
	if (iter == 2) op_solution[iter] = ind + 1;
	else op_solution[iter] = ind;
}


void print_table(vector<vector<drobi>> matrix, vector<drobi> Z, vector<int> op_solution) {
	int k = 0;
	cout << "b/p\t1\t";
	for (int i = 0; i < Z.size() - 1; ++i) {
		cout << "x" << i + 1 << "\t";
	}
	cout << endl;
	for (int i = 0; i < matrix.size(); ++i) {
		if (op_solution[i]<0) cout << "-" << "\t";
		else cout << "x" << op_solution[i] + 1 << "\t";
		cout << matrix[i][matrix[i].size() - 1] << "\t";
		for (int j = 0; j < matrix[i].size() - 1; ++j) {
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
	cout << "Z\t" << Z[Z.size() - 1] << "\t";
	for (int i = 0; i < Z.size() - 1; ++i) {
		cout << Z[i] << "\t";
	}
	cout << endl << endl;
}

void simplex(vector<vector<drobi>>& matrix,vector<drobi>& Z, vector<int>& op_solution) {
	for (int i = 0; i < op_solution.size(); ++i) {
		int ind = op_solution[i];
		if (matrix[i][ind] != drobi(1) || ind == 4) {    //|| ind == 4
			drobi koef = matrix[i][ind];
			for (int j = 0; j < matrix[i].size(); ++j) {
				matrix[i][j] = matrix[i][j] / koef;
			}
			for (int k = 0; k < op_solution.size(); ++k) {
				if (k != i) {
					koef = matrix[k][ind];
					for (int h = 0; h < matrix[i].size(); ++h) {
						matrix[k][h] = matrix[k][h] - (matrix[i][h] * koef);
					}
				}
			}

			if (Z[ind] != drobi(0)) {
				koef = Z[ind];
				for (int j = 0; j < Z.size(); ++j) {
					Z[j] = Z[j] - (matrix[i][j] * koef);
				}
				if (Z[Z.size() - 1] <drobi(0)) {
					for (int j = 0; j < Z.size(); ++j) {
						Z[j] = Z[j] * drobi(-1);
					}
				}
			}
		}
	}
}

int op_isvalid(vector<vector<drobi>> matrix, vector<int> op_solution) {
	int k = 0;
	for (int i = 0; i < matrix.size(); i++) {
		if (matrix[i].back() > 0) k++;
	}
	return k;
}



void basis(vector<vector<drobi>>& matrix, vector<drobi>& Z, vector<int> op_solution_need) {
	for (int j = 0; j < Z.size(); j++) {
		Z[j] = drobi(0) - Z[j];
	}
	for (int i = 0; i < op_solution_need.size(); ++i) {
		int ind = op_solution_need[i];
		if ((matrix[i][ind] != drobi(1)) || ind == 4) {
			drobi koef = matrix[i][ind];
			for (int j = 0; j < matrix[i].size(); ++j) {
				matrix[i][j] = matrix[i][j] / koef;
			}
			for (int k = 0; k < op_solution_need.size(); ++k) {
				if (k != i) {
					koef = matrix[k][ind];
					for (int h = 0; h < matrix[i].size(); ++h) {
						matrix[k][h] = matrix[k][h] - (matrix[i][h] * koef);
					}
				}
			}

			if (Z[ind] != drobi(0)) {
				koef = Z[ind];
				for (int j = 0; j < Z.size(); ++j) {
					Z[j] = Z[j] - (matrix[i][j] * koef);
				}
				if (Z[Z.size() - 1] < drobi(0)) {
					for (int j = 0; j < Z.size(); ++j) {
						Z[j] = Z[j] * drobi(-1);
					}
				}
			}
		}
	}
}

void search_basis(vector<vector<drobi>>& matrix, vector<int>& op_solution, vector<drobi>& Z, vector<int> op_solution_need) {
	vector<drobi> co;
	int i = 0;
	for (int j = 0; j < Z.size(); j++) {
		Z[j] = drobi(0) - Z[j];
	}
	for (int l = 0; l < op_solution.size(); ++l) {
		i = l;
		basisCO(matrix, Z, op_solution, i, l);
		int ind = op_solution_need[l];
		if (matrix[i][ind] != drobi(1)) {
			drobi koef = matrix[i][ind];
			for (int j = 0; j < matrix[i].size(); ++j) {
				matrix[i][j] = matrix[i][j] / koef;
			}
			for (int k = 0; k < op_solution.size(); ++k) {
				if (k != i) {
					koef = matrix[k][ind];
					for (int h = 0; h < matrix[i].size(); ++h) {
						matrix[k][h] = matrix[k][h] - (matrix[i][h] * koef);
					}
				}
			}
			if (Z[ind] != drobi(0)) {
				koef = Z[ind];
				for (int j = 0; j < Z.size(); ++j) {
					Z[j] = Z[j] - (matrix[i][j] * koef);
				}
				if (Z[Z.size() - 1] < drobi(0)) {
					for (int j = 0; j < Z.size(); ++j) {
						Z[j] = Z[j] * drobi(-1);
					}
				}
			}
		}
		print_table(matrix, Z, op_solution);
	}
}

bool needSimplex(vector<drobi> Z,int& ind) {
	drobi min(0);
	int minInd = -1;
	for (int i = 0; i < Z.size()-1; i++) {
		if (Z[i] > min) {
			min = Z[i];
			cout << Z[i] << " ";
			minInd = i;
		}
	}
	cout << endl;
	if (minInd == -1) return false;
	else {
		ind = minInd;
		return true;
	}
}

int CO(vector<vector<drobi>> matrix, vector<drobi> Z, vector<int>& op_solution, int ind) {
	vector<drobi> co;
	drobi min = drobi(INT_MAX);
	int index = -1;
	for (int i = 0; i < matrix.size(); i++) {
		co.push_back(matrix[i].back() / matrix[i][ind]);
		cout << co[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < co.size(); i++) {
		if (co[i] > drobi(0) && co[i] ^ min) {
			min = co[i];
			index = i;
		}
	}
	if (index == -1) {
		cout << "No solution";
		return -1;
	}
	cout << "CO\t" <<  min << endl;
	cout << "Target: x" << op_solution[index] + 1 << "change x" << ind + 1 << endl;
	op_solution[index] = ind;
	cout << endl;
	return 1;
}

void search_op(vector<vector<drobi>>& matrix, vector<drobi>& Z, vector<int>& op_solution) {
	int ind = 0;
	if (!op_isvalid(matrix, op_solution)) {
		CO(matrix, Z, op_solution,ind);
		simplex(matrix, Z, op_solution);
	}
}

void print_simplex(vector<vector<drobi>>& matrix, vector<drobi>& Z, vector<int>& op_solution) {
	int ind = 0;
	while (needSimplex(Z, ind)) {
		if (CO(matrix, Z, op_solution, ind) ==-1) break;
		simplex(matrix, Z, op_solution);
		cout << "b/p\t1\t";
		for (int i = 0; i < Z.size() - 1; ++i) {
			cout << "x" << i + 1 << "\t";
		}
		cout << endl;
		for (int i = 0; i < matrix.size(); ++i) {
			cout << "x" << op_solution[i] + 1 << "\t";
			cout << matrix[i][matrix[i].size() - 1] << "\t";
			for (int j = 0; j < matrix[i].size() - 1; ++j) {
				cout << matrix[i][j] << "\t";
			}
			cout << endl;
		}
		cout << "Z\t" << Z[Z.size() - 1] << "\t";
		for (int i = 0; i < Z.size() - 1; ++i) {
			cout << Z[i] << "\t";
		}
		cout << endl << endl;
	}
}

void print_answer(vector<vector<drobi>> matrix,vector<drobi> Z,vector<int> op_solution, vector<drobi> Z1) {
	vector<drobi> ans;
	drobi answer;
	for (int i = 0; i < Z.size() - 1; i++) {
		ans.push_back(drobi(0));
	}

	cout << "X = (";
	for (int i = 0; i < op_solution.size(); i++) {
		ans[op_solution[i]] = matrix[i].back();
	}
	for (int i = 0; i < ans.size(); i++) {
		cout << ans[i] << " ";
	}
	cout << ")";
	cout << endl;
	cout << "Z = ";
	for (int i = 0; i < Z1.size() - 1; i++) {
		answer = answer + Z1[i] * ans[i];
	}
	cout << answer;
}






int main()
{
	vector<int> op_solution = { -1,-1,-1 };
	vector<int> op_solution_need = { 1,0,3 };
	vector<drobi> Z = { 6,1,0,0,0,0 };
	vector<vector<drobi>> matrix = {
	    {3,7,-1,0,0,36},
		{3,1,0,-1,0,14},
		{2,1,0,0,-1,13}
	};
	print_table(matrix, Z, op_solution);
	//search_basis(matrix, op_solution, Z, op_solution_need);
	basis(matrix, Z, op_solution_need);
	print_table(matrix, Z, op_solution_need);
	vector<drobi> Z1 = { 6,1,0,0,0,0 };
	print_simplex(matrix, Z, op_solution_need);
	print_answer(matrix, Z, op_solution_need, Z1);
}

		