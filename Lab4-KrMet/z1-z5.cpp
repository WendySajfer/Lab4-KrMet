#include <iostream>
#include "windows.h"
#include <vector>
#include "string"
#include <algorithm>
#include <cmath>

using namespace std;

class Translate {
private:
	//vector<char> Alfabet = { 'А', 'Б', 'В', 'Г', 'Д', 'Е', 'Ё', 'Ж', 'З', 'И', 'Й', 'К', 'Л', 'М', 'Н', 'О', 'П', 'Р', 'С', 'Т', 'У', 'Ф', 'Х', 'Ц', 'Ч', 'Ш', 'Щ', 'Ъ', 'Ы', 'Ь', 'Э', 'Ю', 'Я', 'а', 'б', в, 'г', 'д', 'е', 'ё', 'ж', 'з', 'и', 'й', 'к', 'л', 'м', 'н', 'о', 'п', 'р', 'с', 'т', 'у', 'ф', 'х', 'ц', 'ч', 'ш', 'щ', 'ъ', 'ы', 'ь', 'э', 'ю', 'я', ' ', '.', ',', '!', '?'};
	vector<char> Alfabet = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', ' '};
	int Alfabet_size = Alfabet.size();
	vector<vector<int>> Blocks_Text;
	vector<vector<int>> Cipher_Blocks_Text;
	vector<vector<char>> Cipher_Text_char;
	vector<vector<char>> Blocks_Text_char;
	vector<vector<int>> Composition_module(vector<vector<int>> A, vector<vector<int>> B) {
		if (A[0].size() != B.size()) return {};
		vector<vector<int>> C;
		vector<int>::iterator it;
		for (int i = 0; i < A.size(); i++) {
			vector<int> comp;
			for (int j = 0; j < B[0].size(); j++) {
				int buf_ij = 0;
				for (int k = 0; k < A[0].size(); k++) {
					buf_ij += A[i][k] * B[k][j]; // c_ij = E(a_ik*b_kj); (E - сумма)
				}
				comp.push_back(buf_ij);
			}
			C.push_back(comp);
		}
		for (int i = 0; i < C.size(); i++) {
			for (int j = 0; j < C[0].size(); j++) {
				C[i][j] = C[i][j] % Alfabet_size;
			}
		}
		return C;
	}
public:
	string Encryption(string Text, vector<vector<int>> Key) {
		int size_key = Key.size();
		int blocks = Text.size() / size_key;
		int end_block = Text.size() % size_key;
		string Cipher_Text = "";
		for (int i = 0; i < blocks; i++) {
			vector<char> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(Text[i * (size_key)+j]);
			}
			Blocks_Text_char.push_back(buf_str);
		}
		if (end_block != 0) {
			vector<char> buf_str;
			for (int j = 0; j < end_block; j++) {
				buf_str.push_back(Text[Text.size() - end_block + j]);
			}
			for (int j = 0; j < size_key - end_block; j++) {
				buf_str.push_back(' ');
			}
			Blocks_Text_char.push_back(buf_str);
		}
		int buf = 0;
		for (int i = 0; i < Blocks_Text_char.size(); i++) {
			vector<int> buf_int;
			for (int j = 0; j < size_key; j++) {
				buf = find(Alfabet.begin(), Alfabet.end(), Blocks_Text_char[i][j]) - Alfabet.begin();
				buf_int.push_back(buf);
			}
			Blocks_Text.push_back(buf_int);
		}
		for (int i = 0; i < Blocks_Text.size(); i++) {
			vector<vector<int>> buf;
			for (int j = 0; j < size_key; j++) {
				buf.push_back({ Blocks_Text[i][j] });
			}
			vector<vector<int>> buf_matrix = Composition_module(Key, buf);
			vector<int> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(buf_matrix[j][0]);
			}
			Cipher_Blocks_Text.push_back(buf_str);
		}

		for (int i = 0; i < Cipher_Blocks_Text.size(); i++) {
			vector<char> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(Alfabet[Cipher_Blocks_Text[i][j]]);
			}
			Cipher_Text_char.push_back(buf_str);
		}
		
		for (int i = 0; i < Cipher_Text_char.size(); i++) {
			for (int j = 0; j < size_key; j++) {
				Cipher_Text = Cipher_Text + Cipher_Text_char[i][j];
			}
		}
		cout << "Cipher_Text: " << Cipher_Text << endl;
		return Cipher_Text;
	}
	string Decryption(string Cipher_Text, vector<vector<int>> Key_Inv) {
		int size_key = Key_Inv.size();
		int blocks = Cipher_Text.size() / size_key;
		string Text = "";
		for (int i = 0; i < blocks; i++) {
			vector<char> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(Cipher_Text[i * (size_key)+j]);
			}
			Cipher_Text_char.push_back(buf_str);
		}
		int buf = 0;
		for (int i = 0; i < Cipher_Text_char.size(); i++) {
			vector<int> buf_int;
			for (int j = 0; j < size_key; j++) {
				buf = find(Alfabet.begin(), Alfabet.end(), Cipher_Text_char[i][j]) - Alfabet.begin();
				buf_int.push_back(buf);
			}
			Cipher_Blocks_Text.push_back(buf_int);
		}
		for (int i = 0; i < Cipher_Blocks_Text.size(); i++) {
			vector<vector<int>> buf;
			for (int j = 0; j < size_key; j++) {
				buf.push_back({ Cipher_Blocks_Text[i][j] });
			}
			vector<vector<int>> buf_matrix = Composition_module(Key_Inv, buf);
			vector<int> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(buf_matrix[j][0]);
			}
			Blocks_Text.push_back(buf_str);
		}

		for (int i = 0; i < Blocks_Text.size(); i++) {
			vector<char> buf_str;
			for (int j = 0; j < size_key; j++) {
				buf_str.push_back(Alfabet[Blocks_Text[i][j]]);
			}
			Blocks_Text_char.push_back(buf_str);
		}

		for (int i = 0; i < Blocks_Text_char.size(); i++) {
			for (int j = 0; j < size_key; j++) {
				Text = Text + Blocks_Text_char[i][j];
			}
		}
		cout << "Text: " << Text << endl;
		return Text;
	}
};
class Functions {
private:
	vector<int> Simple_Numbers;
	vector<int> R_numbers;
	vector<int> Q_numbers;
	vector<int> u_v = { 0,0 };
	int f_Eiler = 0;
	int GCD = 1;
public:
	vector<int> Simple_numbers_Decomposition(int n) {
		int n1 = n;
		for (int i = 2; i <= sqrt(n1); i++) {
			while (n % i == 0) {
				Simple_Numbers.push_back(i);
				n /= i;
			}
		}
		if (n != 1) Simple_Numbers.push_back(n);
		return Simple_Numbers;
	}
	int NodEuclidean(int a, int b) {
		if (a == 0 || b == 0) {
			return 0;
		}
		int D, d, r, q;
		D = a;
		d = b;
		while (true) {
			r = D % d;
			q = (D - r) / d;
			if (r < 0) { //проверка отрицательных a, b
				if (d > 0) {
					q--;
					r += d;
				}
				if (d < 0) {
					q++;
					r -= d;
				}
			}
			Q_numbers.push_back(q);
			if (r == 0) break;
			R_numbers.push_back(r);
			D = d;
			d = r;
		}
		GCD = d;
		return GCD;
	}

	void Recursion_Bezu(int ux, int vx, int uy, int vy, int q_index) {
		//Основано на формуле: r(n) = a(u(n-2) - u(n-1)q(n)) + b(v(n-2) - v(n-1)q(n)) (для n>1)
		int buf_ux = ux, buf_vx = vx;
		ux = uy;
		vx = vy;
		uy = buf_ux - ux * Q_numbers[q_index];
		vy = buf_vx - vx * Q_numbers[q_index];
		q_index++;
		if (q_index + 1 == Q_numbers.size()) {
			u_v = { uy, vy };
		}
		else Recursion_Bezu(ux, vx, uy, vy, q_index); //рекурсия
	}
	vector<int> IdentityBezu(int a, int b) {
		int a1 = 1, u = 1, v = 1;
		if (GCD != b) {
			if (GCD == a) v = 0;
			else {
				v *= -(Q_numbers[Q_numbers.size() - 2]);
				if (GCD != a * u + b * v) {
					v = -(Q_numbers[0]);
					Recursion_Bezu(0, 1, 1, v, 1);
					u = u_v[0];
					v = u_v[1];
				}
			}
		}
		return { u, v };
	}
	int Eiler_Method4(int n) {
		Simple_numbers_Decomposition(n);
		double f = n, buf;
		for (int i = 0; i < Simple_Numbers.size(); i++) {
			buf = Simple_Numbers[i];
			f *= (1 - (1 / buf));
			if (i < Simple_Numbers.size() - 1) {
				while (true && i + 1 < Simple_Numbers.size()) {
					if (Simple_Numbers[i] == Simple_Numbers[i + 1]) Simple_Numbers.erase(Simple_Numbers.begin() + i + 1);
					else break;
				}
			}
		}
		f_Eiler = f;
		return f_Eiler;
	}
};

int InverseBezu(int m, int n) {
	int inv_n = 0;
	Functions F1;
	int GCD = F1.NodEuclidean(m, n);
	if (GCD != 1) return 0;
	vector<int> U_V = F1.IdentityBezu(m, n);
	inv_n = U_V[1] % m;
	if (inv_n < 0) inv_n = m + inv_n;
	return inv_n;
}

class Matrix {
private:
	vector<vector<int>> matrix;
	vector<int> width_form;
	int m = 0, n = 0;

	void Format() {
		width_form.clear();
		int width, buf_width, width_null;
		string str_width;
		width = 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				str_width = to_string(matrix[j][i]);
				for (int i = str_width.size() - 1; i >= 0; i--) {
					if (str_width[i] == '0' || str_width[i] == ',') {
						str_width.erase(i);
					}
					else break;
				}
				buf_width = str_width.size();
				if (width < buf_width) width = buf_width;
			}
			width_form.push_back(width);
		}
	}
	void Size() {
		m = matrix.size();
		n = matrix[0].size();
		for (int i = 1; i < m; i++) {
			if (matrix[i].size() != n) {
				n = 0;
				break;
			}
		}
	}
public:
	void Input_Matrix(int m, int n)
	{
		int number;
		vector<int> str_matrix;
		cout << "Enter the [" << m << "," << n << "] matrix:" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> number;
				str_matrix.push_back(number);
			}
			matrix.push_back(str_matrix);
			str_matrix.clear();
		}
	}
	void Output_Matrix() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	vector<vector<int>> get_Matrix() {
		return matrix;
	}
	void set_matrix(vector<vector<int>> matrix_copy) {
		matrix = matrix_copy;
	}
};

class Determinant {
private:
	vector<vector<int>> matrix_det;
	int det = 0;
	int Det_2x(vector<vector<int>> matrix_2x) {
		int right = matrix_2x[0][0] * matrix_2x[1][1];
		int left = matrix_2x[0][1] * matrix_2x[1][0];
		return right - left;
	}
	int Det_3x(vector<vector<int>> matrix_3x) {
		int right = 0;
		int buf = 1, buf1 = 1;
		for (int i = 0; i < 3; i++) {
			buf *= matrix_3x[i][i];
		}
		right += buf;
		buf = 1;
		for (int i = 0; i < 3; i++) {
			if (i + 1 < 3) {
				buf *= matrix_3x[i][i + 1];
				buf1 *= matrix_3x[i + 1][i];
			}
			else {
				buf *= matrix_3x[i][0];
				buf1 *= matrix_3x[0][i];
			}
		}
		right += buf;
		right += buf1;
		int left = 0;
		buf = matrix_3x[0][2] * matrix_3x[2][0] * matrix_3x[1][1];
		left += buf;
		buf = matrix_3x[0][1] * matrix_3x[1][0] * matrix_3x[2][2];
		left += buf;
		buf = matrix_3x[1][2] * matrix_3x[2][1] * matrix_3x[0][0];
		left += buf;
		return right - left;
	}
	int Det_mx(vector<vector<int>> matrix_mx) {
		if (matrix_mx.size() == 1) return matrix_mx[0][0];
		if (matrix_mx.size() == 2) return Det_2x(matrix_mx);
		else if (matrix_mx.size() == 3) return Det_3x(matrix_mx);
		else {
			int det_mx = 0;
			for (int i = 0; i < matrix_mx.size(); i++) {
				vector<vector<int>> buf_m = matrix_mx;
				buf_m.erase(buf_m.begin());
				int buf_size = buf_m.size();
				for (int buf_index = 0; buf_index < buf_size; buf_index++) {
					buf_m[buf_index].erase(buf_m[buf_index].begin() + i);
				}
				if (i % 2 == 0) {
					det_mx += matrix_mx[0][i] * Det_mx(buf_m);
				}
				else det_mx -= matrix_mx[0][i] * Det_mx(buf_m);
			}
			return det_mx;
		}
	}
public:
	void Input_Det_Matrix(vector<vector<int>> matrix) {
		matrix_det = matrix;
	}
	void Decision() {
		det = Det_mx(matrix_det);
	}
	double get_det() {
		return det;
	}
};

class Inverse_of_module {
private:
	vector<vector<int>> matrix;
	Matrix Inv_m;
	vector<vector<int>> algebra_add;
	vector<vector<int>> inv_matrix;
	int inverse_det = 0, module = 2;
	void Create_Algebra() {
		Determinant D;
		int algebra = 0;
		for (int i = 0; i < matrix.size(); i++) {
			vector<int> buf_str_algebra;
			for (int j = 0; j < matrix[i].size(); j++) {
				vector<vector<int>> buf = matrix;
				for (int k = 0; k < buf.size(); k++) {
					vector<int>::iterator it = buf[k].begin() + j;
					buf[k].erase(it);
				}
				buf.erase(buf.begin() + i);
				D.Input_Det_Matrix(buf);
				D.Decision();
				algebra = D.get_det();
				if ((i + j) % 2 == 1) {
					algebra = -algebra;
				}
				buf_str_algebra.push_back(algebra);
			}
			algebra_add.push_back(buf_str_algebra);
		}
	}
	void Inverse_Matrix() {
		bool flag_deny = false;
		for (int i = 0; i < algebra_add.size(); i++) { //1.Делим по модулю матрицу алгебраических дополнений(!не теряя знаки)
			for (int j = 0; j < algebra_add[i].size(); j++) {
				flag_deny = (algebra_add[i][j] < 0);
				algebra_add[i][j] = algebra_add[i][j] % module;
			}
		}
		for (int i = 0; i < algebra_add.size(); i++) { //2.Умножаем матрицу алгебраических дополнений на обратный детерминанту элемент.
			for (int j = 0; j < algebra_add[i].size(); j++) {
				algebra_add[i][j] *= inverse_det;
			}
		}
		for (int i = 0; i < algebra_add.size(); i++) { //3.Делим данную матрицу по модулю(!не теряя знаки)
			for (int j = 0; j < algebra_add[i].size(); j++) {
				flag_deny = (algebra_add[i][j] < 0);
				algebra_add[i][j] = algebra_add[i][j] % module;
			}
		}
		inv_matrix = algebra_add;
		for (int i = 0; i < inv_matrix.size(); i++) { //4.Транспонируем ее
			for (int j = 0; j < inv_matrix[i].size(); j++) {
				inv_matrix[i][j] = algebra_add[j][i];
			}
		}
		for (int i = 0; i < inv_matrix.size(); i++) { //5.Отрицательные эементы меняем на противопоожные положительные
			for (int j = 0; j < inv_matrix[i].size(); j++) {
				flag_deny = (inv_matrix[i][j] < 0);
				if (flag_deny) inv_matrix[i][j] += module;
			}
		}
		//Если перемножить матрицу ключа и эту матрицу, а потом результат разделить по модулю, мы получим единичную матрицу
	}
public:
	void Input(vector<vector<int>> matrix_copy, int inv_det, int m) {
		matrix = matrix_copy;
		inverse_det = inv_det;
		module = m;
	}

	void Inverse_Decision() {
		Create_Algebra();
		Inverse_Matrix();
		Inv_m.set_matrix(inv_matrix);
		cout << "Key^(-1) = " << endl;
		Inv_m.Output_Matrix();
	}
	vector<vector<int>> get_inv() {
		return inv_matrix;
	}
};

class Swap_Group {
private:
	vector<int> index_swap;
	vector<vector<int>> begin_swap;
	vector<vector<int>> e_f;
	vector<vector<vector<int>>> Group;

	vector<vector<int>> Composition_Swap(vector<vector<int>> A_swap, vector<vector<int>> B_swap) {
		vector<vector<int>> C_swap;
		vector<int> comp_swap;
		vector<int>::iterator it;
		int buf, buf_index;
		C_swap.push_back(index_swap);
		for (int i = 0; i < index_swap.size(); i++) {
			buf = B_swap[1][i];
			it = find(A_swap[0].begin(), A_swap[0].end(), buf);
			buf_index = it - A_swap[0].begin();
			buf = A_swap[1][buf_index];
			comp_swap.push_back(buf);
		}
		C_swap.push_back(comp_swap);
		return C_swap;
	}
public:
	void Input(vector<vector<int>> copy_swap) {
		begin_swap = copy_swap;
		index_swap = copy_swap[0];
	}
	void Output(vector<vector<int>> out_swap) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < index_swap.size(); j++) {
				cout << out_swap[i][j] << " ";
			}
			cout << endl;
		}
	}
	void Group_Create() {
		e_f.push_back(index_swap);
		e_f.push_back(index_swap);
		Group.push_back(begin_swap);
		cout << "f" << endl;
		Output(Group[0]);
		int index = 0;
		while (true) {
			Group.push_back(Composition_Swap(Group[index], begin_swap));
			index++;
			cout << "f^" << index + 1 << endl;
			Output(Group[index]);
			if (Group[index] == e_f) break;
		}
	}
	void Result() {
		cout << "<f> = {e, f";
		for (int i = 2; i <= Group.size() - 1; i++) {
			cout << ", f^" << i;
		}
		cout << "}" << endl;
		cout << "|<f>| = " << Group.size() << endl;
		cout << "f^(-1) = f^" << Group.size() - 1 << endl;
		Output(*(Group.end() - 2));
	}
};



class Tasks {
private:
	int f_Eiler = 0, inv_det = 0, mod = 26;//71
	Functions F;
	vector<int> Multi_elements;
	vector<int> Deg_elements;
	vector<int> Antiderivative_elements;
	vector<vector<int>> Key, Key_inv;
	bool Inverse_determinant() {
		Determinant Det;
		Det.Input_Det_Matrix(Key);
		Det.Decision();
		int det = Det.get_det();
		if (det == 0) return false;
		inv_det = InverseBezu(mod, det);
		if (inv_det == 0) return false;
		return true;
	}
	bool Inverse_Key() {
		bool flag = true;
		int m = Key.size();
		for (int i = 0; i < m; i++) {
			if (Key[i].size() != m) {
				flag = false;
				break;
			}
		}
		if (!flag) {
			cout << "Error! Key isn't the a square matrix." << endl;
			return false;
		}
		flag = Inverse_determinant();
		if (flag == 0) {
			cout << "Error! Deteminant of key haven't the inverse element." << endl;
			return false;
		}
		Inverse_of_module Inv;
		Inv.Input(Key, inv_det, mod);
		Inv.Inverse_Decision();
		Key_inv = Inv.get_inv();
		return true;
	}
	bool Group_degrees(int m, int g) {
		int buf_g = g;
		bool group_flag = false;
		Deg_elements.push_back(g);
		while (true) {
			buf_g *= g;
			buf_g = buf_g % m;
			if (buf_g != g) {
				Deg_elements.push_back(buf_g);
			}
			if (buf_g == 1) {
				group_flag = true;
				break;
			}
			if (buf_g == 0) break;
			if (buf_g == Deg_elements[0]) break;
		}
		return group_flag;
	}
	bool Antideriv_el(int m, int g) {
		bool group_flag = Group_degrees(m, g);
		if (Deg_elements.size() == f_Eiler && group_flag) {
			Deg_elements.clear();
			Antiderivative_elements.push_back(g);
			return true;
		}
		Deg_elements.clear();
		if (group_flag) return true;
		else return false;
	}

	int Multiplicative_elements_size(int m) {
		int numbers_multi_size = 0;
		for (int i = 1; i < m; i++) {
			int GCD = F.NodEuclidean(m, i);
			if (GCD == 1) {
				Multi_elements.push_back(i);
				numbers_multi_size++;
			}
		}
		return numbers_multi_size;
	}
public:
	void Multiplicative_elements(int m) {
		if (m > 46341 || m < 2) { //m < sqrt(int_max) + 1 for *
			cout << "Invalid m." << endl;
			return;
		}
		f_Eiler = F.Eiler_Method4(m);
		int numbers_multi_el = Multiplicative_elements_size(m);
		cout << "Elements:";
		for (int i = 0; i < numbers_multi_el; i++) {
			cout << ' ' << Multi_elements[i];
		}
		cout << endl;
		if (numbers_multi_el == f_Eiler) {
			cout << numbers_multi_el << "==" << f_Eiler << endl;
		}
		else cout << numbers_multi_el << "!=" << f_Eiler << endl;
	}
	void Degrees_elements(int m, int g) {
		if (m > 46341 || m < 1) { //m < sqrt(int_max) + 1 for *
			cout << "Invalid m." << endl;
			return;
		}
		if (g >= m || g < 0) {
			cout << "Invalid element g." << endl;
			return;
		}
		f_Eiler = F.Eiler_Method4(m);
		bool group_flag = Group_degrees(m, g);
		int size_G = Deg_elements.size();
		cout << "Elements_degree:";
		for (int i = 0; i < size_G; i++) {
			cout << ' ' << Deg_elements[i];
		}
		cout << endl;
		if (group_flag) {
			cout << "|G| = " << size_G << endl;
			if (size_G == f_Eiler) {
				cout << g << " - is the antiderivative root module " << m << endl;
			}
		}
	}
	void Cyclic_Group(int m) {
		if (m > 46341 || m < 1) { //m < sqrt(int_max) + 1 for *
			cout << "Invalid m." << endl;
			return;
		}
		f_Eiler = F.Eiler_Method4(m);
		Multiplicative_elements_size(m);
		bool cyclic_flag = true;
		for (int i = 0; i < Multi_elements.size(); i++) {
			if (cyclic_flag) cyclic_flag = Antideriv_el(m, Multi_elements[i]);
			else break;
		}
		if (Antiderivative_elements.size() == 0) cyclic_flag = false;
		if (cyclic_flag) {
			cout << "The group is cyclic." << endl;
			cout << "Antiderivative elements: ";
			for (int i = 0; i < Antiderivative_elements.size(); i++) {
				cout << ' ' << Antiderivative_elements[i];
			}
			cout << endl;
		}
		else cout << "The group isn't cyclic." << endl;
	}
	void Group_of_Swap(int n) {
		if (n < 2) {
			cout << "Invalid n." << endl;
			return;
		}
		Swap_Group S_G;
		vector<vector<int>> swap;
		vector<int> buf;
		vector<int> sort_vector;
		cout << "Input swap vector:" << endl;
		int buf_number = 0;
		for (int i = 0; i < n; i++) {
			cin >> buf_number;
			buf.push_back(buf_number);
		}
		sort_vector = buf;
		sort(sort_vector.begin(), sort_vector.end());
		swap.push_back(sort_vector);
		swap.push_back(buf);
		S_G.Input(swap);
		S_G.Group_Create();
		S_G.Result();
	}
	void Encryption_of_Hill(int n) {
		Matrix Key_Hill;
		Translate Hill;
		string Text;
		Key_Hill.Input_Matrix(n, n);
		Key = Key_Hill.get_Matrix();
		bool flag = true;
		int m = Key.size();
		for (int i = 0; i < m; i++) {
			if (Key[i].size() != m) {
				flag = false;
				break;
			}
		}
		if (!flag) {
			cout << "Error! Key isn't the a square matrix." << endl;
			return;
		}
		flag = Inverse_determinant();
		if (flag == 0) {
			cout << "Error! Deteminant of key haven't the inverse element." << endl;
			return;
		}
		cin.ignore(1000, '\n');
		cout << "Input message: " << endl;
		getline(cin, Text);
		Hill.Encryption(Text, Key);
	}
	void Decryption_of_Hill(int n) {
		Matrix Key_Hill;
		Translate Hill;
		string Cipher_Text;
		Key_Hill.Input_Matrix(n, n);
		Key = Key_Hill.get_Matrix();
		bool flag = Inverse_Key();
		if (!flag) return;
		cin.ignore(1000, '\n');
		cout << "Input message: " << endl;
		getline(cin, Cipher_Text);
		if (Cipher_Text.size() % n != 0) {
			cout << "Error! Cipher must be divided into blocks of n." << endl;
			return;
		}
		Hill.Decryption(Cipher_Text, Key_inv);
	}
};


int main() {
	setlocale(LC_ALL, "Rus");
	SetConsoleCP(1251);
	int m, g, n, task, exit = 1;
	while (exit == 1) {
		Tasks Math;
		cout << "1.Multiplicative elements" << endl << "2.Degrees elements" << endl << "3.Cyclic Group" << endl << "4.Group of Swap" << endl << "5.Encryption of Hill" << endl << "6.Decryption of Hill" << endl << "7.Exit" << endl << "Choose a way:" << endl;
		cin >> task;
		switch (task)
		{
		case 1:
			cout << "Input module:" << endl << "m = ";
			cin >> m;
			Math.Multiplicative_elements(m);
			break;
		case 2:
			cout << "Input module:" << endl << "m = ";
			cin >> m;
			cout << "Input element:" << endl << "g = ";
			cin >> g;
			Math.Degrees_elements(m, g);
			break;
		case 3:
			cout << "Input module:" << endl << "m = ";
			cin >> m;
			Math.Cyclic_Group(m);
			break;
		case 4:
			cout << "Input numbers of swap:" << endl << "n = ";
			cin >> n;
			Math.Group_of_Swap(n);
			break;
		case 5:
			cout << "Input size of key:" << endl << "n = ";
			cin >> n;
			Math.Encryption_of_Hill(n);
			break;
		case 6:
			cout << "Input size of key:" << endl << "n = ";
			cin >> n;
			Math.Decryption_of_Hill(n);
			break;
		case 7:
			exit = 0;
			break;
		default:
			cout << "This task does not exist" << endl;
			break;
		}
	}
	system("pause");
	return 0;
}