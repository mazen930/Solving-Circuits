#include<iostream>
#include<fstream>
#include<Eigen\Eigen>
using namespace std;
using namespace Eigen;
#include <complex>   
#include <string>

#define pi 3.14159265358979

void printpolar(complex<double>x)
{
	double c = arg(x) * 180 / pi;
	cout << abs(x) << "/_" << c << endl;
}
complex<double> zind(double w, double l)
{
	double z = w*l;
	complex<double> ind(0.0, z);
	return ind;
}
complex<double> zcap(double w, double c)
{
	double z = -1 / (w*c);
	complex<double> cap(0.0, z);
	return cap;
}
complex<double> returncomplex(double ind, double pol)
{
	complex<double> complex = polar(ind, pol*pi / 180);
	return complex;

}

int main()
{


	struct Branch
	{
		string Name;
		int start;
		int end;
		complex<double> Imedance;
	};
	struct Branch branches[20];
	complex<double> Impedance[30];
	ifstream input_file;
	int i = 0; int node = 0; int extra_eq = 0; int a = 0; int q = 0; int w = 4000;
	string name;
	cout << "Enter the name of the file" << endl;
	cin >> name;
	double temp1 = 0, temp2 = 0;
	bool check, flag = false, flag_2nd = true;
	input_file.open(name);
	while (!input_file.eof())
	{
		input_file >> branches[i].Name;

		if (branches[i].Name[0] == 'V' || branches[i].Name[0] == 'I')
		{
			input_file >> branches[i].start;
			input_file >> branches[i].end;
			input_file >> temp1;
			input_file >> temp2;
			branches[i].Imedance = returncomplex(temp1, temp2);

		}
		else
		{
			if (branches[i].Name[0] == 'R')
			{
				input_file >> branches[i].start;
				input_file >> branches[i].end;
				input_file >> temp1;
				branches[i].Imedance = temp1;
			}
			else if (branches[i].Name[0] == 'C')
			{
				input_file >> branches[i].start;
				input_file >> branches[i].end;
				input_file >> temp1;
				temp2 = -1 / (temp1 * w);
				branches[i].Imedance = complex<double>(0.0, temp2);
			}
			else if (branches[i].Name[0] == 'L')
			{
				input_file >> branches[i].start;
				input_file >> branches[i].end;
				input_file >> temp1;
				temp2 = temp1 * w;
				branches[i].Imedance = complex<double>(0.0, temp2);
			}


		}
		i++;
	}
	int g = 0;
	for (int i = 0; i < 20; i++)
	{
		if (branches[i].Name[0] == 'V')
		{
			struct Branch temp = branches[i];
			branches[i] = branches[g];
			branches[g] = temp;
			g++;

		}
	}
	for (int i = 0; i < 20; i++)
	{
		if (branches[i].Name[0] == 'I')
		{
			struct Branch temp = branches[i];
			branches[i] = branches[g];
			branches[g] = temp;
			g++;

		}
	}



	input_file.close();
	for (int j = 0; j < i; j++)
	{
		if (branches[j].start > node)
			node = branches[j].start;
		if (branches[j].end > node)
			node = branches[j].end;
	}
	//cout << node + 1 << endl;

	for (int j = 0; j < i; j++)//to increase number of equation in matrix
	{
		if (branches[j].Name[0] == 'V')
		{
			//if ((branches[j].start != 0) && (branches[j].end != 0))
			extra_eq++;
		}
	}
	int z = node + extra_eq;
	MatrixXcd A(z, z);
	MatrixXcd B(z, 1);
	MatrixXcd X(z, 1);
	for (int k = 0; k < z; k++)
	{
		for (int l = 0; l < z; l++)
		{
			A(k, l) = 0;

		}
	}
	for (int k = 0; k < z; k++)
	{
		B(k, 0) = 0;
		X(k, 0) = 0;
	}

	for (int k = 0; k < z; k++)
	{
		if (flag)
			break;
		flag_2nd = true;
		for (int l = 0; l < z; l++)
		{
			if (flag)
				break;
			if (!flag_2nd)
				break;
			for (int s = 0; s < i; s++)
			{
				if (k == l)
				{
					if (((branches[s].start == l + 1) || (branches[s].end == l + 1)) && (branches[s].Name[0] == 'R' || branches[s].Name[0] == 'L' || branches[s].Name[0] == 'C'))
					{
						A(k, l) += 1.0 / branches[s].Imedance;
					}
				}
				else if (((branches[s].start == k + 1) && (branches[s].end == l + 1)) && (branches[s].Name[0] == 'R' || branches[s].Name[0] == 'L' || branches[s].Name[0] == 'C'))
				{
					A(k, l) = A(l, k) = -(1.0 / branches[s].Imedance);
					break;
				}
				if (l >= node)
				{
					if (branches[s].Name[0] == 'V')
					{

						if (branches[s].start == (k + 1))
						{
							int h = branches[s].start - 1;
							int p = l + s;
							A(k, p) = -1; flag_2nd = false;
							a++;
						}
						if (branches[s].end == (k + 1))
						{
							int h = branches[s].end - 1;
							int p = l + s;
							A(k, p) = 1; flag_2nd = false;
							a++;
						}
					}
				}
				if (k >= node)
				{
					if (branches[s].Name[0] == 'V')
					{
						if ((branches[s].end == 0) || branches[s].start == 0)
						{
							int h;
							int p = k + s;
							if (branches[s].end == 0)
							{
								h = branches[s].start - 1;
								B(p, 0) = branches[s].Imedance;
							}
							if (branches[s].start == 0)
							{
								h = branches[s].end - 1;
								B(p, 0) = -branches[s].Imedance;
							}
							A(p, h) = 1;
							flag = true;
						}
						else
						{
							int f, g;
							f = branches[s].start - 1;
							g = branches[s].end - 1;
							int p = k + s;
							A(p, f) = 1;
							A(p, g) = -1;
							B(p, 0) = branches[s].Imedance;
							flag = true;
						}
					}
				}

			}

		}
	}

	for (int k = 0; k < i; k++)
	{
		if (branches[k].Name[0] == 'I')
		{
			if (branches[k].start)
			{
				int g = branches[k].start - 1;
				B(g, 0) += branches[k].Imedance;
			}
			if (branches[k].end)
			{
				int y = branches[k].end - 1;
				B(y, 0) -= branches[k].Imedance;
			}
		}

	}
	X = A.lu().solve(B);
	complex<double>Current[50];
	int m, r = 1;
	cout << "Voltage at Node (0) = 0" << endl;
	for (m = 0; m < (z - g); m++)
		cout << "Voltage at Node (" << m + 1 << ") = " << X(m, 0) << endl;
	for (int k = m; k < z; k++)
	{
		cout << "The current pass through Voltage Source " << r << " = " << X(k, 0) << endl;
		r++;
	}
	for (int k = 0; k < i; k++)
	{
		if (branches[k].Name[0] == 'R' || branches[k].Name[0] == 'C' || branches[k].Name[0] == 'L')
		{
			int u, v;
			if (branches[k].start == 0)
			{
				u = branches[k].end - 1;
				Current[k] = X(u, 0) / branches[k].Imedance;
			}
			else if (branches[k].end == 0)
			{
				v = branches[k].start - 1;
				Current[k] = -X(v, 0) / branches[k].Imedance;
			}
			else
			{
				u = branches[k].end - 1;
				v = branches[k].start - 1;
				Current[k] = (X(u, 0) - X(v, 0)) / branches[k].Imedance;
			}


			cout << "I(" << branches[k].end << " , " << branches[k].start << ")" << Current[k] << endl;
		}
		if (branches[k].Name[0] == 'I')
		{
			cout << "I(" << branches[k].end << " , " << branches[k].start << ")" << branches[k].Imedance << endl;
		}
	}
	/*for (int f = 0; f < z; f++)
	{
	for (int k = 0; k < z; k++)
	cout << A(f, k) << " ";
	cout << endl;
	}*/
	/*cout << "=============" << endl;
	cout << B << "\n\n";
	cout << "=============" << endl;*/
	return 0;

}