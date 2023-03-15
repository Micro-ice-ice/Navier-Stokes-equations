#include <iostream>
#include "vars.hpp"
#include "cell.hpp"
#include <vector>
#include <fstream>
#include <Eigen/Sparse>
#include "f.hpp"

using namespace std;

void InitMattrix(Eigen::SparseMatrix<double> &A, vector<Cell> &cells)
{
    double a_e = HT, a_w = HT, a_s = HT, a_n = HT;

    // внутриние ячейки
    for (int j = NX; j < SIZE - NX; j += NX)
    {
        for (int i = j + 1; i < j + NX - 1; ++i)
        {

            A.coeffRef(i, cells[i].NorthIndex) = -a_n;
            A.coeffRef(i, cells[i].WestIndex) = -a_w;
            A.coeffRef(i, i) = a_w + a_e + a_n + a_s;

            A.coeffRef(i, cells[i].EastIndex) = -a_e;

            A.coeffRef(i, cells[i].SouthIndex) = -a_s;
        }
    }

    // north
    for (int i = 1; i < NX - 1; ++i)
    {
        A.coeffRef(i, i) = a_w + a_e + a_s;
        A.coeffRef(i, cells[i].WestIndex) = -a_w;
        A.coeffRef(i, cells[i].EastIndex) = -a_e;
        // A.coeffRef(i, cells[i].North) = a_n;
        A.coeffRef(i, cells[i].SouthIndex) = -a_s;
    }

    // south
    for (int i = SIZE - NX + 1; i < SIZE - 1; ++i)
    {
        A.coeffRef(i, i) = a_w + a_e + a_n;
        A.coeffRef(i, cells[i].WestIndex) = -a_w;
        A.coeffRef(i, cells[i].EastIndex) = -a_e;
        A.coeffRef(i, cells[i].NorthIndex) = -a_n;
        // A.coeffRef(i, cells[i].South) = a_s;
    }

    // west
    for (int i = NX; i < SIZE - NX; i += NX)
    {
        A.coeffRef(i, i) = a_e + a_n + a_s;
        // A.coeffRef(i, cells[i].West) = a_w;
        A.coeffRef(i, cells[i].EastIndex) = -a_e;
        A.coeffRef(i, cells[i].NorthIndex) = -a_n;
        A.coeffRef(i, cells[i].SouthIndex) = -a_s;
    }

    // east
    for (int i = 2 * NX - 1; i < SIZE - NX; i += NX)
    {
        A.coeffRef(i, i) = a_w + a_n + a_s;
        A.coeffRef(i, cells[i].WestIndex) = -a_w;
        // A.coeffRef(i, cells[i].East) = a_e;
        A.coeffRef(i, cells[i].NorthIndex) = -a_n;
        A.coeffRef(i, cells[i].SouthIndex) = -a_s;
    }

    // углы
    {
        int i = 0;
        A.coeffRef(i, i) = a_e + a_s;
        A.coeffRef(i, i + 1) = -a_e;
        A.coeffRef(i, i + NX) = -a_s;
    }

    {
        int i = NX - 1;
        A.coeffRef(i, i) = a_w + a_s;
        A.coeffRef(i, i - 1) = -a_w;
        A.coeffRef(i, i + NX) = -a_s;
    }

    {
        int i = SIZE - NX;
        A.coeffRef(i, i) = a_e + a_n;
        A.coeffRef(i, i + 1) = -a_e;
        A.coeffRef(i, i - NX) = -a_n;
    }

    {
        int i = SIZE - 1;
        A.coeffRef(i, i) = a_w + a_n;
        A.coeffRef(i, i - 1) = -a_w;
        A.coeffRef(i, i - NX) = -a_n;
    }

    A.finalize();
    string filename = "A.txt";
    ofstream newfile(filename);

    if (newfile.is_open())
    {
        newfile << A;
    }
}

void Iteration(vector<Cell> &cells, vector<Cell> &cells_next, Eigen::SparseMatrix<double> &A)
{
    int t = 0;
    bool eps_flag = false;
    Eigen::VectorXd u_prev(SIZE);

    for (int i = 0; i < SIZE; ++i)
    {
        u_prev.coeffRef(i) = 0;
    }
    Eigen::VectorXd u(SIZE);
    while (eps_flag != true)
    {
        // cout << t;
        while (true)
        {
            // граничные условия на верхней грани
            for (int i = 0; i < NX; ++i)
            {
                Cell &cell = cells[i];
                cells_next[i].U_e = 1;
                cells_next[i].U_w = 1;
                cells_next[i].V_n = 0;
                cells_next[i].V_s = cell.V_s - HT * (VGradV(cell, cells) - MU * LaplasV(cell, cells) + GradPy(cell, cells));
            }

            // внутриние ячейки
            for (int i = NX; i < SIZE; ++i)
            {
                Cell &cell = cells[i];
                cells_next[i].U_e = cell.East ? cell.U_e - HT * (UGradU(cell, cells) - MU * LaplasU(cell, cells) + GradPx(cell, cells)) : 0;
                cells_next[i].V_s = cell.South ? cell.V_s - HT * (VGradV(cell, cells) - MU * LaplasV(cell, cells) + GradPy(cell, cells)) : 0;
                cells_next[i].U_w = cell.West ? cells_next[cells_next[i].WestIndex].U_e : 0;
                cells_next[i].V_n = cells_next[cells_next[i].NorthIndex].V_s;
            }

            Eigen::VectorXd b(SIZE);
            for (int i = 0; i < SIZE; ++i)
            {
                Cell &cell_prev = cells[i];
                Cell &cell = cells_next[i];
                b.coeffRef(i) = -((cell.U_e - cell.U_w) * HX +
                                  (cell.V_n - cell.V_s) * HY);
            }

            // условие выхода из итерационного процесса
            if (b.lpNorm<2>() < 1e-08)
            {

                // условие стабилизации системы
                for (int i = 0; i < SIZE; ++i)
                {
                    cells_next[i].P = cells[i].P;
                    u.coeffRef(i) = (cells_next[i].U_e + cells_next[i].U_w) / 2;
                }

                eps_flag = (u - u_prev).lpNorm<2>() < 1e-06;
                u_prev = u;

                // запись среза
                if (t % SLISE == 0 or eps_flag == true)
                {
                    ofstream newfile("./results/u_" + to_string(t) + ".txt");

                    if (newfile.is_open())
                    {
                        for (int i = 0; i < SIZE; ++i)
                        {

                            newfile << (cells_next[i].U_e + cells_next[i].U_w) / 2 << ' ';
                        }
                    }

                    newfile.close();

                    newfile.open("./results/v_" + to_string(t) + ".txt");

                    if (newfile.is_open())
                    {
                        for (int i = 0; i < SIZE; ++i)
                        {

                            newfile << (cells_next[i].V_n + cells_next[i].V_s) / 2 << ' ';
                        }
                    }

                    newfile.close();

                    newfile.open("./results/p_" + to_string(t) + ".txt");

                    if (newfile.is_open())
                    {
                        for (int i = 0; i < SIZE; ++i)
                        {

                            newfile << (cells_next[i].P) << ' ';
                        }
                    }

                    newfile.close();
                }

                // переход на следующий шаг
                t++;
                vector<Cell> temp = cells;
                cells = cells_next;
                cells_next = temp;
                break;
            }

            // решеине СЛАУ

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(A);
            solver.factorize(A);

            if (solver.info() != Eigen::Success)
            {
                cout << "ERROR" << endl;
            }
            Eigen::VectorXd p_add = solver.solve(b);

            for (int i = 0; i < SIZE; ++i)
            {
                Cell &cell = cells[i];
                cell.P += ALPHA_P * p_add[i];
            }
        }
    }
}

int main(int, char **)
{
    // установка параметров сетки
    Cell::SetNx(NX);
    Cell::SetNy(NY);
    Cell::SetCellCount(NX * NY);

    // массив ячеек (сетка) в 2 момента времени
    vector<Cell> cells;
    vector<Cell> cells_next;

    // матрица коэффициентов поправки к давлениям
    Eigen::SparseMatrix<double> A(SIZE, SIZE);

    for (int i = 0; i < SIZE; ++i)
    {
        cells_next.push_back(Cell());
    }

    Cell::ResetIndexator();

    for (int i = 0; i < SIZE; ++i)
    {
        cells.push_back(Cell());
    }

    // заполнение матрицы
    InitMattrix(A, cells);

    // запуск основного алгоритма
    Iteration(cells, cells_next, A);

    string filename = "values.txt";
    ofstream newfile("./params.txt");

    // файл параметров
    if (newfile.is_open())
    {
        newfile << L << " " << NX << " " << NY << " " << RE << " " << HT;
    }

    // std::cout << "Hello, world!\n";
}
