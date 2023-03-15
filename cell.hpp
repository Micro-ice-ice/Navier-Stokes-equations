#pragma once

class Cell
{

private:
    static int Indexator;

    static int Nx;

    static int Ny;

    static int CellCount;

public:
    double P = 0;

    double U_w = 0; // Скорость по X слева

    double U_e = 0; // Скорость по X справа

    double V_n = 0; // Скоростб по Y сверху

    double V_s = 0; // Скоростб по Y снизу

    bool West = false; // наличие левого соседа (запад)

    bool East = false; // наличие праваого соседа (восток)

    bool North = false; // наличие верхнего соседа (север)

    bool South = false; // наличие нижнего соседа (юг)

    int Index; // индекс ячейки

    int WestIndex; // индекс левого соседа (запад)

    int EastIndex; // индекс правого соседа (восток)

    int NorthIndex; // индекс верхнего соседа (север)

    int SouthIndex; // индекс нижнего соседа (юг)

    static void ResetIndexator()
    {

        Indexator = 0;
    }

    static void SetNx(int nx) // Размерность сетки по оси X
    {

        Nx = nx;
    }

    static void SetNy(int ny) // Размерность сетки по оси Y
    {

        Ny = ny;
    }

    static void SetCellCount(int count) // Размерность сетки (количество ячеек)
    {

        CellCount = count;
    }

    Cell(/* args */);
    Cell(bool virt);
    ~Cell();
};

Cell::Cell()
{

    Index = Indexator;
    Indexator++;

    if (Index % Nx != 0)
    {

        WestIndex = Index - 1;
        West = true;
    }

    if ((Index + 1) % Nx != 0)
    {

        EastIndex = Index + 1;
        East = true;
    }

    if (Index - Nx >= 0)
    {
        NorthIndex = Index - Nx;
        North = true;
    }
    else
    {
        U_e = 1;
        U_w = 1;
    }

    if (Index + Nx < CellCount)
    {

        SouthIndex = Index + Nx;
        South = true;
    }
}

Cell::Cell(bool virt)
{
}

Cell::~Cell()
{
}

int Cell::Indexator = 0;
int Cell::Nx = 0;
int Cell::Ny = 0;
int Cell::CellCount = 0;