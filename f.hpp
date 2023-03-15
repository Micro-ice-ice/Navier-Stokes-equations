#pragma once
#include "cell.hpp"
#include "vars.hpp"
#include <vector>

using namespace std;

double UGradU(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return cell.U_e * (cell_e.U_e - cell_w.U_e) / 2 / HX + (cell.V_n + cell.V_s + cell_e.V_n + cell_e.V_s) / 4 * (cell_n.U_e - cell_s.U_e) / 2 / HY;
}

double LaplasU(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return (cell_e.U_e - 2 * cell.U_e + cell_w.U_e) / HX / HX + (cell_n.U_e - 2 * cell.U_e + cell_s.U_e) / HY / HY;
}

double GradPx(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return (cell_e.P - cell.P) / HX;
}

double VGradV(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return (cell.U_w + cell.U_e + cell_s.U_w + cell_s.U_e) / 4 * (cell_e.V_s - cell_w.V_s) / 2 / HX + cell.V_s * (cell_n.V_s - cell_s.V_s) / 2 / HY;
}

double LaplasV(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return (cell_e.V_s - 2 * cell.V_s + cell_w.V_s) / HX / HX + (cell_n.V_s - 2 * cell.V_s + cell_s.V_s) / HY / HY;
}

double GradPy(Cell &cell, vector<Cell> &cells)
{
    Cell cell_virtual = Cell(true);
    Cell &cell_w = cell.West ? cells[cell.WestIndex] : cell_virtual;
    Cell &cell_e = cell.East ? cells[cell.EastIndex] : cell_virtual;
    Cell &cell_n = cell.North ? cells[cell.NorthIndex] : cell_virtual;
    Cell &cell_s = cell.South ? cells[cell.SouthIndex] : cell_virtual;

    return (cell.P - cell_s.P) / HY;
}