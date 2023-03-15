#pragma once

const double L = 1; // длина квадрата

int NX = 40; // количество ячеек по оси X

int NY = 40; // количество ячеек по оси Y

double HX = L / (double)NX; // шаг по оси X

double HY = L / (double)NY; // шаг по оси Y

int SIZE = NX * NY; // количество ячеек в сетке

double RE = 1; // число Ренольца

double MU = 1 / RE;

double HT = 0.00001; // шаг по Времени

double ALPHA = 1; // параметр релаксации для скорости

double ALPHA_P = 0.8; // параметр релаксации для давления

int SLISE = 100000;