# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 17:20:10 2020

@author: 40742
"""

# Modules to import models of
# Reaction-Diffusion pattern formation

from models import GrayScott
from models import Giraffe
from models import GrayScott_Pearson
from models import Leopard

# Modules to import all numerical solvers
from solvers import Crank_Nicolson
from solvers import Crank_NicolsonG
from solvers import Euler
from solvers import Euler1

# Module to check performance
from performance import MyTimer


def sample1_GrayScott_Euler():
    u = Euler(GrayScott(0.022, 0.061), 256, 10000, 1)
    return u.configPlot()


def sample2_GrayScott_Euler():
    u = Euler(GrayScott(0.0660, 0.0630), 256, 10000, 1)
    return u.configPlot()


def sample_GrayScott_CN():
    u = Crank_Nicolson(GrayScott(0.022, 0.061), 245, 100, 1)
    return u.configPlot()


def CN_Giraffe():
    u = Crank_NicolsonG(Giraffe(), 120, 60000, 25)
    return u.configPlot()


def Euler_Giraffe():
    u = Euler1(Giraffe(), 120, 2500, 1 / 4)
    return u.configPlot()


def CN_Leopard():
    u = Crank_NicolsonG(Leopard(), 150, 100000, 50)
    return u.configPlot()


def Euler_Leopard():
    u = Euler1(Leopard(), 150, 2500, 1 / 4)
    return u.configPlot()


def main():
    Euler_Leopard()

    with MyTimer():
        Euler1(Giraffe(), 120, 0, 1 / 4)


if __name__ == "__main__":
    main()