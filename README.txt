Existem 4 scripts nesta pasta:

1) writer.py - funções auxiliares que escrevem arquivos de entrada para o OpenFOAM

2) misc.py - funções para testes

3) main.py - funções que geram o contorno do aerofólio e avaliam os coeficientes

4) airfoil.py - (RODE ESTE AQUI) onde você define os parâmetros e chama todas as outras

A função objetivo lê eval/feval.in e escreve os coeficientes em eval/feval.out

Para rodar o script, digite "python3 airfoil.py"

Para rodar o algencan, digite "gfortran -O3 air.f90 -L$algencan-3.1.1/lib -lalgencan -o airtest" e depois "./airtest"
