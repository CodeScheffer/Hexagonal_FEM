# Hexagonal_FEM
Análise linear com elemento Hexagonal 2D de 6 nós
(Manual para Windows)

Caso o rodar.bat não funcionar:

gfortran -c Precisao.f90
gfortran -c Dados.f90
gfortran -c Material.f90
gfortran -c Elemento.f90
gfortran -c Elemento_tri.f90
gfortran -c Corpo.f90
gfortran -c Corpotri.f90
gfortran -c Contorno.f90
gfortran -c Contornotri.f90
gfortran -c Tensao_hex.f90
gfortran -c Tensaotri.f90
gfortran -c Global.f90
gfortran -c Apoios.f90
gfortran -c Main.f90

gfortran -o Hexagonal Apoios.o Contorno.o Contornotri.o Corpo.o Corpotri.o Dados.o Elemento.o Elemento_tri.o Global.o Main.o Material.o Precisao.o Tensao_hex.o Tensaotri.o

Para executar o programa:

.\Hexagonal

Vai gerar um arquivo "resultados.pos" na pasta, alterar a extensão para ".txt", rodar o arquivo "POSPROCES.py", o qual vai formatar as tensões e deformação do arquivo txt para o formato no qual o arquivo plot consegue ler.

vai gerar um arquivo "por_gp_formatado.txt", copiar as informações dele, excluir no arquivo "reultados.txt" todo os dados de tensões, deformações no pontos de Gauss ( excluir esta linha " # Tensoes e Deformacoes nos Pontos de Gauss" e as linhas abaixo) e colar os valores do txt "por_gp_formatado.txt"

Rodar "plot.py" para visualização de resultados.
