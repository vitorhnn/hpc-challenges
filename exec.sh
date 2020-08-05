#!/bin/bash
#Execute aqui pelo script.
#Parâmetros
# ./mdf 100.0 1.0 1.0 1.0 100.0 0
#        |     |   |   |    |   |-> sava log (1 salva/0 não salva)
#        |     |   |   |    |-----> Temperatura da borda
#        |     |   |   |----------> Tamanho em metros do eixo Z
#        |     |   |--------------> Tamanho em metros do eixo Y
#        |     |------------------> Tamanho em metros do eixo X
#        |------------------------> tempo em segundos
time ./mdf 10.0 42.0 42.0 42.0 100.0 1
